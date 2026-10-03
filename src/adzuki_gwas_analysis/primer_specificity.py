"""Offline, declared-rule amplicon search and optional pinned Primer3 thermodynamics."""

from __future__ import annotations

import argparse
import importlib
import json
import math
from importlib import metadata
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.arms import PRIMER_FIELDS
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    output_transaction,
    validate_provenance,
    write_json,
)
from adzuki_gwas_analysis.reference import (
    ReferenceBundle,
    load_reference_bundle,
    reverse_complement,
)
from adzuki_gwas_analysis.tables import read_tsv, write_tsv

ENGINE = "ungapped-terminal-seed-amplicon-v1"
PASS = "passed_within_declared_reference_and_rules"


def load_settings(path: Path) -> dict[str, Any]:
    settings: dict[str, Any] = json.loads(path.read_text())
    expected = {
        "schema_version",
        "max_mismatches",
        "exact_3prime_bases",
        "min_product_bp",
        "max_product_bp",
        "max_hits_per_primer",
        "max_candidate_windows",
        "max_products_per_reaction",
        "max_reference_bases",
        "max_contig_bases",
        "max_designs",
        "thermodynamics",
    }
    if (
        not isinstance(settings, dict)
        or set(settings) != expected
        or type(settings["schema_version"]) is not int
        or settings["schema_version"] != 1
    ):
        raise ValueError("incomplete/unknown specificity settings")
    for field in expected - {"schema_version", "thermodynamics"}:
        if type(settings[field]) is not int or settings[field] < (
            0 if field == "max_mismatches" else 1
        ):
            raise ValueError("search limits require explicit nonnegative/positive integers")
    if settings["max_mismatches"] > 3 or not 1 <= settings["exact_3prime_bases"] <= 20:
        raise ValueError("mismatch/terminal-seed rule is outside supported range")
    if settings["min_product_bp"] > settings["max_product_bp"]:
        raise ValueError("invalid product length range")
    thermo = settings["thermodynamics"]
    fields = {
        "enabled",
        "mv_conc",
        "dv_conc",
        "dntp_conc",
        "dna_conc",
        "temp_c",
        "max_loop",
        "max_structure_tm_C",
    }
    if not isinstance(thermo, dict) or set(thermo) != fields or type(thermo["enabled"]) is not bool:
        raise ValueError("explicit thermodynamic conditions are required")
    for field in fields - {"enabled", "max_loop"}:
        if type(thermo[field]) not in (int, float) or not math.isfinite(thermo[field]):
            raise ValueError("thermodynamic conditions must be finite numbers")
    if (
        any(thermo[key] < 0 for key in ("mv_conc", "dv_conc", "dntp_conc"))
        or thermo["dna_conc"] <= 0
        or type(thermo["max_loop"]) is not int
        or not 0 <= thermo["max_loop"] <= 30
        or not 0 <= thermo["temp_c"] <= 100
    ):
        raise ValueError("thermodynamic composition/temperature is outside supported range")
    return settings


def find_hits(
    sequence: str, primer: str, *, chrom: str, settings: dict[str, Any]
) -> list[dict[str, Any]]:
    """Enumerate ungapped hits with an exact primer 3-prime seed on either strand."""
    seed_length = settings["exact_3prime_bases"]
    if (
        not 2 <= len(primer) <= 60
        or len(primer) < seed_length
        or any(c not in "ACGT" for c in primer)
    ):
        raise ValueError("unsupported primer length/bases or terminal seed length")
    result = []
    candidate_windows = 0
    for strand, query in (("+", primer), ("-", reverse_complement(primer))):
        seed_offset = len(query) - seed_length if strand == "+" else 0
        seed = query[seed_offset : seed_offset + seed_length]
        cursor = sequence.find(seed)
        while cursor >= 0:
            candidate_windows += 1
            if candidate_windows > settings["max_candidate_windows"]:
                raise ValueError("candidate window budget exceeded; specificity is unresolved")
            start = cursor - seed_offset
            end = start + len(query)
            if start >= 0 and end <= len(sequence):
                window = sequence[start:end]
                if all(base in "ACGT" for base in window):
                    mismatches = sum(a != b for a, b in zip(window, query, strict=True))
                    if mismatches <= settings["max_mismatches"]:
                        result.append(
                            dict(
                                chr=chrom,
                                start=start + 1,
                                end=end,
                                strand=strand,
                                mismatches=mismatches,
                            )
                        )
                        if len(result) > settings["max_hits_per_primer"]:
                            raise ValueError("primer hit limit exceeded; specificity is unresolved")
            cursor = sequence.find(seed, cursor + 1)
    return result


def amplicons(
    first: list[dict[str, Any]],
    second: list[dict[str, Any]],
    settings: dict[str, Any],
    *,
    first_role: str = "specific",
    second_role: str = "common",
) -> list[dict[str, Any]]:
    products = []
    seen = set()
    for a in first:
        for b in second:
            if a["chr"] != b["chr"]:
                continue
            plus, minus = (a, b) if a["strand"] == "+" else (b, a)
            if plus["strand"] != "+" or minus["strand"] != "-" or plus["end"] >= minus["start"]:
                continue
            length = minus["end"] - plus["start"] + 1
            if settings["min_product_bp"] <= length <= settings["max_product_bp"]:
                first_hit, second_hit = (plus, minus) if first_role == second_role else (a, b)
                identity = (
                    a["chr"],
                    first_hit["start"],
                    first_hit["strand"],
                    second_hit["start"],
                    second_hit["strand"],
                )
                if identity in seen:
                    continue
                seen.add(identity)
                products.append(
                    dict(
                        chr=a["chr"],
                        start=plus["start"],
                        end=minus["end"],
                        product_bp=length,
                        primer1_role=first_role,
                        primer2_role=second_role,
                        primer1_start=first_hit["start"],
                        primer1_strand=first_hit["strand"],
                        primer2_start=second_hit["start"],
                        primer2_strand=second_hit["strand"],
                        primer1_mismatches=first_hit["mismatches"],
                        primer2_mismatches=second_hit["mismatches"],
                    )
                )
                if len(products) > settings["max_products_per_reaction"]:
                    raise ValueError("amplicon budget exceeded; specificity is unresolved")
    return products


def assess_design(
    design: dict[str, str], bundle: ReferenceBundle, settings: dict[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    bundle.check_dataset(design["dataset_id"], design["reference"], design["assembly_id"])
    bundle.check_snp(design["chr"], int(design["pos"]), design["ref"], design["alt"])
    if design["strand"] not in ("+", "-"):
        raise ValueError("invalid design strand")
    for role in ("specific", "common"):
        start, end = int(design[f"{role}_start"]), int(design[f"{role}_end"])
        sequences = (
            [design[f"specific_{allele}_5to3"] for allele in ("ref", "alt")]
            if role == "specific"
            else [design["common_5to3"]]
        )
        if not 1 <= start <= end <= bundle.contig(design["chr"]).length or any(
            len(sequence) != end - start + 1 for sequence in sequences
        ):
            raise ValueError("design primer coordinates/length are inconsistent")
    specific_terminal = "specific_end" if design["strand"] == "+" else "specific_start"
    if int(design[specific_terminal]) != int(design["pos"]):
        raise ValueError("allele-specific primer 3-prime coordinate does not match SNP")
    for allele in ("ref", "alt"):
        terminal_base = (
            design[allele] if design["strand"] == "+" else reverse_complement(design[allele])
        )
        if not design[f"specific_{allele}_5to3"].endswith(terminal_base):
            raise ValueError("allele-specific primer 3-prime base does not match declared allele")
    if sum(c.length for c in bundle.contigs) > settings["max_reference_bases"] or any(
        c.length > settings["max_contig_bases"] for c in bundle.contigs
    ):
        raise ValueError("reference exceeds declared operating limit")
    expected_start = min(int(design["specific_start"]), int(design["common_start"]))
    expected_end = max(int(design["specific_end"]), int(design["common_end"]))
    if expected_end - expected_start + 1 != int(design["product_bp"]):
        raise ValueError("design amplicon interval is inconsistent")
    all_hits: list[dict[str, Any]] = []
    all_products: list[dict[str, Any]] = []
    assessments: list[dict[str, Any]] = []
    for allele in ("ref", "alt"):
        products, ambiguous = [], False
        counts = {"specific": 0, "common": 0}
        for contig in bundle.contigs:
            sequence = bundle.sequence(contig.name, 1, contig.length)
            ambiguous |= any(base not in "ACGT" for base in sequence)
            # ALT primer's intended target is tested on the explicitly declared single-SNP
            # alternate haplotype, not incorrectly required to match the reference 3' base.
            if allele == "alt" and contig.name == design["chr"]:
                offset = int(design["pos"]) - 1
                sequence = sequence[:offset] + design["alt"] + sequence[offset + 1 :]
            specific = find_hits(
                sequence, design[f"specific_{allele}_5to3"], chrom=contig.name, settings=settings
            )
            common = find_hits(
                sequence, design["common_5to3"], chrom=contig.name, settings=settings
            )
            for role, found in (("specific", specific), ("common", common)):
                counts[role] += len(found)
                if counts[role] > settings["max_hits_per_primer"]:
                    raise ValueError("genome-wide primer hit limit exceeded")
            all_hits.extend(
                {"design_id": design["design_id"], "reaction": allele, "primer_role": role, **hit}
                for role, hits in (("specific", specific), ("common", common))
                for hit in hits
            )
            products.extend(amplicons(specific, common, settings))
            products.extend(
                amplicons(
                    specific, specific, settings, first_role="specific", second_role="specific"
                )
            )
            products.extend(
                amplicons(common, common, settings, first_role="common", second_role="common")
            )
            if len(products) > settings["max_products_per_reaction"]:
                raise ValueError("genome-wide amplicon budget exceeded")
        expected = [
            p
            for p in products
            if (p["chr"], p["start"], p["end"]) == (design["chr"], expected_start, expected_end)
            and (p["primer1_role"], p["primer2_role"]) == ("specific", "common")
            and p["primer1_start"] == int(design["specific_start"])
            and p["primer2_start"] == int(design["common_start"])
            and p["primer1_strand"] == design["strand"]
            and p["primer2_strand"] == ("-" if design["strand"] == "+" else "+")
        ]
        off_target = len(products) - len(expected)
        if len(expected) != 1 or off_target:
            status, reason = "failed", "intended_amplicon_missing_or_multiple_products"
        elif ambiguous:
            status, reason = "ambiguous", "reference_contains_unresolved_bases"
        else:
            status, reason = PASS, "one_intended_product_within_declared_ungapped_rules"
        assessments.append(
            dict(
                design_id=design["design_id"],
                primer_version=design["primer_version"],
                assembly_id=design["assembly_id"],
                data_scope=bundle.data_scope,
                reference_search_scope="supplied_reference_contigs_only",
                reference_completeness="not_assessed",
                reaction=allele,
                specificity=status,
                reason=reason,
                n_products=len(products),
                n_off_target_products=off_target,
                cohort_binding_check=design["cohort_binding_check"],
                assay_validation="not_assessed",
            )
        )
        all_products.extend(
            {"design_id": design["design_id"], "reaction": allele, **p} for p in products
        )
    return all_hits, all_products, assessments


def assess_thermodynamics(design: dict[str, str], settings: dict[str, Any]) -> dict[str, Any]:
    configuration = settings["thermodynamics"]
    if not configuration["enabled"]:
        return {
            "design_id": design["design_id"],
            "status": "not_checked",
            "reason": "not_requested",
        }
    try:
        primer3 = importlib.import_module("primer3")
        version = metadata.version("primer3-py")
    except (ImportError, metadata.PackageNotFoundError) as exc:
        raise ValueError(
            "thermodynamics requested; install the locked thermodynamics extra"
        ) from exc
    if version != "2.3.1":
        raise ValueError("thermodynamics requires locked primer3-py 2.3.1")
    kwargs = {
        key: configuration[key]
        for key in ("mv_conc", "dv_conc", "dntp_conc", "dna_conc", "temp_c", "max_loop")
    }
    oligos = {
        "ref": design["specific_ref_5to3"],
        "alt": design["specific_alt_5to3"],
        "common": design["common_5to3"],
    }
    calculations = {
        f"{role}_{kind}": getattr(primer3, "calc_" + kind)(sequence, **kwargs)
        for role, sequence in oligos.items()
        for kind in ("hairpin", "homodimer")
    }
    for role in ("ref", "alt"):
        calculations[f"{role}_common_heterodimer"] = primer3.calc_heterodimer(
            oligos[role], oligos["common"], **kwargs
        )
    results = {}
    for name, result in calculations.items():
        values = {key: float(getattr(result, key)) for key in ("tm", "dg", "dh", "ds")}
        if any(not math.isfinite(value) for value in values.values()):
            raise ValueError("nonfinite thermodynamic result")
        results[name] = {**values, "structure_found": bool(result.structure_found)}
    failed = any(
        r["structure_found"] and r["tm"] > configuration["max_structure_tm_C"]
        for r in results.values()
    )
    return {
        "design_id": design["design_id"],
        "status": "failed" if failed else "passed_within_declared_conditions",
        "engine": "primer3-py",
        "version": version,
        "license": "GPL-2.0-or-later",
        "conditions": configuration,
        "metrics": results,
        "assay_validation": "not_assessed",
    }


def run_specificity(
    *, design_dir: Path, context_dir: Path, bundle_path: Path, settings_path: Path, output_dir: Path
) -> None:
    record = validate_provenance(design_dir, required=True)
    if record is None or record["operation"] != "arms_design":
        raise ValueError("specificity requires verified ARMS design provenance")
    settings, bundle = load_settings(settings_path), load_reference_bundle(bundle_path)
    context = validate_provenance(context_dir, required=True)
    if (
        context is None
        or context["operation"] != "candidate_context"
        or record["parameters"].get("context_run_id") != context["run_id"]
        or context["parameters"]["reference"]["asset_sha256"]["fasta"] != bundle.checksums["fasta"]
        or context["parameters"]["reference"]["assembly_id"] != bundle.assembly_id
    ):
        raise ValueError("design context/reference identity mismatch")
    inputs = {
        "design/" + p.relative_to(design_dir).as_posix(): p
        for p in design_dir.rglob("*")
        if p.is_file()
    }
    inputs.update(
        settings=settings_path, reference_bundle=bundle_path, reference_fasta=bundle.fasta
    )
    inputs.update(
        {
            "context/" + p.relative_to(context_dir).as_posix(): p
            for p in context_dir.rglob("*")
            if p.is_file()
        }
    )
    initial, environment = input_checksums(inputs), generation_environment()
    designs = list(read_tsv(design_dir / "primer_candidates.tsv", PRIMER_FIELDS))
    if len(designs) > settings["max_designs"] or len({d["design_id"] for d in designs}) != len(
        designs
    ):
        raise ValueError("design limit exceeded or duplicate design ID")
    hits, products, assessments, thermodynamics = [], [], [], []
    for design in designs:
        found, predicted, assessed = assess_design(design, bundle, settings)
        hits.extend(found)
        products.extend(predicted)
        assessments.extend(assessed)
        thermodynamics.append(assess_thermodynamics(design, settings))
    with output_transaction(output_dir) as stage:
        write_tsv(
            stage / "primer_hits.tsv",
            ("design_id", "reaction", "primer_role", "chr", "start", "end", "strand", "mismatches"),
            hits,
        )
        write_tsv(
            stage / "amplicons.tsv",
            (
                "design_id",
                "reaction",
                "chr",
                "start",
                "end",
                "product_bp",
                "primer1_role",
                "primer2_role",
                "primer1_start",
                "primer1_strand",
                "primer2_start",
                "primer2_strand",
                "primer1_mismatches",
                "primer2_mismatches",
            ),
            products,
        )
        write_tsv(
            stage / "specificity_review.tsv",
            (
                "design_id",
                "primer_version",
                "assembly_id",
                "data_scope",
                "reference_search_scope",
                "reference_completeness",
                "reaction",
                "specificity",
                "reason",
                "n_products",
                "n_off_target_products",
                "cohort_binding_check",
                "assay_validation",
            ),
            assessments,
        )
        write_json(stage / "thermodynamics.json", {"schema_version": 1, "designs": thermodynamics})
        (stage / "specificity_report.md").write_text(
            "# Computational primer confirmation\n\n"
            "Search is ungapped, bounded by mismatch count and an exact 3-prime terminal seed. "
            "Pass means one intended amplicon under those rules in the supplied reference. "
            "Reference completeness is not assessed: a supplied sequence subset or short "
            "synthetic reference is not evidence of whole-genome specificity. "
            "Unresolved reference bases prevent pass; absent assemblies, structural variation, "
            "mismatched-terminal priming and unprovided cohort variants "
            "remain outside this check.\n\n"
            "REF and ALT reactions are evaluated separately; the ALT target uses the declared "
            "single-SNP alternate haplotype. This is not evidence of allele discrimination in PCR. "
            "Thermodynamics uses explicit chemical conditions when requested. "
            "Neither computational check promotes assay or trait-utility validation. "
            "Give specificity_review.tsv and thermodynamics.json to the laboratory with the "
            "original design/version; do not replace assay acceptance with these results.\n",
            encoding="utf-8",
        )
        finish_provenance(
            stage,
            operation="primer_specificity",
            inputs=inputs,
            initial=initial,
            parameters={
                "engine": ENGINE,
                "settings": settings,
                "reference": bundle.metadata(),
                "design_run_id": record["run_id"],
            },
            environment=environment,
            families=record["families"],
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("design-dir", "context-dir", "bundle-path", "settings-path", "output-dir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    try:
        run_specificity(**vars(parser.parse_args(argv)))
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(1, f"primer confirmation failed: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
