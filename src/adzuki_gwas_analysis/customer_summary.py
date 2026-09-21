"""Explicit plant summary-statistics adapter with log-space multiple testing."""

from __future__ import annotations

import argparse
import hashlib
import math
import tomllib
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    json_bytes,
    output_transaction,
    write_json,
)
from adzuki_gwas_analysis.reference import (
    ReferenceBundle,
    load_reference_bundle,
    required_text,
    reverse_complement,
)
from adzuki_gwas_analysis.tables import read_tsv

REQUIRED_METADATA = (
    "dataset_id",
    "cohort_id",
    "analysis_id",
    "reference",
    "assembly_id",
    "species",
    "trait",
    "trait_unit",
    "trait_coding",
    "test",
    "effect_type",
    "effect_scale",
    "standard_error_scale",
    "effect_strand",
    "target_effect_allele",
    "sample_size_definition",
)
REQUIRED_COLUMNS = (
    "chr",
    "pos",
    "ref",
    "alt",
    "effect_allele",
    "other_allele",
    "effect",
    "standard_error",
    "pvalue",
    "sample_size",
)


def load_metadata(path: Path, bundle: ReferenceBundle) -> dict[str, Any]:
    document = tomllib.loads(path.read_text(encoding="utf-8"))
    if type(document.get("schema_version")) is not int or document["schema_version"] != 1:
        raise ValueError("unsupported customer summary schema")
    for key in REQUIRED_METADATA:
        required_text(document, key)
    if (
        set(document)
        - set(REQUIRED_METADATA)
        - {"schema_version", "columns", "invalid_row_policy", "chi_square_df", "chi_square_basis"}
    ):
        raise ValueError("unknown customer summary metadata field")
    bundle.check_dataset(document["dataset_id"], document["reference"], document["assembly_id"])
    if document["species"] != bundle.species:
        raise ValueError("species metadata does not match reference bundle")
    choices = {
        "effect_type": ("beta", "odds_ratio"),
        "effect_strand": ("forward", "unknown"),
        "target_effect_allele": ("reported", "ALT"),
        "invalid_row_policy": ("error", "exclude"),
    }
    document.setdefault("invalid_row_policy", "error")
    for key, allowed in choices.items():
        if document[key] not in allowed:
            raise ValueError(f"unsupported {key}")
    if document["effect_type"] == "odds_ratio" and (
        document["standard_error_scale"] != "log_odds" or document["effect_scale"] != "log_odds"
    ):
        raise ValueError("OR inputs require SE on the log-odds scale")
    if (
        document["effect_type"] == "beta"
        and document["standard_error_scale"] != document["effect_scale"]
    ):
        raise ValueError("beta and SE must use the same scale")
    mapping = document.get("columns")
    if not isinstance(mapping, dict) or not set(REQUIRED_COLUMNS) <= set(mapping):
        raise ValueError("columns must explicitly map every required field")
    allowed_columns = set(REQUIRED_COLUMNS) | {"effect_allele_frequency", "neg_log10_pvalue"}
    if set(mapping) - allowed_columns or any(
        not isinstance(x, str) or not x for x in mapping.values()
    ):
        raise ValueError("invalid/unknown column mapping")
    df = document.get("chi_square_df")
    if df is not None:
        if type(df) is not int or df <= 0:
            raise ValueError("chi_square_df must be a positive integer")
        required_text(document, "chi_square_basis")
        if document["test"] == "unknown":
            raise ValueError("an unknown test cannot establish chi-square degrees of freedom")
    return document


def _number(value: str, name: str) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise ValueError(f"invalid {name}") from exc
    if not math.isfinite(result):
        raise ValueError(f"nonfinite {name}")
    return result


def neg_log10_pvalue(value: str, log_value: str | None = None) -> tuple[float, str]:
    try:
        pvalue = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError("invalid pvalue") from exc
    if not pvalue.is_finite() or not 0 <= pvalue <= 1:
        raise ValueError("pvalue outside [0,1]")
    if pvalue == 0:
        if log_value is None:
            raise ValueError("zero pvalue requires an original neg_log10_pvalue")
        result = _number(log_value, "neg_log10_pvalue")
        if result <= 0:
            raise ValueError("rounded zero pvalue needs positive neg_log10_pvalue")
        return result, "original_neg_log10_for_rounded_zero"
    result = float(-pvalue.log10())
    if log_value is not None and not math.isclose(
        result, _number(log_value, "neg_log10_pvalue"), rel_tol=1e-10, abs_tol=1e-10
    ):
        raise ValueError("pvalue and neg_log10_pvalue disagree")
    if not math.isfinite(result):
        raise ValueError("log pvalue exceeds the supported numeric range")
    return result, "decimal_pvalue"


def normalize_row(
    row: dict[str, str],
    meta: dict[str, Any],
    bundle: ReferenceBundle,
    *,
    row_number: int,
    source_sha256: str,
) -> dict[str, Any]:
    for key in (
        "dataset_id",
        "cohort_id",
        "analysis_id",
        "reference",
        "assembly_id",
        "trait",
        "test",
    ):
        if key in row and row[key] != meta[key]:
            raise ValueError("row identity disagrees with declared " + key)
    values = {role: row[column] for role, column in meta["columns"].items()}
    chrom, pos = values["chr"], int(values["pos"])
    ref, alt = values["ref"], values["alt"]
    bundle.check_snp(chrom, pos, ref, alt)
    ea, oa = values["effect_allele"], values["other_allele"]
    if len(ea) != 1 or len(oa) != 1 or ea == oa or any(c not in "ACGT" for c in ea + oa):
        raise ValueError("invalid effect/other allele")
    if meta["effect_strand"] == "forward" and {ea, oa} != {ref, alt}:
        raise ValueError("forward effect alleles disagree with REF/ALT")
    if meta["effect_strand"] == "unknown" and {ea, oa} not in (
        {ref, alt},
        {reverse_complement(ref), reverse_complement(alt)},
    ):
        raise ValueError("effect alleles cannot be reconciled with REF/ALT on either strand")
    n = int(values["sample_size"])
    if n <= 0:
        raise ValueError("sample_size must be a positive integer")
    effect = _number(values["effect"], "effect")
    se = _number(values["standard_error"], "standard_error")
    if se <= 0:
        raise ValueError("standard_error must be positive")
    if meta["effect_type"] == "odds_ratio":
        if effect <= 0:
            raise ValueError("odds ratio must be positive")
        beta = math.log(effect)
    else:
        beta = effect
    frequency: float | None = None
    if "effect_allele_frequency" in values:
        frequency = _number(values["effect_allele_frequency"], "effect_allele_frequency")
        if not 0 <= frequency <= 1:
            raise ValueError("effect allele frequency outside [0,1]")
    logp, p_source = neg_log10_pvalue(values["pvalue"], values.get("neg_log10_pvalue"))
    orientation = "unresolved_strand" if meta["effect_strand"] == "unknown" else "as_reported"
    normalized_ea, normalized_oa = ea, oa
    if orientation != "unresolved_strand" and meta["target_effect_allele"] == "ALT":
        orientation = "aligned_to_ALT"
        if ea == ref:
            beta = -beta
            frequency = None if frequency is None else 1 - frequency
        normalized_ea, normalized_oa = alt, ref
    identity = {"assembly_id": bundle.assembly_id, "chr": chrom, "pos": pos, "ref": ref, "alt": alt}
    return {
        "schema_version": 1,
        **{
            k: meta[k]
            for k in (
                "dataset_id",
                "cohort_id",
                "analysis_id",
                "reference",
                "assembly_id",
                "trait",
                "trait_unit",
                "trait_coding",
                "test",
                "effect_scale",
            )
        },
        **identity,
        "candidate_id": "variant-" + hashlib.sha256(json_bytes(identity)).hexdigest()[:20],
        "effect_allele": normalized_ea,
        "other_allele": normalized_oa,
        "effect_orientation": orientation,
        "source_effect_allele": ea,
        "source_other_allele": oa,
        "source_effect": effect,
        "beta": beta,
        "se": se,
        "effect_allele_frequency": frequency,
        "source_pvalue": values["pvalue"],
        "neg_log10_pvalue": logp,
        "pvalue_source": p_source,
        "sample_size": n,
        "source_row": row_number,
        "source_sha256": source_sha256,
        "variant_identity_source": "explicit_REF_ALT_checked_against_bundle",
        "favorable_allele": "not_inferred",
    }


def adjust_and_cluster(
    rows: list[dict[str, Any]],
    *,
    alpha: float,
    fdr_level: float,
    clustering_distance: int,
    contig_order: list[str],
) -> list[dict[str, Any]]:
    if not rows or not 0 < alpha < 1 or not 0 < fdr_level < 1 or clustering_distance < 0:
        raise ValueError("no valid rows or invalid correction/clustering settings")
    count = len(rows)
    by_p = sorted(rows, key=lambda row: (-row["neg_log10_pvalue"], row["candidate_id"]))
    running = 0.0
    for index in range(count - 1, -1, -1):
        row = by_p[index]
        running = max(running, row["neg_log10_pvalue"] - math.log10(count / (index + 1)))
        row["neg_log10_p_bh"] = running
        row["neg_log10_p_bonferroni"] = max(0.0, row["neg_log10_pvalue"] - math.log10(count))
        row["bonferroni_significant"] = row["neg_log10_p_bonferroni"] >= -math.log10(alpha)
        row["bh_significant"] = running >= -math.log10(fdr_level)
    candidates = [
        row.copy() for row in rows if row["bonferroni_significant"] or row["bh_significant"]
    ]
    order = {contig: index for index, contig in enumerate(contig_order)}
    ordered = sorted(
        candidates, key=lambda row: (order[row["chr"]], row["pos"], row["candidate_id"])
    )
    last_chrom, last_pos, signal = None, 0, 0
    for row in ordered:
        if row["chr"] != last_chrom or row["pos"] - last_pos > clustering_distance:
            signal += 1
        row["signal_id"] = f"signal-{signal:05d}"
        last_chrom, last_pos = row["chr"], row["pos"]
    ordered.sort(
        key=lambda row: (
            -row["neg_log10_pvalue"],
            order[row["chr"]],
            row["pos"],
            row["candidate_id"],
        )
    )
    for rank, row in enumerate(ordered, 1):
        row["candidate_rank"] = rank
        row["priority_tier"] = "Bonferroni" if row["bonferroni_significant"] else "BH_only"
    return ordered


def run_customer_summary(
    *,
    summary: Path,
    metadata_path: Path,
    bundle_path: Path,
    output_dir: Path,
    clustering_distance: int,
    alpha: float = 0.05,
    fdr_level: float = 0.05,
    max_rows: int = 1_000_000,
) -> None:
    from adzuki_gwas_analysis.summary_reporting import write_summary_artifacts

    if max_rows <= 0:
        raise ValueError("max_rows must be positive")
    bundle = load_reference_bundle(bundle_path)
    meta = load_metadata(metadata_path, bundle)
    inputs = {
        "summary": summary,
        "metadata": metadata_path,
        "reference_bundle": bundle_path,
        "reference_fasta": bundle.fasta,
    }
    initial, environment = input_checksums(inputs), generation_environment()
    rows: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    identities: set[str] = set()
    for row_number, raw in enumerate(read_tsv(summary, meta["columns"].values()), 2):
        if row_number - 1 > max_rows:
            raise ValueError("summary exceeds declared max_rows; no partial analysis is published")
        try:
            row = normalize_row(
                raw, meta, bundle, row_number=row_number, source_sha256=initial["summary"]
            )
        except (ValueError, KeyError) as exc:
            if meta["invalid_row_policy"] == "error":
                raise ValueError(f"summary row {row_number}: {exc}") from exc
            excluded.append({"source_row": row_number, "reason": str(exc)})
            continue
        if row["candidate_id"] in identities:
            raise ValueError("duplicate variant identity within the analysis family")
        identities.add(row["candidate_id"])
        rows.append(row)
    contigs = [contig.name for contig in bundle.contigs]
    candidates = adjust_and_cluster(
        rows,
        alpha=alpha,
        fdr_level=fdr_level,
        clustering_distance=clustering_distance,
        contig_order=contigs,
    )
    with output_transaction(output_dir) as stage:
        diagnostics = write_summary_artifacts(
            stage,
            rows=rows,
            candidates=candidates,
            excluded=excluded,
            meta=meta,
            contig_order=contigs,
            alpha=alpha,
            fdr_level=fdr_level,
        )
        write_json(stage / "statistical_diagnostics.json", diagnostics)
        finish_provenance(
            stage,
            operation="customer_summary",
            inputs=inputs,
            initial=initial,
            environment=environment,
            families=[diagnostics["family"]],
            parameters={
                "metadata": meta,
                "alpha": alpha,
                "fdr_level": fdr_level,
                "clustering_distance": clustering_distance,
                "max_rows": max_rows,
                "n_input": len(rows) + len(excluded),
                "n_tests": len(rows),
                "n_excluded": len(excluded),
            },
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("summary", "metadata-path", "bundle-path", "output-dir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--clustering-distance", type=int, required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--fdr-level", type=float, default=0.05)
    parser.add_argument("--max-rows", type=int, default=1_000_000)
    args = parser.parse_args(argv)
    try:
        run_customer_summary(**vars(args))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"customer summary failed: {type(exc).__name__}: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
