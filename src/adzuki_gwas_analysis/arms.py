"""Computational two-reaction ARMS candidates with explicit screening assumptions."""

from __future__ import annotations

import argparse
import hashlib
import math
import tomllib
from collections import Counter
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any
from urllib.parse import unquote

from adzuki_gwas_analysis.primer_metrics import complementary_run, primer_metrics
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    json_bytes,
    output_transaction,
    validate_provenance,
)
from adzuki_gwas_analysis.reference import reverse_complement
from adzuki_gwas_analysis.tables import read_tsv, write_tsv

ENGINE = "adzuki-arms-screen-v1"


@dataclass(frozen=True)
class PrimerSettings:
    min_length: int = 18
    max_length: int = 25
    min_product_bp: int = 75
    max_product_bp: int = 200
    min_tm_C: float = 45
    max_tm_C: float = 70
    target_tm_C: float = 60
    max_tm_difference_C: float = 5
    min_gc_fraction: float = 0.35
    max_gc_fraction: float = 0.65
    max_homopolymer: int = 4
    max_self_match: int = 8
    max_pair_match: int = 8
    sodium_mM: float = 50
    strand_concentration_nM: float = 25
    alternatives: int = 3

    def validate(self) -> None:
        integer_keys = (
            "min_length",
            "max_length",
            "min_product_bp",
            "max_product_bp",
            "max_homopolymer",
            "max_self_match",
            "max_pair_match",
            "alternatives",
        )
        if any(
            type(getattr(self, key)) is not int or getattr(self, key) <= 0 for key in integer_keys
        ):
            raise ValueError("primer count/length settings must be positive integers")
        if not all(
            isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
            for v in asdict(self).values()
        ):
            raise ValueError("primer settings must be finite numeric values")
        if not (
            2 <= self.min_length <= self.max_length <= 60
            and self.min_length < self.min_product_bp <= self.max_product_bp <= 2000
            and self.min_tm_C <= self.target_tm_C <= self.max_tm_C
            and self.max_tm_difference_C >= 0
            and 0 <= self.min_gc_fraction <= self.max_gc_fraction <= 1
            and 0 < self.sodium_mM <= 1000
            and 0 < self.strand_concentration_nM <= 10000
            and self.alternatives <= 100
        ):
            raise ValueError("primer settings are inconsistent or outside supported bounds")


def _fasta(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    name: str | None = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = unquote(line[1:].split()[0])
            if name in result:
                raise ValueError("duplicate candidate FASTA ID")
            result[name] = ""
        elif name is None:
            raise ValueError("missing FASTA header")
        else:
            result[name] += line
    return result


def _metric_reasons(metrics: dict[str, float], settings: PrimerSettings) -> list[str]:
    reasons = []
    if not settings.min_tm_C <= metrics["tm_C"] <= settings.max_tm_C:
        reasons.append("tm_out_of_range")
    if not settings.min_gc_fraction <= metrics["gc_fraction"] <= settings.max_gc_fraction:
        reasons.append("gc_out_of_range")
    if metrics["max_homopolymer"] > settings.max_homopolymer:
        reasons.append("homopolymer")
    if metrics["self_complementary_run"] > settings.max_self_match:
        reasons.append("self_complementarity_screen")
    return reasons


def _genomic_interval(candidate: dict[str, str], start: int, end: int) -> tuple[int, int]:
    if candidate["strand"] == "+":
        return int(candidate["interval_start"]) + start, int(candidate["interval_start"]) + end - 1
    return int(candidate["interval_end"]) - end + 1, int(candidate["interval_end"]) - start


def design_candidate(
    candidate: dict[str, str],
    sequence: str,
    neighbors: list[dict[str, str]],
    settings: PrimerSettings,
    *,
    cohort_available: bool,
    neighbor_window: int,
) -> tuple[list[dict[str, Any]], Counter[str]]:
    reasons: Counter[str] = Counter()
    designs: list[dict[str, Any]] = []
    offset, pos = int(candidate["snp_offset_0based"]), int(candidate["pos"])
    ref, alt = candidate["ref"], candidate["alt"]
    if len(ref) != 1 or len(alt) != 1 or ref == alt or any(c not in "ACGT" for c in ref + alt):
        return [], Counter({"not_biallelic_snp": 1})
    oriented_ref, oriented_alt = (
        (ref, alt)
        if candidate["strand"] == "+"
        else (reverse_complement(ref), reverse_complement(alt))
    )
    if not 0 <= offset < len(sequence) or sequence[offset] != oriented_ref:
        raise ValueError("context allele/offset does not match sequence")
    if any(c not in "ACGT" for c in sequence):
        return [], Counter({"ambiguous_flanking_bases": 1})
    if offset < settings.min_length - 1 or len(sequence) - offset < settings.min_length:
        return [], Counter({"insufficient_flanks": 1})
    cache: dict[str, dict[str, float]] = {}

    def metrics(oligo: str) -> dict[str, float]:
        if oligo not in cache:
            cache[oligo] = primer_metrics(
                oligo,
                sodium_mM=settings.sodium_mM,
                strand_concentration_nM=settings.strand_concentration_nM,
            )
        return cache[oligo]

    for length in range(settings.min_length, settings.max_length + 1):
        left = offset - length + 1
        if left < 0:
            continue
        specific_ref = sequence[left : offset + 1]
        specific_alt = specific_ref[:-1] + oriented_alt
        ref_metrics, alt_metrics = metrics(specific_ref), metrics(specific_alt)
        failures = _metric_reasons(ref_metrics, settings) + _metric_reasons(alt_metrics, settings)
        if failures:
            reasons.update(failures)
            continue
        for product in range(settings.min_product_bp, settings.max_product_bp + 1):
            right_end = left + product
            if right_end > len(sequence):
                break
            for common_length in range(settings.min_length, settings.max_length + 1):
                right_start = right_end - common_length
                if right_start <= offset:
                    continue
                common = reverse_complement(sequence[right_start:right_end])
                common_metrics = metrics(common)
                failures = _metric_reasons(common_metrics, settings)
                tm_values = [m["tm_C"] for m in (ref_metrics, alt_metrics, common_metrics)]
                difference = max(tm_values) - min(tm_values)
                if difference > settings.max_tm_difference_C:
                    failures.append("tm_difference")
                if failures:
                    reasons.update(failures)
                    continue
                pair_match = max(
                    complementary_run(specific_ref, common), complementary_run(specific_alt, common)
                )
                if pair_match > settings.max_pair_match:
                    reasons["primer_dimer_screen"] += 1
                    continue
                intervals = (
                    _genomic_interval(candidate, left, offset + 1),
                    _genomic_interval(candidate, right_start, right_end),
                )
                known_overlaps = sum(
                    any(
                        int(n["pos"]) <= end and int(n["pos"]) + len(n["ref"]) - 1 >= start
                        for start, end in intervals
                    )
                    for n in neighbors
                    if n["source"] == "cohort_vcf"
                )
                if known_overlaps:
                    reasons["cohort_variant_in_binding_site"] += 1
                    continue
                coverage = "unavailable"
                if cohort_available:
                    coverage = (
                        "checked_within_context_window"
                        if all(
                            pos - neighbor_window <= start <= end <= pos + neighbor_window
                            for start, end in intervals
                        )
                        else "partial_window_unverified"
                    )
                score = sum(abs(tm - settings.target_tm_C) for tm in tm_values) + difference
                design_id = (
                    "design-"
                    + hashlib.sha256(
                        json_bytes(
                            {
                                "candidate_id": candidate["candidate_id"],
                                "assembly": candidate["assembly_id"],
                                "sequences": [specific_ref, specific_alt, common],
                                "intervals": intervals,
                                "settings": asdict(settings),
                                "engine": ENGINE,
                            }
                        )
                    ).hexdigest()[:20]
                )
                row: dict[str, Any] = {
                    **{
                        key: candidate[key]
                        for key in (
                            "candidate_id",
                            "dataset_id",
                            "reference",
                            "assembly_id",
                            "chr",
                            "pos",
                            "ref",
                            "alt",
                            "strand",
                        )
                    },
                    "design_id": design_id,
                    "primer_version": design_id,
                    "score": score,
                    "specific_ref_5to3": specific_ref,
                    "specific_alt_5to3": specific_alt,
                    "common_5to3": common,
                    "product_bp": product,
                    "pair_complementary_run": pair_match,
                    "specific_start": intervals[0][0],
                    "specific_end": intervals[0][1],
                    "common_start": intervals[1][0],
                    "common_end": intervals[1][1],
                    "cohort_binding_check": coverage,
                    "specificity": "unverified",
                    "validation_state": "computational_candidate",
                    "engine": ENGINE,
                }
                for prefix, values in (
                    ("ref", ref_metrics),
                    ("alt", alt_metrics),
                    ("common", common_metrics),
                ):
                    row.update({prefix + "_" + key: value for key, value in values.items()})
                designs.append(row)
                designs.sort(
                    key=lambda item: (item["score"], item["product_bp"], item["design_id"])
                )
                del designs[settings.alternatives :]
    designs.sort(key=lambda row: (row["score"], row["product_bp"], row["design_id"]))
    for rank, row in enumerate(designs[: settings.alternatives], 1):
        row["rank"] = rank
        row["selection_reason"] = (
            "passed configured screens; rank by Tm score, product size, stable ID"
        )
    if not designs and not reasons:
        reasons["insufficient_product_interval"] += 1
    return designs[: settings.alternatives], reasons


PRIMER_FIELDS = (
    "candidate_id",
    "dataset_id",
    "reference",
    "assembly_id",
    "chr",
    "pos",
    "ref",
    "alt",
    "strand",
    "design_id",
    "primer_version",
    "rank",
    "score",
    "specific_ref_5to3",
    "specific_alt_5to3",
    "common_5to3",
    "product_bp",
    "specific_start",
    "specific_end",
    "common_start",
    "common_end",
    "pair_complementary_run",
    "cohort_binding_check",
    "specificity",
    "validation_state",
    "engine",
    "selection_reason",
) + tuple(
    prefix + "_" + metric
    for prefix in ("ref", "alt", "common")
    for metric in ("length", "gc_fraction", "tm_C", "max_homopolymer", "self_complementary_run")
)


def run_arms(*, context_dir: Path, config: Path, output_dir: Path) -> None:
    record = validate_provenance(context_dir, required=True)
    if record is None or record["operation"] != "candidate_context":
        raise ValueError("ARMS design requires verified candidate_context provenance")
    configuration = tomllib.loads(config.read_text())
    if configuration.get("schema_version") != 1 or not isinstance(
        configuration.get("primer"), dict
    ):
        raise ValueError("primer config requires schema_version=1 and [primer]")
    overrides = configuration["primer"]
    if set(overrides) - {field.name for field in fields(PrimerSettings)}:
        raise ValueError("unknown primer setting")
    settings = PrimerSettings(**overrides)
    settings.validate()
    inputs = {"settings": config, **{p.name: p for p in context_dir.iterdir() if p.is_file()}}
    initial, environment = input_checksums(inputs), generation_environment()
    sequences = _fasta(context_dir / "flanking_sequences.fasta")
    candidates = list(
        read_tsv(
            context_dir / "candidate_context.tsv",
            ("candidate_id", "ref", "alt", "snp_offset_0based"),
        )
    )
    if set(sequences) != {c["candidate_id"] for c in candidates} or len(sequences) != len(
        candidates
    ):
        raise ValueError("candidate/FASTA ID accounting mismatch")
    neighbors = list(read_tsv(context_dir / "neighboring_variants.tsv"))
    designs: list[dict[str, Any]] = []
    markers: list[dict[str, Any]] = []
    for candidate in candidates:
        sequence = sequences[candidate["candidate_id"]]
        if hashlib.sha256(sequence.encode()).hexdigest() != candidate["sequence_sha256"]:
            raise ValueError("candidate sequence checksum mismatch")
        selected, reasons = design_candidate(
            candidate,
            sequence,
            [n for n in neighbors if n["candidate_id"] == candidate["candidate_id"]],
            settings,
            cohort_available="cohort_vcf" in record["inputs"],
            neighbor_window=int(record["parameters"]["neighbor_window_bp"]),
        )
        designs.extend(selected)
        markers.append(
            {
                **candidate,
                "design_count": len(selected),
                "validation_state": "computational_candidate" if selected else "not_designed",
                "rejection_reasons": ";".join(
                    f"{key}:{value}" for key, value in sorted(reasons.items())
                ),
                "specificity": "unverified",
                "experimental_validation": "not_performed",
            }
        )
    with output_transaction(output_dir) as stage:
        write_tsv(stage / "primer_candidates.tsv", PRIMER_FIELDS, designs)
        write_tsv(
            stage / "arms_marker_candidates.tsv",
            (
                "candidate_id",
                "dataset_id",
                "reference",
                "assembly_id",
                "chr",
                "pos",
                "ref",
                "alt",
                "strand",
                "design_count",
                "validation_state",
                "rejection_reasons",
                "mask_overlap",
                "specificity",
                "experimental_validation",
            ),
            markers,
        )
        (stage / "marker_design_report.md").write_text(
            "# ARMS computational candidate report\n\n"
            f"Candidates assessed: {len(markers)}. Primer designs returned: {len(designs)}.\n\n"
            "Each design describes two separate reactions: "
            "REF-specific + common and ALT-specific + common. "
            "The discriminating base is the allele-specific primer's 3-prime terminal base. "
            "No intentional secondary mismatch is introduced.\n\n"
            "Tm: perfect-match DNA/DNA nearest-neighbor model, SantaLucia 1998, "
            "monovalent salt only. Mg, dNTP, DMSO, mismatches and PCR kinetics are not modeled. "
            "Self/dimer values are contiguous "
            "complementarity screens, not hairpin/dimer free energies.\n\n"
            "Genome-wide specificity is unverified. "
            "Cohort checks cover only the declared context window; "
            "unavailable or partial coverage is not proof that binding sites lack variants.\n\n"
            "Review alternatives and repeat/variant risks in the TSVs. "
            "Confirm specificity, optimize PCR conditions, include controls and replicates, "
            "and measure call rate/concordance in the intended "
            "population before any validated-marker claim. No experiments have been performed.\n",
            encoding="utf-8",
        )
        finish_provenance(
            stage,
            operation="arms_design",
            inputs=inputs,
            initial=initial,
            environment=environment,
            families=record["families"],
            parameters={
                "primer": asdict(settings),
                "engine": ENGINE,
                "engine_license": "MIT (project)",
                "context_run_id": record["run_id"],
                "model_reference": "https://doi.org/10.1073/pnas.95.4.1460",
            },
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("context-dir", "config", "output-dir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        run_arms(**vars(args))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"ARMS design failed: {type(exc).__name__}: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
