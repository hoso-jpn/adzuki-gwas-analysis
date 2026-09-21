"""Strict laboratory return contract bound to a computational primer version."""

from __future__ import annotations

import math
import tomllib
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.reference import required_text
from adzuki_gwas_analysis.tables import read_tsv

IDENTITY_FIELDS = ("candidate_id", "design_id", "primer_version", "assembly_id", "ref", "alt")
RESULT_FIELDS = (
    *IDENTITY_FIELDS,
    "allele_encoding",
    "conditions_id",
    "plate_id",
    "batch_id",
    "sample_id",
    "replicate_id",
    "control_type",
    "call",
    "no_call_reason",
    "expected_genotype",
    "reference_genotype",
)
STATES = ("computational_candidate", "assay_validated", "population_validated")
TEXT_METADATA = (
    "validation_dataset_id",
    "conditions_id",
    "measurement_conditions",
    "target_population",
    "reference_method",
    "analyst_review_reference",
    "reviewer_id",
    "population_evidence_reference",
)
INT_METADATA = ("minimum_samples", "minimum_comparable_samples", "minimum_replicated_samples")
RATE_METADATA = ("minimum_call_rate", "minimum_concordance")


def load_assay_metadata(path: Path) -> dict[str, Any]:
    meta = tomllib.loads(path.read_text(encoding="utf-8"))
    expected = set(TEXT_METADATA + INT_METADATA + RATE_METADATA) | {
        "schema_version",
        "data_scope",
        "requested_state",
        "analyst_approved",
        "population_scope_confirmed",
    }
    if (
        set(meta) != expected
        or type(meta.get("schema_version")) is not int
        or meta["schema_version"] != 1
    ):
        raise ValueError("incomplete/unknown assay metadata contract")
    for key in TEXT_METADATA:
        required_text(meta, key)
    if (
        meta["data_scope"] not in ("synthetic", "customer", "public")
        or meta["requested_state"] not in STATES
    ):
        raise ValueError("unknown assay scope or requested state")
    for key in ("analyst_approved", "population_scope_confirmed"):
        if type(meta[key]) is not bool:
            raise ValueError("review flags must be explicit booleans")
    for key in INT_METADATA:
        if type(meta[key]) is not int or meta[key] < (
            1 if key == "minimum_replicated_samples" else 2
        ):
            raise ValueError(
                "sample thresholds require at least two samples and one replicate group"
            )
    for key in RATE_METADATA:
        value = meta[key]
        if (
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(value)
            or not 0 < value <= 1
        ):
            raise ValueError("rate thresholds must be finite fractions in (0,1]")
    return meta


def load_assay_results(
    path: Path, designs: dict[str, dict[str, str]], meta: dict[str, Any], *, max_rows: int
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()
    sample_definitions: dict[tuple[str, str], tuple[str, ...]] = {}
    for row in read_tsv(path, RESULT_FIELDS):
        if len(rows) >= max_rows or set(row) != set(RESULT_FIELDS):
            raise ValueError("result limit exceeded or unknown result columns")
        design = designs.get(row["design_id"])
        if design is None or any(row[key] != design[key] for key in IDENTITY_FIELDS):
            raise ValueError("unknown design, old primer version, assembly or allele mismatch")
        if (
            row["allele_encoding"] != "ALT_0_1_2_forward"
            or row["conditions_id"] != meta["conditions_id"]
        ):
            raise ValueError("allele direction or measurement conditions mismatch")
        if any(not row[key] for key in ("sample_id", "replicate_id", "plate_id", "batch_id")):
            raise ValueError("plate/batch/sample/replicate IDs must be explicit")
        if any(
            row[key] not in ("0", "1", "2", "NA")
            for key in ("call", "expected_genotype", "reference_genotype")
        ):
            raise ValueError("calls must be literal 0/1/2/NA ALT dosages")
        control = row["control_type"]
        if control not in ("sample", "positive", "negative"):
            raise ValueError("unknown control type")
        if control == "positive" and (
            row["expected_genotype"] == "NA" or row["reference_genotype"] != "NA"
        ):
            raise ValueError(
                "positive controls require an expected genotype, separate from sample truth"
            )
        if control == "negative" and (
            row["expected_genotype"] != "NA" or row["reference_genotype"] != "NA"
        ):
            raise ValueError("negative controls must expect no call")
        if control == "sample" and row["expected_genotype"] != "NA":
            raise ValueError(
                "samples use reference_genotype; expected_genotype is reserved for controls"
            )
        if row["call"] == "NA" and control != "negative" and not row["no_call_reason"]:
            raise ValueError("a no-call reason is required")
        if row["call"] != "NA" and row["no_call_reason"]:
            raise ValueError("a called genotype cannot carry a no-call reason")
        identity = tuple(
            row[key] for key in ("design_id", "batch_id", "plate_id", "sample_id", "replicate_id")
        )
        if identity in seen:
            raise ValueError("duplicate assay result ID")
        seen.add(identity)
        sample_key = (row["candidate_id"], row["sample_id"])
        definition = tuple(
            row[key] for key in ("control_type", "expected_genotype", "reference_genotype")
        )
        if sample_key in sample_definitions and sample_definitions[sample_key] != definition:
            raise ValueError("sample truth/control definition changes between replicates")
        sample_definitions[sample_key] = definition
        rows.append(row)
    return rows
