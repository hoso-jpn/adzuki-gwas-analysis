"""Review independent quantitative-trait evidence separately from genotyping accuracy.

This is an evidence consumer, not a trial-analysis engine. Numerical effects and
intervals must come from the declared, reviewed external analysis.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.assay_contract import IDENTITY_FIELDS
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    output_transaction,
    validate_provenance,
    write_json,
)
from adzuki_gwas_analysis.tables import read_tsv, write_tsv

ASSAY_MAPPING = {
    "computational_candidate": "not_validated",
    "assay_validated": "assay_conditions_validated",
    "population_validated": "genotyping_validated_in_declared_population",
}
TEXT_FIELDS = (
    *IDENTITY_FIELDS,
    "trait",
    "trait_unit",
    "target_population",
    "environment",
    "discovery_dataset_id",
    "validation_dataset_id",
    "analysis_plan_reference",
    "analysis_engine",
    "analysis_version",
    "analysis_model",
    "analysis_run_sha256",
    "validation_data_sha256",
    "discovery_data_sha256",
    "reviewer_id",
    "review_reference",
    "baseline_strategy",
    "relatedness_review_reference",
)


def _read_units(path: Path, *, validation: bool) -> dict[str, dict[str, str]]:
    required: tuple[str, ...] = ("sample_id", "family_id")
    if validation:
        required += ("population", "environment", "trait", "unit", "value")
    units: dict[str, dict[str, str]] = {}
    for row in read_tsv(path, required):
        if set(row) != set(required) or any(not row[key] for key in required):
            raise ValueError("incomplete/unknown biological-unit contract")
        if row["sample_id"] in units:
            raise ValueError(
                "duplicate biological sample; technical replicates are not independent"
            )
        if len(units) >= 100000:
            raise ValueError("trait evidence unit limit exceeded")
        if validation and not math.isfinite(float(row["value"])):
            raise ValueError("nonfinite validation phenotype")
        units[row["sample_id"]] = row
    if not units:
        raise ValueError("empty biological-unit evidence")
    return units


def _effect(value: Any, unit: str) -> tuple[float, float, float]:
    if not isinstance(value, dict) or set(value) != {"estimate", "lower95", "upper95", "unit"}:
        raise ValueError("effect requires estimate, 95% interval and unit")
    if value["unit"] != unit:
        raise ValueError("effect/phenotype unit mismatch")
    numbers = [value[key] for key in ("estimate", "lower95", "upper95")]
    if any(type(v) not in (int, float) or not math.isfinite(v) for v in numbers):
        raise ValueError("effect interval must contain finite numbers")
    estimate, lower, upper = map(float, numbers)
    if not lower <= estimate <= upper or lower >= upper:
        raise ValueError("effect interval is inconsistent")
    return estimate, lower, upper


def assess_evidence(
    meta: dict[str, Any],
    design: dict[str, str],
    discovery: dict[str, dict[str, str]],
    validation: dict[str, dict[str, str]],
    hashes: dict[str, str],
) -> dict[str, Any]:
    expected = set(TEXT_FIELDS) | {
        "schema_version",
        "data_scope",
        "desired_direction",
        "effect_allele",
        "analyst_approved",
        "independence_basis",
        "minimum_samples",
        "minimum_effect",
        "minimum_utility",
        "association_effect",
        "selection_contrast",
    }
    if (
        set(meta) != expected
        or type(meta["schema_version"]) is not int
        or meta["schema_version"] != 1
    ):
        raise ValueError("incomplete/unknown trait evidence contract")
    for field in TEXT_FIELDS:
        value = meta[field]
        if not isinstance(value, str) or not value or any(c in value for c in "\t\r\n"):
            raise ValueError("trait evidence metadata must be explicit single-line text")
    if any(meta[field] != design[field] for field in IDENTITY_FIELDS):
        raise ValueError("trait evidence uses a different candidate/design/assembly/allele version")
    if meta["effect_allele"] != design["alt"]:
        raise ValueError("trait evidence effect must be explicitly normalized to genomic ALT")
    if meta["data_scope"] not in ("synthetic", "customer", "public"):
        raise ValueError("unknown evidence scope")
    if type(meta["analyst_approved"]) is not bool or meta["desired_direction"] not in (
        "increase",
        "decrease",
    ):
        raise ValueError("explicit reviewer approval and trait direction are required")
    if meta["independence_basis"] != "independent_families":
        raise ValueError(
            "v1 accepts independent-family validation only; related trials need separate review"
        )
    if type(meta["minimum_samples"]) is not int or meta["minimum_samples"] < 4:
        raise ValueError("at least four independent biological samples must be prespecified")
    for field in ("minimum_effect", "minimum_utility"):
        if (
            type(meta[field]) not in (int, float)
            or not math.isfinite(meta[field])
            or meta[field] < 0
        ):
            raise ValueError("acceptance effects must be finite nonnegative thresholds")
    for field in ("analysis_run_sha256", "validation_data_sha256", "discovery_data_sha256"):
        if not re.fullmatch("[0-9a-f]{64}", meta[field]):
            raise ValueError("analysis and input evidence require SHA-256 identities")
    if any(meta[f"{role}_data_sha256"] != hashes[role] for role in ("validation", "discovery")):
        raise ValueError("trait evidence does not match the supplied input data")
    if meta["analysis_run_sha256"] != hashes["analysis_run"]:
        raise ValueError("external analysis run record checksum mismatch")
    association = _effect(meta["association_effect"], meta["trait_unit"])
    utility = _effect(meta["selection_contrast"], meta["trait_unit"])
    reasons = []
    if meta["discovery_dataset_id"] == meta["validation_dataset_id"] or set(discovery) & set(
        validation
    ):
        reasons.append("discovery_validation_reuse")
    discovered_families = {row["family_id"] for row in discovery.values()}
    validation_families = [row["family_id"] for row in validation.values()]
    if discovered_families & set(validation_families) or len(set(validation_families)) != len(
        validation
    ):
        reasons.append("independent_family_scope_not_met")
    if len(validation) < meta["minimum_samples"]:
        reasons.append("insufficient_biological_samples")
    if any(
        row[key] != meta[field]
        for row in validation.values()
        for key, field in (
            ("population", "target_population"),
            ("environment", "environment"),
            ("trait", "trait"),
            ("unit", "trait_unit"),
        )
    ):
        reasons.append("validation_trait_population_environment_mismatch")
    if meta["target_population"] != design["target_population"]:
        reasons.append("assay_population_mismatch")
    if not meta["analyst_approved"] or any(
        meta[key] in ("unknown", "not_applicable")
        for key in (
            "reviewer_id",
            "review_reference",
            "analysis_plan_reference",
            "relatedness_review_reference",
            "baseline_strategy",
        )
    ):
        reasons.append("review_or_prespecified_plan_missing")
    sign = 1 if meta["desired_direction"] == "increase" else -1
    effect_lower = association[1] if sign == 1 else -association[2]
    utility_lower = utility[1] if sign == 1 else -utility[2]
    effect_supported = effect_lower > meta["minimum_effect"]
    replication = "supported" if effect_supported and not reasons else "not_supported"
    if not effect_supported:
        reasons.append("association_effect_criterion_not_met")
    if utility_lower <= meta["minimum_utility"]:
        reasons.append("selection_utility_criterion_not_met")
    if design["state"] != "population_validated":
        reasons.append("population_genotyping_validation_required")
    proposed = "supported_within_declared_scope" if not reasons else "not_supported"
    synthetic = meta["data_scope"] == "synthetic" or design["data_scope"] == "synthetic"
    return {
        "genotyping_validation": ASSAY_MAPPING[design["state"]],
        "association_replication": "simulation_only" if synthetic else replication,
        "trait_utility": "not_assessed" if synthetic else proposed,
        "simulated_trait_utility": proposed if synthetic else "",
        "reasons": ";".join(reasons)
        or (
            "synthetic_not_operational_validation"
            if synthetic
            else "criteria_met_within_declared_scope"
        ),
        "n_biological_samples": len(validation),
        "n_families": len(set(validation_families)),
        "association_effect": meta["association_effect"],
        "selection_contrast": meta["selection_contrast"],
    }


def run_trait_review(
    *,
    assay_dir: Path,
    output_dir: Path,
    evidence_path: Path | None = None,
    discovery_path: Path | None = None,
    validation_path: Path | None = None,
    analysis_run_path: Path | None = None,
) -> None:
    supplied = (evidence_path, discovery_path, validation_path, analysis_run_path)
    if any(p is not None for p in supplied) and not all(p is not None for p in supplied):
        raise ValueError("trait evidence, discovery and validation must be supplied together")
    record = validate_provenance(assay_dir, required=True)
    if record is None or record["operation"] != "assay_review":
        raise ValueError("trait review requires a verified assay_review bundle")
    inputs = {
        "assay/" + p.relative_to(assay_dir).as_posix(): p
        for p in assay_dir.rglob("*")
        if p.is_file()
    }
    for role, path in zip(
        ("evidence", "discovery", "validation", "analysis_run"), supplied, strict=True
    ):
        if path is not None:
            inputs[role] = path
    initial, environment = input_checksums(inputs), generation_environment()
    designs = list(
        read_tsv(
            assay_dir / "marker_states.tsv",
            (*IDENTITY_FIELDS, "state", "data_scope", "target_population"),
        )
    )
    if len({d["design_id"] for d in designs}) != len(designs) or any(
        d["state"] not in ASSAY_MAPPING for d in designs
    ):
        raise ValueError("invalid assay states or duplicate design identity")
    evidence: dict[str, Any] | None = None
    assessed = None
    if evidence_path is not None and discovery_path is not None and validation_path is not None:
        evidence = json.loads(evidence_path.read_text())
        if not isinstance(evidence, dict):
            raise ValueError("trait evidence must be an object")
        matches = [d for d in designs if d["design_id"] == evidence.get("design_id")]
        if len(matches) != 1:
            raise ValueError("trait evidence design is absent from assay bundle")
        assessed = assess_evidence(
            evidence,
            matches[0],
            _read_units(discovery_path, validation=False),
            _read_units(validation_path, validation=True),
            initial,
        )
    states = []
    for design in designs:
        state = {
            **{key: design[key] for key in IDENTITY_FIELDS},
            "genotyping_validation": ASSAY_MAPPING[design["state"]],
            "association_replication": "not_assessed",
            "trait_utility": "not_assessed",
            "simulated_trait_utility": "",
            "reasons": "no_trait_evidence",
            "target_population": design["target_population"],
            "trait": "",
            "environment": "",
        }
        if (
            evidence is not None
            and design["design_id"] == evidence["design_id"]
            and assessed is not None
        ):
            state.update(
                {
                    key: assessed[key]
                    for key in (
                        "genotyping_validation",
                        "association_replication",
                        "trait_utility",
                        "simulated_trait_utility",
                        "reasons",
                    )
                }
            )
            state.update(trait=evidence["trait"], environment=evidence["environment"])
        states.append(state)
    with output_transaction(output_dir) as stage:
        fields = (
            *IDENTITY_FIELDS,
            "genotyping_validation",
            "association_replication",
            "trait_utility",
            "simulated_trait_utility",
            "reasons",
            "target_population",
            "trait",
            "environment",
        )
        write_tsv(stage / "evidence_axes.tsv", fields, states)
        write_tsv(
            stage / "trait_supported_panel.tsv",
            fields,
            [s for s in states if s["trait_utility"] == "supported_within_declared_scope"],
        )
        write_json(
            stage / "trait_review.json",
            {
                "schema_version": 1,
                "evidence": evidence,
                "assessment": assessed,
                "assay_run_id": record["run_id"],
            },
        )
        (stage / "trait_review_report.md").write_text(
            "# Genotyping and trait-selection evidence\n\n"
            "Genotyping accuracy is independent from trait-selection utility. "
            "Legacy population_validated means genotyping validation "
            "in the declared population.\n\n"
            "Trait endpoints and intervals are consumed from the declared external analysis; "
            "this tool does not recompute or authenticate that analysis. Review its original run, "
            "phenotype preparation, prespecified selection rule and comparator. "
            "Support is limited to the named trait, population, environment and analysis. "
            "It does not establish causality or general yield gain.\n\n"
            "See evidence_axes.tsv, trait_review.json and the checksum-bound source evidence. "
            "Synthetic evidence never promotes operational trait utility.\n",
            encoding="utf-8",
        )
        finish_provenance(
            stage,
            operation="trait_utility_review",
            inputs=inputs,
            initial=initial,
            parameters={"assay_run_id": record["run_id"], "evidence": evidence},
            environment=environment,
            families=[],
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in (
        "assay-dir",
        "output-dir",
        "evidence-path",
        "discovery-path",
        "validation-path",
        "analysis-run-path",
    ):
        parser.add_argument("--" + flag, type=Path, required=flag in ("assay-dir", "output-dir"))
    try:
        run_trait_review(**vars(parser.parse_args(argv)))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"trait evidence review failed: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
