"""Observed technical call metrics and explicit evidence-dependent review states."""

from __future__ import annotations

from typing import Any


def assay_metrics(rows: list[dict[str, str]]) -> dict[str, Any]:
    samples = [row for row in rows if row["control_type"] == "sample"]
    called = [row for row in samples if row["call"] != "NA"]
    comparable = [row for row in called if row["reference_genotype"] != "NA"]
    matched = sum(row["call"] == row["reference_genotype"] for row in comparable)
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in samples:
        grouped.setdefault(row["sample_id"], []).append(row)
    replicate_groups = [group for group in grouped.values() if len(group) > 1]
    comparable_replicates = sum(
        sum(row["call"] != "NA" for row in group) >= 2 for group in replicate_groups
    )
    mismatches = sum(
        len({row["call"] for row in group if row["call"] != "NA"}) > 1 for group in replicate_groups
    )
    failed_controls = sum(
        (row["control_type"] == "positive" and row["call"] != row["expected_genotype"])
        or (row["control_type"] == "negative" and row["call"] != "NA")
        for row in rows
    )
    plates = {(row["batch_id"], row["plate_id"]) for row in rows}
    incomplete_plates = 0
    for batch, plate in plates:
        controls = [row for row in rows if row["batch_id"] == batch and row["plate_id"] == plate]
        expected = {
            row["expected_genotype"] for row in controls if row["control_type"] == "positive"
        }
        if expected != {"0", "1", "2"} or not any(
            row["control_type"] == "negative" for row in controls
        ):
            incomplete_plates += 1
    complete_samples = sum(all(row["call"] != "NA" for row in group) for group in grouped.values())
    return {
        "n_sample_attempts": len(samples),
        "n_called": len(called),
        "n_no_calls": len(samples) - len(called),
        "call_rate": len(called) / len(samples) if samples else None,
        "n_comparable_calls": len(comparable),
        "n_concordant": matched,
        "concordance": matched / len(comparable) if comparable else None,
        "n_samples": len(grouped),
        "n_comparable_samples": len({row["sample_id"] for row in comparable}),
        "n_complete_samples": complete_samples,
        "complete_sample_call_rate": complete_samples / len(grouped) if grouped else None,
        "n_replicated_samples": len(replicate_groups),
        "n_comparable_replicated_samples": comparable_replicates,
        "n_replicate_groups_with_no_call": sum(
            any(row["call"] == "NA" for row in group) for group in replicate_groups
        ),
        "n_replicate_disagreements": mismatches,
        "n_controls": len(rows) - len(samples),
        "n_failed_controls": failed_controls,
        "n_plates": len(plates),
        "n_incomplete_control_plates": incomplete_plates,
    }


def review_state(
    metrics: dict[str, Any], meta: dict[str, Any] | None, previous: dict[str, str] | None
) -> dict[str, Any]:
    reasons = []
    state = "computational_candidate"
    if meta is None or metrics["n_sample_attempts"] == 0:
        return {
            "state": state,
            "simulated_state": "",
            "acceptance": "not_assessed",
            "reasons": "no_sample_measurements",
        }
    if not meta["analyst_approved"] or any(
        meta[key] in ("unknown", "not_applicable")
        for key in ("analyst_review_reference", "reviewer_id")
    ):
        reasons.append("analyst_review_not_approved")
    if meta["measurement_conditions"] in ("unknown", "not_applicable"):
        reasons.append("measurement_conditions_unverified")
    for metric, threshold in (
        ("n_samples", "minimum_samples"),
        ("n_comparable_samples", "minimum_comparable_samples"),
        ("n_comparable_replicated_samples", "minimum_replicated_samples"),
        ("call_rate", "minimum_call_rate"),
        ("complete_sample_call_rate", "minimum_call_rate"),
        ("concordance", "minimum_concordance"),
    ):
        if metrics[metric] is None or metrics[metric] < meta[threshold]:
            reasons.append("insufficient_" + metric)
    for key in ("n_replicate_disagreements", "n_failed_controls", "n_incomplete_control_plates"):
        if metrics[key]:
            reasons.append(key)
    if meta["reference_method"] in ("unknown", "not_applicable"):
        reasons.append("reference_method_unverified")
    if not reasons and meta["requested_state"] != state:
        state = "assay_validated"
        if meta["requested_state"] == "population_validated":
            previous_state = previous.get("state") if previous else None
            if (
                meta["data_scope"] == "synthetic"
                and previous
                and previous.get("data_scope") == "synthetic"
            ):
                previous_state = previous.get("simulated_state")
            if previous_state not in ("assay_validated", "population_validated"):
                reasons.append("prior_assay_validation_required")
            if (
                not meta["population_scope_confirmed"]
                or meta["population_evidence_reference"] in ("unknown", "not_applicable")
                or meta["target_population"] in ("unknown", "not_applicable")
            ):
                reasons.append("population_scope_evidence_required")
            if previous and previous.get("validation_dataset_id") == meta["validation_dataset_id"]:
                reasons.append("distinct_population_validation_dataset_required")
            if not reasons:
                state = "population_validated"
    result = {
        "state": state,
        "simulated_state": "",
        "acceptance": "accepted" if not reasons else "rejected_requested_state",
        "reasons": ";".join(reasons) or "criteria_met_within_declared_scope",
    }
    if meta["requested_state"] == "computational_candidate":
        result.update(
            acceptance="promotion_not_requested",
            reasons=";".join(reasons) or "promotion_not_requested",
        )
    if meta["data_scope"] == "synthetic":
        result.update(
            simulated_state=state,
            state="computational_candidate",
            reasons=result["reasons"] + ";synthetic_only_not_operational_validation",
        )
    return result
