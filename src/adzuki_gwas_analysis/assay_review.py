"""Prepare laboratory handoff and review returned assay evidence without external actions."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.arms import PRIMER_FIELDS
from adzuki_gwas_analysis.assay_contract import (
    IDENTITY_FIELDS,
    RESULT_FIELDS,
    load_assay_metadata,
    load_assay_results,
)
from adzuki_gwas_analysis.assay_metrics import assay_metrics, review_state
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    output_transaction,
    validate_provenance,
    write_json,
)
from adzuki_gwas_analysis.tables import read_tsv, write_tsv

HISTORY_FIELDS = (
    "candidate_id",
    "prior_design_id",
    "replacement_design_id",
    "reason",
    "evidence_reference",
    "analyst_id",
)
STATE_FIELDS = (
    *IDENTITY_FIELDS,
    "data_scope",
    "validation_dataset_id",
    "conditions_id",
    "target_population",
    "state",
    "simulated_state",
    "acceptance",
    "reasons",
)


def _record_inputs(root: Path, prefix: str) -> dict[str, Path]:
    return {
        prefix + "/" + p.relative_to(root).as_posix(): p for p in root.rglob("*") if p.is_file()
    }


def run_assay_review(
    *,
    design_dir: Path,
    output_dir: Path,
    results: Path | None = None,
    metadata_path: Path | None = None,
    previous_review: Path | None = None,
    redesign_history: Path | None = None,
    max_rows: int = 100_000,
) -> None:
    if (results is None) != (metadata_path is None) or max_rows <= 0:
        raise ValueError(
            "assay results and metadata must be supplied together with a valid row limit"
        )
    record = validate_provenance(design_dir, required=True)
    if record is None or record["operation"] != "arms_design":
        raise ValueError("assay review requires verified ARMS design provenance")
    inputs = _record_inputs(design_dir, "design")
    designs: dict[str, dict[str, str]] = {}
    for row in read_tsv(design_dir / "primer_candidates.tsv", PRIMER_FIELDS):
        if row["design_id"] in designs or any(not row[key] for key in IDENTITY_FIELDS):
            raise ValueError("duplicate/incomplete primer identity")
        designs[row["design_id"]] = row
    markers = list(
        read_tsv(
            design_dir / "arms_marker_candidates.tsv",
            ("candidate_id", "design_count", "rejection_reasons"),
        )
    )
    previous_rows: dict[str, dict[str, str]] = {}
    previous_record: dict[str, Any] | None = None
    if previous_review is not None:
        previous_record = validate_provenance(previous_review, required=True)
        if previous_record is None or previous_record["operation"] != "assay_review":
            raise ValueError("previous review must be a verified assay_review bundle")
        inputs.update(_record_inputs(previous_review, "previous_review"))
        previous_rows = {
            row["design_id"]: row
            for row in read_tsv(previous_review / "marker_states.tsv", STATE_FIELDS)
        }
        for identifier, design in designs.items():
            if identifier in previous_rows and any(
                design[key] != previous_rows[identifier][key] for key in IDENTITY_FIELDS
            ):
                raise ValueError(
                    "previous evidence uses a different primer version/allele contract"
                )
    history: list[dict[str, str]] = []
    if redesign_history is not None:
        if previous_record is None:
            raise ValueError("redesign history requires the verified previous review")
        inputs["redesign_history"] = redesign_history
        edges: set[tuple[str, str]] = set()
        for row in read_tsv(redesign_history, HISTORY_FIELDS):
            if len(history) >= max_rows:
                raise ValueError("redesign history exceeds declared row limit")
            old, new = (
                previous_rows.get(row["prior_design_id"]),
                designs.get(row["replacement_design_id"]),
            )
            if (
                old is None
                or new is None
                or old["design_id"] == new["design_id"]
                or any(not row[key] for key in HISTORY_FIELDS)
            ):
                raise ValueError(
                    "redesign must link a previous design to a distinct "
                    "known replacement with evidence"
                )
            if (
                any(old[key] != new[key] for key in ("candidate_id", "assembly_id", "ref", "alt"))
                or row["candidate_id"] != new["candidate_id"]
            ):
                raise ValueError("redesign history changes candidate or allele identity")
            edge = (old["design_id"], new["design_id"])
            if edge in edges:
                raise ValueError("duplicate redesign history edge")
            edges.add(edge)
            history.append(row)
    meta = load_assay_metadata(metadata_path) if metadata_path is not None else None
    if results is not None and metadata_path is not None:
        inputs.update(results=results, metadata=metadata_path)
    initial, environment = input_checksums(inputs), generation_environment()
    rows = (
        load_assay_results(results, designs, meta, max_rows=max_rows)
        if results is not None and meta is not None
        else []
    )
    by_design: dict[str, list[dict[str, str]]] = {identifier: [] for identifier in designs}
    for row in rows:
        by_design[row["design_id"]].append(row)
    states: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    for identifier, design in designs.items():
        values = assay_metrics(by_design[identifier])
        assessment = review_state(values, meta, previous_rows.get(identifier))
        states.append(
            {
                **{key: design[key] for key in IDENTITY_FIELDS},
                **assessment,
                "data_scope": meta["data_scope"] if meta else "no_measurements",
                "validation_dataset_id": meta["validation_dataset_id"] if meta else "",
                "conditions_id": meta["conditions_id"] if meta else "",
                "target_population": meta["target_population"] if meta else "",
            }
        )
        metrics.append({"design_id": identifier, "candidate_id": design["candidate_id"], **values})
    alternatives = [
        {
            "candidate_id": state["candidate_id"],
            "design_id": state["design_id"],
            "alternative_design_id": other["design_id"],
            "alternative_rank": other["rank"],
            "reason": "available_computational_alternative_not_validated",
        }
        for state in states
        if state["state"] == "computational_candidate"
        for other in designs.values()
        if other["candidate_id"] == state["candidate_id"]
        and other["design_id"] != state["design_id"]
    ]
    with output_transaction(output_dir) as stage:
        write_tsv(
            stage / "handoff.tsv",
            (*PRIMER_FIELDS, "allele_encoding", "design_run_id"),
            [
                {**row, "allele_encoding": "ALT_0_1_2_forward", "design_run_id": record["run_id"]}
                for row in designs.values()
            ],
        )
        # A header-only template is never interpreted as an experiment or a measurement.
        write_tsv(stage / "assay_return_template.tsv", RESULT_FIELDS, [])
        write_tsv(stage / "redesign_history_template.tsv", HISTORY_FIELDS, [])
        write_tsv(stage / "assay_results.tsv", RESULT_FIELDS, rows)
        write_tsv(stage / "marker_states.tsv", STATE_FIELDS, states)
        panel_fields = (
            *IDENTITY_FIELDS,
            "specific_ref_5to3",
            "specific_alt_5to3",
            "common_5to3",
            "product_bp",
            "state",
            "conditions_id",
            "target_population",
            "validation_dataset_id",
        )
        write_tsv(
            stage / "population_validated_panel.tsv",
            panel_fields,
            [
                {**designs[row["design_id"]], **row}
                for row in states
                if row["state"] == "population_validated"
            ],
        )
        write_tsv(
            stage / "assay_metrics.tsv", ("design_id", "candidate_id", *assay_metrics([])), metrics
        )
        write_tsv(
            stage / "design_failures.tsv",
            ("candidate_id", "design_count", "rejection_reasons"),
            [row for row in markers if int(row["design_count"]) == 0],
        )
        write_tsv(
            stage / "alternative_designs.tsv",
            ("candidate_id", "design_id", "alternative_design_id", "alternative_rank", "reason"),
            alternatives,
        )
        write_tsv(stage / "redesign_history.tsv", HISTORY_FIELDS, history)
        write_json(
            stage / "assay_review.json",
            {
                "schema_version": 1,
                "design_run_id": record["run_id"],
                "previous_review_run_id": previous_record["run_id"] if previous_record else None,
                "metadata": meta,
                "n_results": len(rows),
                "n_designs": len(designs),
                "states": {
                    state: sum(row["state"] == state for row in states)
                    for state in (
                        "computational_candidate",
                        "assay_validated",
                        "population_validated",
                    )
                },
            },
        )
        report = [
            "# Assay evidence review",
            "",
            f"Primer designs: {len(designs)}. Returned rows: {len(rows)}.",
            "",
            "| Design | State | Scope | Reasons |",
            "|---|---|---|---|",
        ]
        report.extend(
            f"| {row['design_id']} | {row['state']} | {row['data_scope']} | {row['reasons']} |"
            for row in states
        )
        report.extend(
            [
                "",
                "Measurement conditions, population, review evidence, thresholds and dataset "
                "are in assay_review.json. Denominators and controls are in assay_metrics.tsv; "
                "sample returns are in assay_results.tsv. These contain private assay information.",
                "",
                "Call rate = called sample attempts / all sample attempts "
                "(technical replicates included). "
                "Complete-sample call rate requires all attempts for a sample to be called. "
                "Concordance = matching calls / called attempts with known reference genotype. "
                "No-calls and unavailable truth are excluded from the concordance denominator. "
                "Controls are excluded from both rates and checked separately. "
                "Replicate disagreement counts sample groups with two different observed calls; "
                "groups with no-calls are also counted.",
                "",
                "Every design/batch/plate needs positive controls for 0/1/2 "
                "and a negative control. Passing a small experiment, including 100% observed "
                "concordance, is not a general performance guarantee. Synthetic evidence "
                "never promotes the operational state; simulated_state is test-only.",
                "",
                "Computational analysis supplies sequences, candidates and audit records. "
                "The laboratory owns assay execution, raw-call interpretation, "
                "controls and measurement conditions; "
                "an analyst must review acceptance within the declared population and conditions. "
                "No laboratory orders, purchases, messages or experiments "
                "are performed by this command.",
                "",
            ]
        )
        (stage / "assay_review_report.md").write_text("\n".join(report), encoding="utf-8")
        finish_provenance(
            stage,
            operation="assay_review",
            inputs=inputs,
            initial=initial,
            environment=environment,
            families=record["families"],
            parameters={
                "metadata": meta,
                "max_rows": max_rows,
                "design_run_id": record["run_id"],
                "previous_review_run_id": previous_record["run_id"] if previous_record else None,
            },
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("design-dir", "output-dir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    for flag in ("results", "metadata-path", "previous-review", "redesign-history"):
        parser.add_argument("--" + flag, type=Path)
    parser.add_argument("--max-rows", type=int, default=100_000)
    try:
        run_assay_review(**vars(parser.parse_args(argv)))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"assay review failed: {type(exc).__name__}: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
