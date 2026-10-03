"""Predeclared simulation and independent-result comparisons; never commercial certification."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

from adzuki_gwas_analysis.kinship import kinship_for_chromosome, kinship_matrix
from adzuki_gwas_analysis.mixed_model import ENGINE, fit_null, test_marker
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    output_transaction,
    write_json,
)
from adzuki_gwas_analysis.tables import read_tsv, write_tsv

SCENARIOS = ("null_unrelated", "null_related", "local_effect", "missing", "rare", "boundary")
IDENTITY = ("chr", "pos", "ref", "alt")


def wilson_interval(successes: int, total: int) -> tuple[float, float]:
    if total < 1 or not 0 <= successes <= total:
        raise ValueError("a binomial interval requires a positive independent sample count")
    z = float(stats.norm.ppf(0.975))
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    width = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, center - width), min(1.0, center + width)


def _simulation_plan(path: Path) -> dict[str, Any]:
    plan: dict[str, Any] = json.loads(path.read_text())
    if not isinstance(plan, dict):
        raise ValueError("simulation plan must be a JSON object")
    expected = {
        "schema_version",
        "seed",
        "replicates",
        "n_samples",
        "n_markers",
        "alpha",
        "max_null_upper95",
        "min_power_lower95",
        "kinship_mode",
        "scenarios",
    }
    if (
        set(plan) != expected
        or type(plan["schema_version"]) is not int
        or plan["schema_version"] != 1
    ):
        raise ValueError("unknown/incomplete simulation plan")
    for key, lower, upper in (
        ("seed", 0, 2**32 - 1),
        ("replicates", 3, 10000),
        ("n_samples", 12, 1000),
        ("n_markers", 12, 100000),
    ):
        if type(plan[key]) is not int or not lower <= plan[key] <= upper:
            raise ValueError("invalid simulation size/seed")
    if plan["n_samples"] * plan["n_markers"] > 5000000:
        raise ValueError("simulation exceeds genotype cell budget")
    for key in ("alpha", "max_null_upper95", "min_power_lower95"):
        if (
            type(plan[key]) not in (float, int)
            or not math.isfinite(plan[key])
            or not 0 < plan[key] < 1
        ):
            raise ValueError("invalid simulation rate criterion")
    if plan["kinship_mode"] not in ("global", "loco") or plan["scenarios"] != list(SCENARIOS):
        raise ValueError("declare global/loco and the complete ordered scenario set")
    return plan


def _replicate(plan: dict[str, Any], scenario: str, seed: np.random.SeedSequence) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    n, m = plan["n_samples"], plan["n_markers"]
    group = np.arange(n) % 2
    probabilities = np.full((n, m), 0.35)
    if scenario in ("null_related", "local_effect"):
        probabilities += (group[:, None] - 0.5) * 0.4
    if scenario == "rare":
        probabilities[:, 0] = 0.08
    dosage = rng.binomial(2, probabilities).astype(float)
    chromosomes = [f"chr{j % 3}" for j in range(m)]
    if scenario == "missing":
        dosage[rng.uniform(size=dosage.shape) < 0.05] = np.nan
    frequencies = np.nanmean(dosage, axis=0) / 2
    keep = np.isfinite(frequencies) & (frequencies > 0) & (frequencies < 1)
    if not keep[0]:
        return {"scenario": scenario, "status": "not_tested", "reason": "target_monomorphic"}
    dosage, frequencies = dosage[:, keep], frequencies[keep]
    chromosomes = [c for c, selected in zip(chromosomes, keep, strict=True) if selected]
    dosage = np.where(np.isnan(dosage), 2 * frequencies, dosage)
    global_k = kinship_matrix(dosage, frequencies)
    background_k, _ = kinship_for_chromosome(dosage, frequencies, chromosomes, chromosomes[0])
    generation_k = np.eye(n) if scenario == "null_unrelated" else background_k
    values, vectors = np.linalg.eigh(generation_k)
    genetic = vectors @ (np.sqrt(np.maximum(values, 0)) * rng.normal(size=n))
    if scenario == "boundary":
        genetic *= 0
    true_beta = 1.5 if scenario == "local_effect" else 0.0
    y = true_beta * dosage[:, 0] + genetic + rng.normal(size=n)
    fitted_k = global_k if plan["kinship_mode"] == "global" else background_k
    model = fit_null(y, np.ones((n, 1)), fitted_k)
    beta, se, logp = test_marker(model, dosage[:, 0])
    critical = float(stats.t.ppf(0.975, model.residual_df))
    return {
        "scenario": scenario,
        "status": "tested",
        "reason": "",
        "true_beta": true_beta,
        "beta": beta,
        "standard_error": se,
        "neg_log10_pvalue": logp,
        "reject": logp >= -math.log10(plan["alpha"]),
        "effect_covered": beta - critical * se <= true_beta <= beta + critical * se,
        "null_reml_boundary": model.boundary,
    }


def simulate(plan_path: Path, output_dir: Path) -> None:
    plan = _simulation_plan(plan_path)
    inputs = {"simulation_plan": plan_path}
    initial, environment = input_checksums(inputs), generation_environment()
    seeds = np.random.SeedSequence(plan["seed"]).spawn(len(SCENARIOS) * plan["replicates"])
    rows, summaries = [], []
    for index, scenario in enumerate(SCENARIOS):
        scenario_rows = []
        for repetition in range(plan["replicates"]):
            result = _replicate(plan, scenario, seeds[index * plan["replicates"] + repetition])
            result["replicate"] = repetition
            rows.append(result)
            scenario_rows.append(result)
        tested = [row for row in scenario_rows if row["status"] == "tested"]
        rejected = sum(row["reject"] for row in tested)
        lower, upper = wilson_interval(rejected, len(tested)) if tested else (0.0, 1.0)
        criterion = (
            lower >= plan["min_power_lower95"]
            if scenario == "local_effect"
            else upper <= plan["max_null_upper95"]
        )
        summaries.append(
            {
                "scenario": scenario,
                "independent_replicates": len(tested),
                "not_tested": len(scenario_rows) - len(tested),
                "rejections": rejected,
                "rate": rejected / len(tested) if tested else None,
                "rate_lower95": lower,
                "rate_upper95": upper,
                "criterion_met": bool(criterion and len(tested) == len(scenario_rows)),
                "mean_effect_bias": float(np.mean([r["beta"] - r["true_beta"] for r in tested]))
                if tested
                else None,
                "effect_interval_coverage": sum(r["effect_covered"] for r in tested) / len(tested)
                if tested
                else None,
                "boundary_count": sum(r["null_reml_boundary"] for r in tested),
            }
        )
    with output_transaction(output_dir) as stage:
        write_json(stage / "simulation_plan.json", plan)
        fields = (
            "scenario",
            "replicate",
            "status",
            "reason",
            "true_beta",
            "beta",
            "standard_error",
            "neg_log10_pvalue",
            "reject",
            "effect_covered",
            "null_reml_boundary",
        )
        write_tsv(stage / "replicates.tsv", fields, rows)
        write_json(
            stage / "calibration.json",
            {
                "schema_version": 1,
                "engine": ENGINE,
                "data_scope": "synthetic",
                "scenarios": summaries,
                "simulation_criteria_met": all(s["criterion_met"] for s in summaries),
                "commercial_status": "not_validated",
                "remaining": [
                    "independent_engine_comparison",
                    "real_cohort_reference_validation",
                    "intended_scale_performance",
                    "analyst_review",
                ],
                "scope": "one prespecified marker per independent simulated replicate; "
                "not genome-wide FWER or customer performance",
                "generation": "Gaussian random effects from declared background kinship; "
                "no phenotype permutation",
            },
        )
        finish_provenance(
            stage,
            operation="gwas_calibration",
            inputs=inputs,
            initial=initial,
            parameters=plan,
            environment=environment,
            families=[],
        )


def compare_results(left: Path, right: Path, plan_path: Path, output_dir: Path) -> None:
    """Compare normalized outputs only under an explicitly shared statistical contract."""
    plan: dict[str, Any] = json.loads(plan_path.read_text())
    if not isinstance(plan, dict):
        raise ValueError("comparison plan must be a JSON object")
    expected = {
        "schema_version",
        "left_sha256",
        "right_sha256",
        "left_contract",
        "right_contract",
        "beta_atol",
        "se_atol",
        "logp_atol",
    }
    if (
        set(plan) != expected
        or type(plan["schema_version"]) is not int
        or plan["schema_version"] != 1
    ):
        raise ValueError("unknown/incomplete comparison plan")
    contract_fields = {
        "dataset_id",
        "assembly_id",
        "input_genotypes_sha256",
        "input_phenotypes_sha256",
        "covariates_sha256",
        "model",
        "covariance_strategy",
        "test",
        "effect_encoding",
        "engine",
        "version",
        "license",
        "data_scope",
    }
    for name in ("left_contract", "right_contract"):
        contract = plan[name]
        if (
            not isinstance(contract, dict)
            or set(contract) != contract_fields
            or any(not isinstance(v, str) or not v for v in contract.values())
        ):
            raise ValueError("comparison requires explicit engine, input and model contracts")
        if contract["data_scope"] not in ("synthetic", "public", "customer"):
            raise ValueError("invalid comparison data scope")
    for field in ("beta_atol", "se_atol", "logp_atol"):
        if (
            type(plan[field]) not in (int, float)
            or not math.isfinite(plan[field])
            or plan[field] < 0
        ):
            raise ValueError("comparison tolerances must be finite nonnegative numbers")
    inputs = {"left": left, "right": right, "comparison_plan": plan_path}
    initial, environment = input_checksums(inputs), generation_environment()
    if initial["left"] != plan["left_sha256"] or initial["right"] != plan["right_sha256"]:
        raise ValueError("comparison source changed from predeclared plan")
    datasets = []
    for path in (left, right):
        records = {}
        for row in read_tsv(
            path,
            (
                *IDENTITY,
                "effect_allele",
                "other_allele",
                "beta",
                "se",
                "neg_log10_pvalue",
            ),
        ):
            identity = (row["chr"], int(row["pos"]), row["ref"], row["alt"])
            if identity in records:
                raise ValueError("duplicate comparison variant identity")
            for field in ("beta", "se", "neg_log10_pvalue"):
                if not math.isfinite(float(row[field])):
                    raise ValueError("nonfinite comparison value")
            if float(row["se"]) <= 0 or float(row["neg_log10_pvalue"]) < 0:
                raise ValueError("invalid comparison SE/log-p")
            records[identity] = row
        datasets.append(records)
    a, b = datasets
    if not a or set(a) != set(b):
        raise ValueError("comparison variant sets must be nonempty and identical")
    shared_fields = contract_fields - {"engine", "version", "license"}
    differences = [
        field
        for field in sorted(shared_fields)
        if plan["left_contract"][field] != plan["right_contract"][field]
    ]
    rows = []
    for identity, first in a.items():
        second = b[identity]
        if any(first[field] != second[field] for field in ("effect_allele", "other_allele")):
            raise ValueError("comparison effect orientation mismatch; normalize explicitly first")
        delta = {
            field: abs(float(first[field]) - float(second[field]))
            for field in ("beta", "se", "neg_log10_pvalue")
        }
        rows.append(
            {
                **{field: first[field] for field in IDENTITY},
                **delta,
                "within_tolerance": delta["beta"] <= plan["beta_atol"]
                and delta["se"] <= plan["se_atol"]
                and delta["neg_log10_pvalue"] <= plan["logp_atol"],
            }
        )
    independent = plan["left_contract"]["engine"].strip().casefold() != (
        plan["right_contract"]["engine"].strip().casefold()
    )
    status = (
        "not_comparable"
        if differences
        else "passed"
        if independent and all(row["within_tolerance"] for row in rows)
        else "failed"
    )
    with output_transaction(output_dir) as stage:
        write_tsv(stage / "differences.tsv", rows[0].keys(), rows)
        write_json(
            stage / "comparison.json",
            {
                "schema_version": 1,
                "status": status,
                "independent_engine": independent,
                "contract_differences": differences,
                "n_variants": len(rows),
                "commercial_status": "not_validated",
                "note": "Artifact comparison is not proof of external execution or commercial "
                "release; preserve original run provenance and obtain analyst review.",
            },
        )
        finish_provenance(
            stage,
            operation="gwas_engine_comparison",
            inputs=inputs,
            initial=initial,
            parameters=plan,
            environment=environment,
            families=[],
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("simulate", "compare"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--plan-path", type=Path, required=True)
        sub.add_argument("--output-dir", type=Path, required=True)
        if command == "compare":
            sub.add_argument("--left", type=Path, required=True)
            sub.add_argument("--right", type=Path, required=True)
    args = vars(parser.parse_args(argv))
    command = args.pop("command")
    try:
        (simulate if command == "simulate" else compare_results)(**args)
    except (ValueError, KeyError, OSError, np.linalg.LinAlgError) as exc:
        parser.exit(1, f"GWAS validation failed: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
