"""Run a bounded, dense quantitative-trait LMM on explicitly coded individual data."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import numpy as np

from adzuki_gwas_analysis.customer_summary import (
    REQUIRED_COLUMNS,
    adjust_and_cluster,
    normalize_row,
)
from adzuki_gwas_analysis.individual_inputs import load_config, load_individual_inputs
from adzuki_gwas_analysis.mixed_model import ENGINE, fit_null, test_marker
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    output_transaction,
    write_json,
)
from adzuki_gwas_analysis.reference import load_reference_bundle
from adzuki_gwas_analysis.summary_reporting import write_summary_artifacts
from adzuki_gwas_analysis.tables import write_tsv


def run_individual_gwas(
    *,
    genotypes: Path,
    phenotypes: Path,
    samples: Path,
    variants: Path,
    config_path: Path,
    bundle_path: Path,
    output_dir: Path,
) -> None:
    bundle = load_reference_bundle(bundle_path)
    config = load_config(config_path, bundle)
    inputs = {
        "genotypes": genotypes,
        "phenotypes": phenotypes,
        "samples": samples,
        "variants": variants,
        "config": config_path,
        "reference_bundle": bundle_path,
        "reference_fasta": bundle.fasta,
    }
    initial, environment = input_checksums(inputs), generation_environment()
    data = load_individual_inputs(
        genotypes=genotypes,
        phenotypes=phenotypes,
        samples=samples,
        variants=variants,
        bundle=bundle,
        config=config,
    )
    called = np.sum(~np.isnan(data.dosage), axis=0)
    frequency = np.divide(
        np.nansum(data.dosage, axis=0), 2 * called, out=np.zeros(len(called)), where=called > 0
    )
    missing = 1 - called / len(data.sample_ids)
    retained: list[int] = []
    marker_qc: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for index, marker in enumerate(data.markers):
        reasons = []
        if called[index] == 0:
            reasons.append("no_observed_calls")
        if frequency[index] in (0, 1):
            reasons.append("monomorphic")
        if min(frequency[index], 1 - frequency[index]) < config["min_maf"]:
            reasons.append("below_min_maf")
        if missing[index] > config["max_marker_missing"]:
            reasons.append("marker_missingness")
        row = {
            "marker_id": marker["marker_id"],
            "source_row": index + 2,
            "observed_calls": int(called[index]),
            "missing_fraction": float(missing[index]),
            "alt_frequency": float(frequency[index]),
            "status": "excluded" if reasons else "retained",
            "reason": ";".join(reasons),
        }
        marker_qc.append(row)
        if reasons:
            excluded.append({"source_row": index + 2, "reason": row["reason"]})
        else:
            retained.append(index)
    if len(retained) < 2:
        raise ValueError("fewer than two polymorphic markers remain after QC")
    af = frequency[retained]
    genotype = data.dosage[:, retained].copy()
    genotype = np.where(np.isnan(genotype), 2 * af, genotype)
    centered = genotype - 2 * af
    kinship = centered @ centered.T / float(2 * np.sum(af * (1 - af)))
    covariates = data.covariates
    pcs = np.empty((len(data.sample_ids), 0))
    eigenvalues: list[float] = []
    if config["n_pcs"]:
        values, vectors = np.linalg.eigh(kinship)
        if config["n_pcs"] > np.count_nonzero(values > max(float(values.max()), 1) * 1e-10):
            raise ValueError("requested PCs exceed positive kinship rank")
        pcs = vectors[:, ::-1][:, : config["n_pcs"]].copy()
        # Resolve the arbitrary sign so exported scores are stable on a fixed eigensystem.
        for col in range(pcs.shape[1]):
            if pcs[int(np.argmax(np.abs(pcs[:, col]))), col] < 0:
                pcs[:, col] *= -1
        eigenvalues = [float(v) for v in values[::-1][: config["n_pcs"]]]
        covariates = np.column_stack((covariates, pcs))
    model = fit_null(data.phenotype, covariates, kinship)
    meta = {
        key: config[key]
        for key in (
            "dataset_id",
            "cohort_id",
            "analysis_id",
            "reference",
            "assembly_id",
            "species",
            "trait",
            "trait_unit",
            "trait_coding",
        )
    }
    meta.update(
        test="null_REML_LMM_t_approximation",
        effect_type="beta",
        effect_scale=config["trait_unit"],
        standard_error_scale=config["trait_unit"],
        effect_strand="forward",
        target_effect_allele="ALT",
        sample_size_definition="analyzed_individuals_after_mean_imputation",
        columns={
            key: key for key in (*REQUIRED_COLUMNS, "effect_allele_frequency", "neg_log10_pvalue")
        },
    )
    rows: list[dict[str, Any]] = []
    associations: list[dict[str, Any]] = []
    for offset, index in enumerate(retained):
        marker = data.markers[index]
        try:
            beta, se, logp = test_marker(model, genotype[:, offset])
        except ValueError as exc:
            if str(exc) != "marker is collinear with covariates":
                raise
            marker_qc[index].update(status="excluded", reason="collinear_with_covariates")
            excluded.append({"source_row": index + 2, "reason": "collinear_with_covariates"})
            continue
        association = {
            **marker,
            "effect_allele": marker["alt"],
            "other_allele": marker["ref"],
            "effect": beta,
            "standard_error": se,
            "pvalue": math.pow(10, -logp),
            "neg_log10_pvalue": logp,
            "sample_size": len(data.sample_ids),
            "effect_allele_frequency": float(af[offset]),
        }
        associations.append(association)
        # Underflowed p is explicit: normalize_row requires the original log p and never clamps it.
        row = normalize_row(
            {k: str(v) for k, v in association.items()},
            meta,
            bundle,
            row_number=index + 2,
            source_sha256=initial["genotypes"],
        )
        row["pvalue_source"] = "individual_engine_log_t_tail"
        rows.append(row)
        marker_qc[index]["status"] = "tested"
    contigs = [contig.name for contig in bundle.contigs]
    candidates = adjust_and_cluster(
        rows,
        alpha=config["alpha"],
        fdr_level=config["fdr_level"],
        clustering_distance=config["clustering_distance"],
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
            alpha=config["alpha"],
            fdr_level=config["fdr_level"],
            analysis_origin="individual_gwas",
        )
        write_tsv(stage / "association_results.tsv", associations[0].keys(), associations)
        write_tsv(stage / "marker_qc.tsv", marker_qc[0].keys(), marker_qc)
        write_tsv(stage / "sample_qc.tsv", data.sample_qc[0].keys(), data.sample_qc)
        np.save(stage / "kinship.npy", kinship, allow_pickle=False)
        write_tsv(
            stage / "sample_order.tsv", ("sample_id",), [{"sample_id": i} for i in data.sample_ids]
        )
        write_tsv(
            stage / "pc_scores.tsv",
            ("sample_id", *(f"PC{i + 1}" for i in range(pcs.shape[1]))),
            [
                {
                    "sample_id": identifier,
                    **{f"PC{col + 1}": float(pcs[row, col]) for col in range(pcs.shape[1])},
                }
                for row, identifier in enumerate(data.sample_ids)
            ],
        )
        write_json(
            stage / "model.json",
            {
                "schema_version": 1,
                "engine": ENGINE,
                "model": "y = intercept + covariates + PC + beta*ALT_dosage + u + e",
                "kinship": "global VanRaden method 1: ZZ' / (2 sum(p*(1-p)))",
                "imputation": "per-marker observed-cohort mean ALT dosage",
                "covariance_ratio": model.delta,
                "ratio_bounds_log": [-12, 12],
                "null_reml_boundary": model.boundary,
                "residual_df": model.residual_df,
                "pvalue": "two-sided t; covariance ratio fitted under null and held fixed",
                "covariate_transform": data.covariate_transform,
                "n_pcs": config["n_pcs"],
                "pc_eigenvalues": eigenvalues,
                "seed": "not_applicable_no_stochastic_steps",
                "n_samples": len(data.sample_ids),
                "n_input_markers": len(data.markers),
                "n_kinship_markers": len(retained),
                "n_tests": len(rows),
                "warnings": [
                    "approximate covariance; no per-marker REML or LOCO",
                    "global K can cause proximal contamination",
                    "no real-cohort or Seedcore-01 validation performed",
                ],
                "configuration": config,
            },
        )
        write_json(stage / "statistical_diagnostics.json", diagnostics)
        finish_provenance(
            stage,
            operation="individual_gwas",
            inputs=inputs,
            initial=initial,
            environment=environment,
            families=[diagnostics["family"]],
            parameters={
                "engine": ENGINE,
                "configuration": config,
                "n_input": len(data.markers),
                "n_tests": len(rows),
                "n_excluded": len(excluded),
                "reference": bundle.metadata(),
            },
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in (
        "genotypes",
        "phenotypes",
        "samples",
        "variants",
        "config-path",
        "bundle-path",
        "output-dir",
    ):
        parser.add_argument("--" + flag, type=Path, required=True)
    try:
        run_individual_gwas(**vars(parser.parse_args(argv)))
    except (ValueError, KeyError, OSError, np.linalg.LinAlgError) as exc:
        parser.exit(1, f"individual GWAS failed: {type(exc).__name__}: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
