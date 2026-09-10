"""One quantitative observation per sample, explicit ALT dosage and trial covariates."""

from __future__ import annotations

import itertools
import math
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from adzuki_gwas_analysis.mixed_model import FloatArray
from adzuki_gwas_analysis.reference import ReferenceBundle, required_text
from adzuki_gwas_analysis.tables import read_tsv

TEXT_FIELDS = (
    "dataset_id",
    "cohort_id",
    "analysis_id",
    "reference",
    "assembly_id",
    "species",
    "trait",
    "trait_unit",
    "trait_coding",
    "trait_type",
    "design_note",
    "genotype_encoding",
    "permission_reference",
    "transfer_method",
    "retention_policy",
    "deletion_policy",
    "output_ownership",
)
DEFAULTS: dict[str, Any] = {
    "min_maf": 0.05,
    "max_marker_missing": 0.1,
    "max_sample_missing": 0.1,
    "max_samples": 1000,
    "max_markers": 100000,
    "max_genotype_cells": 5000000,
    "n_pcs": 0,
    "covariates": [],
    "alpha": 0.05,
    "fdr_level": 0.05,
}


def load_config(path: Path, bundle: ReferenceBundle) -> dict[str, Any]:
    config = tomllib.loads(path.read_text())
    if type(config.get("schema_version")) is not int or config["schema_version"] != 1:
        raise ValueError("unsupported individual GWAS config schema")
    if config.get("data_use_confirmed") is not True or config.get("data_scope") not in (
        "synthetic",
        "customer",
        "public",
    ):
        raise ValueError("data scope and confirmed permission are required")
    for key in TEXT_FIELDS:
        required_text(config, key)
    if config["trait_type"] != "quantitative" or config["genotype_encoding"] != "ALT_0_1_2":
        raise ValueError("v1 requires quantitative traits and explicit ALT_0_1_2 dosage")
    bundle.check_dataset(config["dataset_id"], config["reference"], config["assembly_id"])
    if config["species"] != bundle.species:
        raise ValueError("species mismatch")
    if (
        set(config)
        - set(TEXT_FIELDS)
        - set(DEFAULTS)
        - {"schema_version", "data_use_confirmed", "data_scope", "clustering_distance"}
    ):
        raise ValueError("unknown individual GWAS config field")
    config = {**DEFAULTS, **config}
    for key in ("max_samples", "max_markers", "max_genotype_cells", "n_pcs", "clustering_distance"):
        if type(config.get(key)) is not int or config[key] < 0:
            raise ValueError(f"{key} requires an explicit nonnegative integer")
    if config["max_samples"] < 4 or config["max_markers"] < 2:
        raise ValueError("sample/marker limits are too small")
    for key in ("min_maf", "max_marker_missing", "max_sample_missing", "alpha", "fdr_level"):
        value = config[key]
        if (
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(value)
        ):
            raise ValueError(f"invalid {key}")
    if not 0 <= config["min_maf"] <= 0.5 or not all(
        0 <= config[k] < 1 for k in ("max_marker_missing", "max_sample_missing")
    ):
        raise ValueError("QC thresholds are outside supported fractions")
    if not all(0 < config[k] < 1 for k in ("alpha", "fdr_level")):
        raise ValueError("invalid multiple testing thresholds")
    covariates = config["covariates"]
    if (
        not isinstance(covariates, list)
        or any(not isinstance(c, str) or not c for c in covariates)
        or len(set(covariates)) != len(covariates)
    ):
        raise ValueError("covariates must list unique column names")
    return config


def _unique_rows(
    path: Path, key: str, required: tuple[str, ...], limit: int
) -> dict[str, dict[str, str]]:
    result = {}
    for row in read_tsv(path, required):
        identifier = row[key]
        if not identifier or identifier in result:
            raise ValueError(f"empty/duplicate {key}")
        if len(result) >= limit:
            raise ValueError("input exceeds declared operating limit")
        result[identifier] = row
    return result


@dataclass
class IndividualInputs:
    sample_ids: list[str]
    markers: list[dict[str, str]]
    dosage: FloatArray
    phenotype: FloatArray
    covariates: FloatArray
    sample_qc: list[dict[str, Any]]
    covariate_transform: dict[str, dict[str, float]]


def load_individual_inputs(
    *,
    genotypes: Path,
    phenotypes: Path,
    samples: Path,
    variants: Path,
    bundle: ReferenceBundle,
    config: dict[str, Any],
) -> IndividualInputs:
    sample_rows = _unique_rows(
        samples, "sample_id", ("sample_id", "cohort_id"), config["max_samples"]
    )
    sample_ids = list(sample_rows)
    if len(sample_ids) < 4 or any(
        row["cohort_id"] != config["cohort_id"] for row in sample_rows.values()
    ):
        raise ValueError("insufficient samples or mixed cohort IDs")
    phenotype_rows = _unique_rows(
        phenotypes, "sample_id", ("sample_id", "trait", "value", "unit"), config["max_samples"]
    )
    if set(phenotype_rows) != set(sample_ids):
        raise ValueError("genotype/phenotype/sample ID sets must match exactly")
    if any(
        row["trait"] != config["trait"] or row["unit"] != config["trait_unit"]
        for row in phenotype_rows.values()
    ):
        raise ValueError("mixed trait or phenotype unit")
    y = np.array(
        [float(phenotype_rows[identifier]["value"]) for identifier in sample_ids], dtype=float
    )
    if not np.isfinite(y).all() or len(np.unique(y)) < 3:
        raise ValueError("phenotype must contain finite, varying quantitative observations")
    variant_rows = _unique_rows(
        variants,
        "marker_id",
        ("marker_id", "chr", "pos", "ref", "alt", "assembly_id"),
        config["max_markers"],
    )
    reader = read_tsv(genotypes, ("sample_id",))
    first = next(reader, None)
    if first is None:
        raise ValueError("genotype matrix is empty")
    marker_ids = [column for column in first if column != "sample_id"]
    if set(marker_ids) != set(variant_rows) or not marker_ids:
        raise ValueError("genotype marker IDs disagree with the variant metadata")
    markers = [variant_rows[identifier] for identifier in marker_ids]
    if len(sample_ids) * len(marker_ids) > config["max_genotype_cells"]:
        raise ValueError("genotype matrix exceeds declared cell limit")
    loci = set()
    for marker in markers:
        if marker["assembly_id"] != bundle.assembly_id:
            raise ValueError("variant assembly mismatch")
        bundle.check_snp(marker["chr"], int(marker["pos"]), marker["ref"], marker["alt"])
        locus = (marker["chr"], int(marker["pos"]), marker["ref"], marker["alt"])
        if locus in loci:
            raise ValueError("duplicate variant identity")
        loci.add(locus)
    dosage = np.empty((len(sample_ids), len(marker_ids)), dtype=float)
    sample_index = {identifier: index for index, identifier in enumerate(sample_ids)}
    seen = set()
    calls = {"0": 0.0, "1": 1.0, "2": 2.0, "NA": float("nan")}
    for row in itertools.chain((first,), reader):
        identifier = row["sample_id"]
        if identifier in seen or identifier not in sample_index:
            raise ValueError("duplicate or unknown genotype sample ID")
        seen.add(identifier)
        try:
            dosage[sample_index[identifier]] = np.fromiter(
                (calls[row[marker]] for marker in marker_ids), dtype=float
            )
        except KeyError as exc:
            raise ValueError("genotype cells must be literal 0/1/2/NA ALT dosages") from exc
    if seen != set(sample_ids):
        raise ValueError("genotype sample IDs are missing")
    missing = np.isnan(dosage).mean(axis=1)
    if np.any(missing > config["max_sample_missing"]):
        raise ValueError("sample missingness exceeds the declared limit")
    covariates = [np.ones(len(sample_ids), dtype=float)]
    transform = {}
    for name in config["covariates"]:
        values = np.array([float(sample_rows[i][name]) for i in sample_ids], dtype=float)
        if not np.isfinite(values).all() or float(values.std()) == 0:
            raise ValueError("covariate is nonfinite or constant; supply numeric design columns")
        mean, scale = float(values.mean()), float(values.std())
        transform[name] = {"mean": mean, "scale": scale}
        covariates.append((values - mean) / scale)
    return IndividualInputs(
        sample_ids,
        markers,
        dosage,
        y,
        np.column_stack(covariates),
        [{"sample_id": i, "missing_fraction": float(missing[j])} for j, i in enumerate(sample_ids)],
        transform,
    )
