"""Explicit global/leave-one-chromosome-out kinship contracts for diploid dosage."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from adzuki_gwas_analysis.mixed_model import FloatArray


def kinship_matrix(dosage: FloatArray, frequency: FloatArray) -> FloatArray:
    """VanRaden method 1; dosage is already QC-filtered and mean-imputed."""
    if (
        dosage.ndim != 2
        or frequency.shape != (dosage.shape[1],)
        or dosage.shape[1] < 2
        or not np.isfinite(dosage).all()
        or not np.isfinite(frequency).all()
        or np.any((frequency <= 0) | (frequency >= 1))
    ):
        raise ValueError("kinship requires at least two finite polymorphic markers")
    centered = dosage - 2 * frequency
    return centered @ centered.T / float(2 * np.sum(frequency * (1 - frequency)))


def kinship_for_chromosome(
    dosage: FloatArray,
    frequency: FloatArray,
    chromosomes: Sequence[str],
    chromosome: str,
) -> tuple[FloatArray, int]:
    """Exclude every marker on the test chromosome, never silently fall back to global K."""
    if len(chromosomes) != dosage.shape[1] or chromosome not in chromosomes:
        raise ValueError("LOCO chromosome identity does not match genotype columns")
    keep = np.array([chrom != chromosome for chrom in chromosomes], dtype=bool)
    if int(keep.sum()) < 2:
        raise ValueError(f"LOCO for {chromosome} requires two markers on other chromosomes")
    return kinship_matrix(dosage[:, keep], frequency[keep]), int(keep.sum())
