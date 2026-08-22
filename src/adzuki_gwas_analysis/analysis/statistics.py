"""Pure numeric multiple-testing correction and genomic-inflation diagnostics.

Every function here operates on an in-memory ``numpy`` array of p-values and
performs no file I/O. Callers are expected to pass p-values that already
satisfy schema v1's contract (finite, ``0 < p <= 1``), but -- following the
precedent set by :func:`adzuki_gwas_analysis.analysis.qq.compute_qq_points`
-- this module never silently drops or coerces an out-of-contract value; it
raises ``ValueError`` loudly instead, since a violation here indicates a bug
upstream, not a row to quietly exclude.

Bonferroni is a family-wise error rate (FWER) threshold for the multiple-
testing family the caller passes in (see
:mod:`adzuki_gwas_analysis.analysis.diagnostics` for how that family is
defined for this repository) -- it is not a universal genome-wide-
significance threshold, and it is not adjusted for linkage disequilibrium
among markers (it may therefore be conservative for LD-correlated SNPs).

Benjamini-Hochberg (BH) controls the false discovery rate (FDR) via
``scipy.stats.false_discovery_control``. The original Benjamini & Hochberg
(1995) guarantee holds under independence; Benjamini & Yekutieli (2001)
extended it to positive regression dependency on a subset (PRDS) -- a
narrower condition than "any dependence structure." This module does not
claim FDR control under arbitrary dependence, and does not implement the
Benjamini-Yekutieli correction.

The genomic inflation factor (lambda_GC) is reported as a diagnostic value
only: nothing here divides a test statistic or p-value by it, and lambda_GC
alone cannot distinguish population stratification/kinship/batch effects
from polygenicity as the cause of any inflation (Bulik-Sullivan et al. 2015).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


def _validate_pvalues(pvalues: np.ndarray) -> np.ndarray:
    """Return a private float64 copy of ``pvalues``, validated and never aliased.

    A copy (not ``np.asarray``) is made deliberately: if ``pvalues`` is
    already a float64 array, ``np.asarray`` would return the same object,
    and any in-place operation downstream (in this module or inside SciPy)
    could then mutate the caller's original array. Copying here guarantees
    every function in this module leaves its input untouched.
    """
    array = np.array(pvalues, dtype="float64", copy=True)
    if array.size == 0:
        raise ValueError("pvalues must not be empty")
    if not np.isfinite(array).all() or not ((array > 0.0) & (array <= 1.0)).all():
        raise ValueError(
            "pvalues must all be finite and satisfy 0 < p <= 1; this indicates "
            "validated-data guarantees were violated upstream"
        )
    return array


def _validate_unit_interval(value: float, *, name: str) -> None:
    if not (0.0 < value <= 1.0):
        raise ValueError(f"{name} must satisfy 0 < {name} <= 1, got {value}")


@dataclass(frozen=True, slots=True)
class BonferroniResult:
    """Bonferroni FWER correction over one multiple-testing family."""

    m: int
    alpha: float
    threshold: float
    adjusted: np.ndarray
    significant: np.ndarray
    discoveries: int


def compute_bonferroni(pvalues: np.ndarray, *, alpha: float = 0.05) -> BonferroniResult:
    """Compute the Bonferroni threshold, adjusted p-values, and discovery mask.

    ``threshold = alpha / m`` and ``adjusted = min(p * m, 1.0)`` are
    equivalent significance tests (``p <= alpha / m`` iff
    ``min(p * m, 1.0) <= alpha``); both are returned so callers never need to
    re-derive one from the other.
    """
    validated = _validate_pvalues(pvalues)
    _validate_unit_interval(alpha, name="alpha")

    m = validated.size
    threshold = alpha / m
    adjusted = np.minimum(validated * m, 1.0)
    significant = validated <= threshold
    return BonferroniResult(
        m=m,
        alpha=alpha,
        threshold=threshold,
        adjusted=adjusted,
        significant=significant,
        discoveries=int(significant.sum()),
    )


@dataclass(frozen=True, slots=True)
class BhResult:
    """Benjamini-Hochberg FDR correction over one multiple-testing family."""

    m: int
    fdr_level: float
    adjusted: np.ndarray
    significant: np.ndarray
    discoveries: int
    raw_p_cutoff: float | None


def compute_bh(pvalues: np.ndarray, *, fdr_level: float = 0.05) -> BhResult:
    """Compute BH-adjusted p-values, the discovery mask, and the raw-p cutoff.

    Uses :func:`scipy.stats.false_discovery_control` (``method="bh"``), which
    returns adjusted p-values aligned to the original input order (it sorts
    internally for the step-up procedure, then restores input order before
    returning -- callers of this function never see or need to undo any
    reordering). ``raw_p_cutoff`` is the largest raw p-value among rejected
    variants, or ``None`` if there are zero discoveries -- never a crash or a
    silently-omitted field.

    The adjusted p-value returned here is a BH-adjusted p-value, not a
    Storey-style q-value (a different estimator); callers must not relabel
    it as a "q-value".
    """
    validated = _validate_pvalues(pvalues)
    _validate_unit_interval(fdr_level, name="fdr_level")

    adjusted = stats.false_discovery_control(validated, method="bh")
    significant = adjusted <= fdr_level
    discoveries = int(significant.sum())
    raw_p_cutoff = float(validated[significant].max()) if discoveries > 0 else None
    return BhResult(
        m=validated.size,
        fdr_level=fdr_level,
        adjusted=adjusted,
        significant=significant,
        discoveries=discoveries,
        raw_p_cutoff=raw_p_cutoff,
    )


@dataclass(frozen=True, slots=True)
class LambdaGcResult:
    """Genomic inflation factor (lambda_GC), a diagnostic value only.

    Nothing in this repository divides a test statistic or p-value by
    ``lambda_gc``, and this value alone must not be read as identifying (or
    ruling out) population stratification, kinship, or batch effects --
    polygenicity also inflates GWAS test statistics and is not
    distinguishable from confounding by lambda_GC alone.
    """

    df: int
    expected_median: float
    lambda_gc: float


def compute_lambda_gc(pvalues: np.ndarray, *, df: int = 1) -> LambdaGcResult:
    """Compute lambda_GC = median(chi2.isf(p, df)) / chi2.ppf(0.5, df).

    Uses the survival-function side (``isf``) rather than
    ``chi2.ppf(1 - p, df)`` to avoid catastrophic cancellation at small
    p-values, and never substitutes the common ``-2 * log(p)`` approximation
    (which is only equal to the 1-df chi-square inverse survival function
    exactly at ``df=1`` in a different derivation, not a general substitute).
    ``expected_median`` is computed via ``chi2.ppf``, not rounded to a fixed
    literal, so it reflects the ``df`` actually used. ``p == 1`` maps to
    ``chi2 == 0`` (``isf(1.0, df) == 0.0``).
    """
    validated = _validate_pvalues(pvalues)
    if not isinstance(df, int) or isinstance(df, bool) or df <= 0:
        raise ValueError(f"df must be a positive integer, got {df!r}")

    chi2_values = stats.chi2.isf(validated, df=df)
    if not np.isfinite(chi2_values).all():
        raise ValueError(
            "chi2.isf produced a non-finite value from otherwise-valid pvalues; "
            "this indicates a numerical contract violation, not a value to silently drop"
        )
    expected_median = float(stats.chi2.ppf(0.5, df=df))
    lambda_gc = float(np.median(chi2_values) / expected_median)
    return LambdaGcResult(df=df, expected_median=expected_median, lambda_gc=lambda_gc)
