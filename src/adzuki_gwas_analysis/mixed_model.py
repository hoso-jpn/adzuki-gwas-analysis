"""Small quantitative-trait LMM engine with a covariance ratio fitted under the null.

Association tests reuse that ratio (an approximation); this is not GEMMA or an
independently benchmarked replacement for its per-marker likelihood fitting.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy import optimize, stats

FloatArray = NDArray[np.float64]
ENGINE = "numpy-scipy-null-reml-lmm-v1"


@dataclass(frozen=True)
class NullModel:
    rotation: FloatArray
    weights: FloatArray
    covariates: FloatArray
    phenotype: FloatArray
    delta: float
    residual_df: int
    boundary: bool
    covariate_q: FloatArray
    residual_y: FloatArray


def fit_null(phenotype: FloatArray, covariates: FloatArray, kinship: FloatArray) -> NullModel:
    n = len(phenotype)
    if (
        phenotype.ndim != 1
        or covariates.ndim != 2
        or kinship.shape != (n, n)
        or covariates.shape[0] != n
        or not all(np.isfinite(x).all() for x in (phenotype, covariates, kinship))
    ):
        raise ValueError("invalid/nonfinite LMM inputs")
    if not np.allclose(kinship, kinship.T, rtol=1e-10, atol=1e-12):
        raise ValueError("kinship is not symmetric")
    if np.linalg.matrix_rank(covariates) != covariates.shape[1]:
        raise ValueError("covariate design is rank deficient")
    residual_df = n - covariates.shape[1] - 1
    if residual_df < 2 or np.ptp(phenotype) == 0:
        raise ValueError("insufficient residual degrees of freedom or constant phenotype")
    eigenvalues, vectors = np.linalg.eigh(kinship)
    if float(eigenvalues.min()) < -1e-8:
        raise ValueError("kinship is not positive semidefinite")
    eigenvalues = np.maximum(eigenvalues, 0)
    rotated_y, rotated_c = vectors.T @ phenotype, vectors.T @ covariates

    def objective(log_delta: float) -> float:
        variance = eigenvalues + math.exp(log_delta)
        weight = 1 / np.sqrt(variance)
        design, y = rotated_c * weight[:, None], rotated_y * weight
        coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
        residual = y - design @ coefficients
        rss = float(residual @ residual)
        sign, logdet = np.linalg.slogdet(design.T @ design)
        if rss <= 0 or sign <= 0:
            return float("inf")
        df = n - covariates.shape[1]
        return float(df * math.log(rss / df) + np.log(variance).sum() + logdet)

    result = optimize.minimize_scalar(
        objective, method="bounded", bounds=(-12, 12), options={"xatol": 1e-8}
    )
    if not result.success or not math.isfinite(float(result.fun)):
        raise ValueError("null REML optimization failed")
    delta = math.exp(float(result.x))
    weights = 1 / np.sqrt(eigenvalues + delta)
    design, whitened_y = rotated_c * weights[:, None], rotated_y * weights
    q, _ = np.linalg.qr(design, mode="reduced")
    return NullModel(
        rotation=vectors.T,
        weights=weights,
        covariates=design,
        phenotype=whitened_y,
        delta=delta,
        residual_df=residual_df,
        boundary=abs(float(result.x)) > 11.99,
        covariate_q=q,
        residual_y=whitened_y - q @ (q.T @ whitened_y),
    )


def test_marker(model: NullModel, dosage: FloatArray) -> tuple[float, float, float]:
    """Return ALT-dosage beta, SE, -log10(two-sided t p) conditional on null covariance."""
    x = (model.rotation @ dosage) * model.weights
    q = model.covariate_q
    x_residual = x - q @ (q.T @ x)
    y_residual = model.residual_y
    information = float(x_residual @ x_residual)
    if information <= np.finfo(float).eps * max(float(x @ x), np.finfo(float).tiny):
        raise ValueError("marker is collinear with covariates")
    beta = float(x_residual @ y_residual) / information
    residual = y_residual - beta * x_residual
    variance = float(residual @ residual) / model.residual_df
    if variance <= 0:
        raise ValueError("nonpositive marker residual variance")
    se = math.sqrt(variance / information)
    log_p = math.log(2) + float(stats.t.logsf(abs(beta / se), model.residual_df))
    if not math.isfinite(log_p):
        raise ValueError("association p-value tail cannot be represented reliably")
    return beta, se, max(0.0, -log_p / math.log(10))
