"""QQ-plot expected/observed -log10(p) computation.

The original ``scripts/02_qq_plot.py`` silently dropped any ``pval`` outside
``(0, 1]`` before computing quantiles. This module does the opposite: schema
v1 validation already guarantees every ``pval`` in the input is a finite
value in ``(0, 1]`` (see :mod:`adzuki_gwas_analysis.schema`), so
:func:`compute_qq_points` treats a value outside that range as a bug to
surface loudly, not a row to quietly exclude.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class QqPoints:
    """Observed and expected -log10(p) quantiles for a QQ plot."""

    observed: np.ndarray
    expected: np.ndarray
    variant_count: int


def compute_qq_points(pvals: np.ndarray) -> QqPoints:
    """Compute sorted observed vs. expected -log10(p) for ``pvals``.

    Raises ``ValueError`` if any value is not finite or not in ``(0, 1]`` --
    on validated input this should never trigger, but a loud failure here is
    preferable to the silent row-dropping this replaces.
    """
    pvals = np.asarray(pvals, dtype="float64")
    if not np.isfinite(pvals).all() or not ((pvals > 0.0) & (pvals <= 1.0)).all():
        raise ValueError(
            "compute_qq_points requires every value to be finite and in (0, 1]; "
            "this indicates validated-data guarantees were violated upstream"
        )

    n = len(pvals)
    observed = -np.log10(np.sort(pvals))
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    return QqPoints(observed=observed, expected=expected, variant_count=n)
