"""Matplotlib rendering for Manhattan, QQ, and regional plots.

Headless-safe: this module never calls ``matplotlib.use()`` itself (that
would clobber a caller's own backend choice). Callers that need headless
rendering (CI, tests, batch smoke tests) set ``MPLBACKEND=Agg`` in the
environment before matplotlib is imported anywhere in the process.

Every ``plot_*`` function writes to a temporary file in the same directory
as its final destination and then ``os.replace``s it into place, so a
failure partway through rendering never leaves a truncated/corrupt PNG at
the destination path.

Figure sizes intentionally match the original scripts exactly (Manhattan:
matplotlib's default 6.4x4.8in; QQ: 6x6in; regional: 8x4in) -- three
different plot types with three different existing figsizes, not one shared
size.
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from adzuki_gwas_analysis.analysis.chromosomes import ManhattanCoordinates
from adzuki_gwas_analysis.analysis.qq import QqPoints
from adzuki_gwas_analysis.analysis.regions import TopVariant

DEFAULT_DPI = 300


@dataclass(frozen=True, slots=True)
class PlotResult:
    """Metadata about a saved plot, for logging/equivalence checks (not the image itself)."""

    output_path: Path
    width_px: int
    height_px: int
    dpi: int


def _atomic_savefig(fig: Any, output_path: Path, *, dpi: int) -> PlotResult:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=output_path.parent, suffix=".png.tmp")
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        fig.savefig(tmp_path, dpi=dpi, format="png")
        os.replace(tmp_path, output_path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise
    finally:
        width_in, height_in = fig.get_size_inches()
        plt.close(fig)
    return PlotResult(
        output_path=output_path,
        width_px=round(width_in * dpi),
        height_px=round(height_in * dpi),
        dpi=dpi,
    )


def plot_manhattan(
    coords: ManhattanCoordinates,
    pval: pd.Series,
    *,
    threshold: float,
    output_path: Path,
    title: str = "Water Permeability GWAS",
    dpi: int = DEFAULT_DPI,
) -> PlotResult:
    """Render a genome-wide Manhattan plot; matches ``scripts/01_manhattan_plot.py``."""
    fig, ax = plt.subplots()
    ax.scatter(coords.x, -np.log10(pval.to_numpy(dtype="float64")), s=2, alpha=0.6)
    ax.axhline(-np.log10(threshold), linestyle="--")
    ax.set_xticks(coords.tick_positions)
    ax.set_xticklabels(coords.tick_labels, rotation=0)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    fig.tight_layout()
    return _atomic_savefig(fig, output_path, dpi=dpi)


def plot_qq(
    qq_points: QqPoints,
    *,
    output_path: Path,
    title: str = "QQ Plot: Water Permeability GWAS",
    dpi: int = DEFAULT_DPI,
) -> PlotResult:
    """Render a QQ plot; matches ``scripts/02_qq_plot.py``."""
    fig = plt.figure(figsize=(6, 6))
    ax = fig.gca()
    ax.scatter(qq_points.expected, qq_points.observed, s=3, alpha=0.5)

    max_val = max(qq_points.expected.max(), qq_points.observed.max())
    ax.plot([0, max_val], [0, max_val], linestyle="--")

    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title(title)
    fig.tight_layout()
    return _atomic_savefig(fig, output_path, dpi=dpi)


def plot_regional(
    region_df: pd.DataFrame,
    top_variant: TopVariant,
    *,
    threshold: float,
    output_path: Path,
    title: str,
    dpi: int = DEFAULT_DPI,
) -> PlotResult:
    """Render a regional association plot; matches ``scripts/03_regional_plot.py``."""
    fig = plt.figure(figsize=(8, 4))
    ax = fig.gca()
    minuslog10p = -np.log10(region_df["pval"].to_numpy(dtype="float64"))
    ax.scatter(region_df["pos"], minuslog10p, s=8, alpha=0.7)
    ax.axhline(-np.log10(threshold), linestyle="--")
    ax.scatter([top_variant.pos], [-np.log10(top_variant.pval)], s=40, marker="*")

    ax.set_xlabel(f"Position on {top_variant.chrom} (bp)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    fig.tight_layout()
    return _atomic_savefig(fig, output_path, dpi=dpi)
