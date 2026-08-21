"""Regional subsetting and top-variant selection within a post-hoc visualization window.

These windows (see :mod:`adzuki_gwas_analysis.analysis.region_config`) are
visualization-only: they are not independently defined QTL intervals or LD
blocks, and this module makes no such claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import pandas as pd

from adzuki_gwas_analysis.errors import EmptyRegionError


def subset_region(df: pd.DataFrame, *, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Return rows of ``df`` with ``chr == chrom`` and ``start <= pos <= end``.

    Row order is preserved from ``df`` (no sorting), which matters for
    :func:`select_top_variant`'s first-occurrence tie-break.
    """
    return df[(df["chr"] == chrom) & (df["pos"] >= start) & (df["pos"] <= end)]


@dataclass(frozen=True, slots=True)
class TopVariant:
    """The single most-significant variant in a region, by ``pval``."""

    chrom: str
    pos: int
    allele1: str
    allele0: str
    af: float
    beta: float
    pval: float


def select_top_variant(
    df: pd.DataFrame, *, dataset_id: str, region_id: str, chrom: str, start: int, end: int
) -> TopVariant:
    """Return the row with the smallest ``pval`` in ``df``.

    On a tie, ``pandas.Series.idxmin`` returns the *first* occurrence in
    ``df``'s row order -- i.e. the first occurrence in the source file, not
    (for example) the smallest ``pos`` -- matching the original
    ``scripts/03_regional_plot.py`` / ``04_extract_top_variants_by_region.py``
    behavior exactly.

    Raises :class:`~adzuki_gwas_analysis.errors.EmptyRegionError` if ``df`` is
    empty, rather than raising an unlabeled ``ValueError`` from ``idxmin``.
    """
    if df.empty:
        raise EmptyRegionError(
            dataset_id=dataset_id, region_id=region_id, chrom=chrom, start=start, end=end
        )

    idx = df["pval"].idxmin()
    return TopVariant(
        chrom=str(df.at[idx, "chr"]),
        pos=int(cast(Any, df.at[idx, "pos"])),
        allele1=str(df.at[idx, "allele1"]),
        allele0=str(df.at[idx, "allele0"]),
        af=float(cast(Any, df.at[idx, "af"])),
        beta=float(cast(Any, df.at[idx, "beta"])),
        pval=float(cast(Any, df.at[idx, "pval"])),
    )
