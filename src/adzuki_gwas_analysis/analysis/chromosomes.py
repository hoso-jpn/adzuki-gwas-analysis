"""Chromosome ordering and Manhattan-plot x-coordinate computation.

The original ``scripts/01_manhattan_plot.py`` ordered chromosomes with plain
Python ``sorted(df["chr"].unique())`` -- a lexicographic string sort. That
happens to match numeric order for the real data (``Chr01``..``Chr11``, all
zero-padded to 2 digits), but it is a coincidence of that padding, not a
guarantee: lexicographic order would put ``"Chr10"`` before ``"Chr2"``.
:func:`natural_chromosome_key` sorts by the leading non-digit prefix and the
trailing integer instead, so the order is correct regardless of padding. This
module makes no claim about a reference genome's own contig order (e.g. which
chromosome comes "first" biologically) -- only that ``ChrN`` sorts by
numeric ``N``.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

_TRAILING_DIGITS_RE = re.compile(r"^(.*?)(\d+)$")


def natural_chromosome_key(chrom: str) -> tuple[int, str, int, str]:
    """Sort key for chromosome names, numeric-suffix-aware.

    Names matching ``<prefix><digits>`` (e.g. ``"Chr7"``, ``"Chr07"``) sort by
    ``(0, prefix, int(digits), "")``, so ``"Chr2"`` sorts before ``"Chr10"``.
    Names with no trailing digits sort after all of those, lexicographically
    by the full string, so the key is still total and deterministic.
    """
    match = _TRAILING_DIGITS_RE.match(chrom)
    if match is None:
        return (1, "", 0, chrom)
    prefix, digits = match.groups()
    return (0, prefix, int(digits), "")


def order_chromosomes(chroms: Sequence[str]) -> list[str]:
    """Return the distinct values of ``chroms`` in natural chromosome order."""
    return sorted(set(chroms), key=natural_chromosome_key)


@dataclass(frozen=True, slots=True)
class ManhattanCoordinates:
    """Genome-wide x-coordinates for a Manhattan plot, plus axis tick metadata."""

    x: pd.Series
    tick_positions: list[float]
    tick_labels: list[str]
    chrom_order: list[str]
    variant_count: int
    chrom_variant_counts: dict[str, int]


def compute_manhattan_coordinates(df: pd.DataFrame) -> ManhattanCoordinates:
    """Compute genome-wide x-coordinates for ``df`` (must have ``chr``, ``pos``).

    Chromosomes are laid out left-to-right in :func:`natural_chromosome_key`
    order. Each chromosome's x-coordinate is its ``pos`` plus a running
    offset equal to the sum of previous chromosomes' max ``pos`` -- the same
    accumulation used by the original ``scripts/01_manhattan_plot.py``, which
    (for this repository's zero-padded ``Chr01``..``Chr11`` data) natural
    order reproduces exactly.
    """
    chrom_order = order_chromosomes(df["chr"].astype(str).unique().tolist())

    x = pd.Series(index=df.index, dtype="float64")
    tick_positions: list[float] = []
    tick_labels: list[str] = []
    chrom_variant_counts: dict[str, int] = {}

    offset = 0
    for chrom in chrom_order:
        mask = df["chr"] == chrom
        sub_pos = df.loc[mask, "pos"]
        chrom_x = sub_pos + offset
        x.loc[mask] = chrom_x

        tick_positions.append(float(chrom_x.median()))
        tick_labels.append(chrom.replace("Chr", ""))
        chrom_variant_counts[chrom] = int(mask.sum())

        offset += int(sub_pos.max())

    return ManhattanCoordinates(
        x=x,
        tick_positions=tick_positions,
        tick_labels=tick_labels,
        chrom_order=chrom_order,
        variant_count=len(df),
        chrom_variant_counts=chrom_variant_counts,
    )
