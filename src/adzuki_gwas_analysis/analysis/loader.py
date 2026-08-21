"""Load one validated GWAS summary-statistics file into a DataFrame for analysis.

This is a distinct responsibility from
:mod:`adzuki_gwas_analysis.loader`, which streams rows one at a time purely
for validation and never returns a DataFrame. Analysis operations (Manhattan
layout, QQ quantiles, regional subsetting) need the position array and full
p-value column in memory at once, so some in-memory holding is unavoidable
here -- but only for the single dataset being analyzed, never more than one
of the 6 files at a time, and only for the 7 columns analysis actually uses
(not all 15 contract columns).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

#: The only columns analysis needs, out of the 15 required by schema v1.
ANALYSIS_COLUMNS: tuple[str, ...] = ("chr", "pos", "allele1", "allele0", "af", "beta", "pval")

#: Explicit dtypes, chosen to cut memory relative to pandas' default
#: per-column object/int64 inference without changing any value's meaning:
#: ``chr``/``allele1``/``allele0`` have very few distinct values across
#: ~1.7M rows (category avoids repeating full Python string objects per row);
#: ``pos`` fits comfortably in a 32-bit integer for this species' chromosome
#: lengths; ``af``/``beta``/``pval`` stay float64 so their string
#: representation round-trips identically to the original file (needed for
#: byte-for-byte TSV equivalence with pre-existing committed outputs).
ANALYSIS_DTYPES: dict[str, str] = {
    "chr": "category",
    "pos": "int32",
    "allele1": "category",
    "allele0": "category",
    "af": "float64",
    "beta": "float64",
    "pval": "float64",
}


def load_analysis_frame(path: Path) -> pd.DataFrame:
    """Read ``path``'s ``ANALYSIS_COLUMNS`` into a DataFrame with ``ANALYSIS_DTYPES``.

    Callers must run :func:`adzuki_gwas_analysis.validate.validate_dataset`
    against ``path`` first and check its ``success`` -- this function does
    not re-validate; it assumes the file already satisfies the schema v1
    contract (header, column count, value ranges) and simply parses it with
    pandas for numerical/plotting use.
    """
    return pd.read_csv(
        path,
        sep="\t",
        usecols=list(ANALYSIS_COLUMNS),
        dtype=ANALYSIS_DTYPES,
    )
