"""Result dataclass and output-table construction for statistical diagnostics.

The multiple-testing family computed here is fixed and deliberately narrow:
**one dataset_id x the manifest's declared primary p-value column x every
schema-v1-validated variant in that one file.** The 6 Dryad files, the two
references (Miyagi/Shumari), the 3 traits, and post-hoc visualization
regions are never combined into one family, and ``p_wald``/``p_score`` are
never substituted for the primary column -- see
:mod:`adzuki_gwas_analysis.analysis.statistics` for the numeric methods and
their own caveats (FWER vs. FDR, the dependency condition BH actually
requires, and what lambda_GC does and does not establish).
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from adzuki_gwas_analysis.analysis.statistics import BhResult, BonferroniResult, LambdaGcResult

#: Bumped only if the summary TSV's column set/meaning changes.
DIAGNOSTICS_SCHEMA_VERSION = 1

#: What ``pval`` (the manifest's primary column) means, recorded in machine-readable
#: form so the summary TSV is self-describing without cross-referencing the manifest.
PVALUE_SEMANTICS = "likelihood_ratio_test"

#: Column order for significant_variants.tsv (Phase 0 / Issue #7 recommendation).
SIGNIFICANT_VARIANTS_COLUMNS: tuple[str, ...] = (
    "chr",
    "pos",
    "allele1",
    "allele0",
    "af",
    "beta",
    "pval",
    "pval_bonferroni",
    "pval_bh",
    "bonferroni_significant",
    "bh_significant",
)


def format_family_scope(dataset_id: str) -> str:
    """Machine-readable statement of the multiple-testing family for ``dataset_id``."""
    return (
        f"1 dataset_id ({dataset_id}) x manifest pvalue_columns.primary "
        f"x schema-v1-validated variants in that one file"
    )


@dataclass(frozen=True, slots=True)
class DiagnosticsResult:
    """Everything needed to render the two diagnostics output files for one dataset."""

    dataset_id: str
    source_sha256: str
    pvalue_column: str
    n_tests: int
    bonferroni: BonferroniResult
    bh: BhResult
    lambda_gc: LambdaGcResult


def build_summary_table(result: DiagnosticsResult) -> pd.DataFrame:
    """Build the 1-row ``statistical_diagnostics.tsv`` table.

    ``bh_raw_p_cutoff`` is written as an explicit missing value (empty field
    in the TSV) when there were zero BH discoveries, via ``float("nan")`` in
    an otherwise-float64 column -- never a Python ``None`` sentinel, which
    pandas would otherwise infer as an ``object``-dtype column.
    """
    raw_p_cutoff = result.bh.raw_p_cutoff
    row = {
        "schema_version": DIAGNOSTICS_SCHEMA_VERSION,
        "dataset_id": result.dataset_id,
        "source_sha256": result.source_sha256,
        "pvalue_column": result.pvalue_column,
        "pvalue_semantics": PVALUE_SEMANTICS,
        "family_scope": format_family_scope(result.dataset_id),
        "n_tests": result.n_tests,
        "alpha": result.bonferroni.alpha,
        "bonferroni_threshold": result.bonferroni.threshold,
        "bonferroni_discoveries": result.bonferroni.discoveries,
        "fdr_level": result.bh.fdr_level,
        "bh_raw_p_cutoff": raw_p_cutoff if raw_p_cutoff is not None else float("nan"),
        "bh_discoveries": result.bh.discoveries,
        "lambda_gc_df": result.lambda_gc.df,
        "expected_chi2_median": result.lambda_gc.expected_median,
        "lambda_gc": result.lambda_gc.lambda_gc,
    }
    return pd.DataFrame([row])


def build_significant_variants_table(
    variant_df: pd.DataFrame, result: DiagnosticsResult
) -> pd.DataFrame:
    """Build the ``significant_variants.tsv`` table: the union of Bonferroni/BH discoveries.

    ``variant_df`` must be the same, in-order DataFrame (``chr, pos, allele1,
    allele0, af, beta, pval``) that ``result.bonferroni``/``result.bh`` were
    computed from -- row order is preserved (this repository's existing
    ``subset_region``/``select_top_variant`` deliberately never sort either),
    and rows where neither correction flags significance are dropped. Yields
    a header-only, zero-row table (not an omitted file) when nothing is
    significant.
    """
    out = variant_df.loc[:, ["chr", "pos", "allele1", "allele0", "af", "beta", "pval"]].copy()
    out["pval_bonferroni"] = result.bonferroni.adjusted
    out["pval_bh"] = result.bh.adjusted
    out["bonferroni_significant"] = result.bonferroni.significant
    out["bh_significant"] = result.bh.significant

    mask = result.bonferroni.significant | result.bh.significant
    filtered: pd.DataFrame = out.loc[mask, list(SIGNIFICANT_VARIANTS_COLUMNS)]
    return filtered.reset_index(drop=True)
