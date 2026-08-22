"""Unit tests for adzuki_gwas_analysis.analysis.diagnostics.

Pure in-memory tests -- no file I/O.
"""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from adzuki_gwas_analysis.analysis.diagnostics import (
    DIAGNOSTICS_SCHEMA_VERSION,
    SIGNIFICANT_VARIANTS_COLUMNS,
    DiagnosticsResult,
    build_significant_variants_table,
    build_summary_table,
    format_family_scope,
)
from adzuki_gwas_analysis.analysis.statistics import (
    compute_bh,
    compute_bonferroni,
    compute_lambda_gc,
)


def _make_variant_df(pvalues: list[float]) -> pd.DataFrame:
    n = len(pvalues)
    return pd.DataFrame(
        {
            "chr": [f"Chr0{i % 2 + 1}" for i in range(n)],
            "pos": [1000 + i for i in range(n)],
            "allele1": ["A"] * n,
            "allele0": ["G"] * n,
            "af": [0.3] * n,
            "beta": [0.1 * i for i in range(n)],
            "pval": pvalues,
        }
    )


def _make_result(
    pvalues: list[float], *, alpha: float = 0.05, fdr_level: float = 0.05
) -> DiagnosticsResult:
    array = np.array(pvalues)
    return DiagnosticsResult(
        dataset_id="miyagi_water_permeability",
        source_sha256="a" * 64,
        pvalue_column="pval",
        n_tests=len(pvalues),
        bonferroni=compute_bonferroni(array, alpha=alpha),
        bh=compute_bh(array, fdr_level=fdr_level),
        lambda_gc=compute_lambda_gc(array, df=1),
    )


class FormatFamilyScopeTests(unittest.TestCase):
    def test_includes_dataset_id(self) -> None:
        scope = format_family_scope("miyagi_water_permeability")
        self.assertIn("miyagi_water_permeability", scope)

    def test_does_not_mention_other_datasets_or_p_columns(self) -> None:
        scope = format_family_scope("miyagi_water_permeability").lower()
        for forbidden in ("shumari", "p_wald", "p_score", "6 dataset"):
            self.assertNotIn(forbidden, scope)


class BuildSummaryTableTests(unittest.TestCase):
    def test_single_row(self) -> None:
        result = _make_result([0.9, 0.8, 0.7])
        table = build_summary_table(result)
        self.assertEqual(len(table), 1)

    def test_expected_columns_present_in_order(self) -> None:
        result = _make_result([0.9, 0.8, 0.7])
        table = build_summary_table(result)
        expected = [
            "schema_version",
            "dataset_id",
            "source_sha256",
            "pvalue_column",
            "pvalue_semantics",
            "family_scope",
            "n_tests",
            "alpha",
            "bonferroni_threshold",
            "bonferroni_discoveries",
            "fdr_level",
            "bh_raw_p_cutoff",
            "bh_discoveries",
            "lambda_gc_df",
            "expected_chi2_median",
            "lambda_gc",
        ]
        self.assertEqual(list(table.columns), expected)

    def test_schema_version_is_1(self) -> None:
        result = _make_result([0.9, 0.8, 0.7])
        table = build_summary_table(result)
        self.assertEqual(table["schema_version"].iloc[0], DIAGNOSTICS_SCHEMA_VERSION)

    def test_bh_raw_p_cutoff_is_nan_when_zero_discoveries(self) -> None:
        result = _make_result([0.9, 0.8, 0.7])
        self.assertIsNone(result.bh.raw_p_cutoff)
        table = build_summary_table(result)
        self.assertTrue(pd.isna(table["bh_raw_p_cutoff"].iloc[0]))

    def test_bh_raw_p_cutoff_is_set_when_discoveries_exist(self) -> None:
        result = _make_result([1e-10, 1e-12, 1e-15])
        table = build_summary_table(result)
        self.assertFalse(pd.isna(table["bh_raw_p_cutoff"].iloc[0]))

    def test_pvalue_semantics_names_lrt(self) -> None:
        result = _make_result([0.9, 0.8, 0.7])
        table = build_summary_table(result)
        self.assertIn("likelihood_ratio", table["pvalue_semantics"].iloc[0])

    def test_family_scope_matches_format_family_scope(self) -> None:
        result = _make_result([0.9, 0.8, 0.7])
        table = build_summary_table(result)
        self.assertEqual(
            table["family_scope"].iloc[0], format_family_scope("miyagi_water_permeability")
        )


class BuildSignificantVariantsTableTests(unittest.TestCase):
    def test_column_order(self) -> None:
        variant_df = _make_variant_df([1e-10, 0.9, 0.8])
        result = _make_result([1e-10, 0.9, 0.8])
        table = build_significant_variants_table(variant_df, result)
        self.assertEqual(list(table.columns), list(SIGNIFICANT_VARIANTS_COLUMNS))

    def test_zero_discoveries_yields_header_only_table(self) -> None:
        variant_df = _make_variant_df([0.9, 0.8, 0.7])
        result = _make_result([0.9, 0.8, 0.7])
        table = build_significant_variants_table(variant_df, result)
        self.assertEqual(len(table), 0)
        self.assertEqual(list(table.columns), list(SIGNIFICANT_VARIANTS_COLUMNS))

    def test_union_of_bonferroni_and_bh_significant(self) -> None:
        pvalues = [1e-10, 0.9, 0.8, 0.7, 0.6]
        variant_df = _make_variant_df(pvalues)
        result = _make_result(pvalues)
        table = build_significant_variants_table(variant_df, result)
        expected_rows = int((result.bonferroni.significant | result.bh.significant).sum())
        self.assertEqual(len(table), expected_rows)

    def test_preserves_original_input_row_order(self) -> None:
        # Two discoveries with positions that are NOT in ascending pval order --
        # the output must keep the original variant_df row order, not sort by
        # significance or by pval.
        pvalues = [0.9, 1e-12, 0.8, 1e-10, 0.7]
        variant_df = _make_variant_df(pvalues)
        result = _make_result(pvalues)
        table = build_significant_variants_table(variant_df, result)
        self.assertEqual(list(table["pos"]), [1001, 1003])

    def test_does_not_mutate_input_variant_df(self) -> None:
        variant_df = _make_variant_df([1e-10, 0.9, 0.8])
        snapshot = variant_df.copy()
        result = _make_result([1e-10, 0.9, 0.8])
        build_significant_variants_table(variant_df, result)
        pd.testing.assert_frame_equal(variant_df, snapshot)


if __name__ == "__main__":
    unittest.main()
