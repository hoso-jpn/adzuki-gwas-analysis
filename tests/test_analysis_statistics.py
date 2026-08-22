"""Unit tests for adzuki_gwas_analysis.analysis.statistics.

Pure numeric tests -- no file I/O. Reference values for Benjamini-Hochberg are
computed by hand from the standard step-up formula (not by calling SciPy a
second time), and the lambda_GC reference values are well-known chi-square(1)
constants (the p=0.05 critical value 3.8414588... and the median
0.4549364...), so these tests do not merely re-verify a SciPy call with
another SciPy call.
"""

from __future__ import annotations

import unittest

import numpy as np

from adzuki_gwas_analysis.analysis.statistics import (
    compute_bh,
    compute_bonferroni,
    compute_lambda_gc,
)

# chi-square(1) critical value at p=0.05 and the chi-square(1) median --
# standard reference constants (e.g. Devlin & Roeder 1999's genomic-control
# convention), computed independently of this module for use as fixed
# expected values below.
_CHI2_1DF_CRITICAL_AT_P05 = 3.8414588206941285
_CHI2_1DF_MEDIAN = 0.454936423119572


class ValidatePvaluesTests(unittest.TestCase):
    """Shared numeric-contract behavior, exercised via compute_bonferroni."""

    def test_rejects_empty_array(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([]))

    def test_rejects_nan(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([0.1, np.nan]))

    def test_rejects_positive_infinity(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([0.1, np.inf]))

    def test_rejects_negative_infinity(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([0.1, -np.inf]))

    def test_rejects_zero(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([0.1, 0.0]))

    def test_rejects_negative(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([0.1, -0.01]))

    def test_rejects_greater_than_one(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([0.1, 1.5]))

    def test_accepts_boundary_pvalue_of_one(self) -> None:
        result = compute_bonferroni(np.array([1.0, 0.5]))
        self.assertEqual(result.m, 2)

    def test_accepts_single_element(self) -> None:
        result = compute_bonferroni(np.array([0.3]))
        self.assertEqual(result.m, 1)

    def test_does_not_mutate_input_array(self) -> None:
        original = np.array([0.3, 0.01, 0.5, 1.0])
        snapshot = original.copy()
        compute_bonferroni(original)
        compute_bh(original)
        compute_lambda_gc(original)
        np.testing.assert_array_equal(original, snapshot)


class ComputeBonferroniTests(unittest.TestCase):
    def test_threshold_is_alpha_over_m(self) -> None:
        pvalues = np.array([0.5, 0.01, 0.2, 0.001])
        result = compute_bonferroni(pvalues, alpha=0.05)
        self.assertAlmostEqual(result.threshold, 0.05 / 4)

    def test_adjusted_is_min_p_times_m_and_one(self) -> None:
        pvalues = np.array([0.5, 0.01, 0.9])
        result = compute_bonferroni(pvalues, alpha=0.05)
        np.testing.assert_allclose(result.adjusted, np.minimum(pvalues * 3, 1.0))
        self.assertTrue((result.adjusted <= 1.0).all())

    def test_significance_mask_matches_adjusted_p_le_alpha(self) -> None:
        # p <= alpha/m must be exactly equivalent to adjusted (min(p*m,1)) <= alpha.
        pvalues = np.array([0.001, 0.02, 0.5, 0.0125, 0.3])
        result = compute_bonferroni(pvalues, alpha=0.05)
        np.testing.assert_array_equal(result.significant, result.adjusted <= 0.05)

    def test_boundary_p_equal_to_threshold_is_significant(self) -> None:
        m = 4
        alpha = 0.05
        threshold = alpha / m
        pvalues = np.array([threshold, 0.9, 0.8, 0.7])
        result = compute_bonferroni(pvalues, alpha=alpha)
        self.assertTrue(bool(result.significant[0]))

    def test_zero_discoveries(self) -> None:
        pvalues = np.array([0.9, 0.8, 0.7])
        result = compute_bonferroni(pvalues, alpha=0.05)
        self.assertEqual(result.discoveries, 0)

    def test_all_discoveries(self) -> None:
        pvalues = np.array([1e-10, 1e-12, 1e-15])
        result = compute_bonferroni(pvalues, alpha=0.05)
        self.assertEqual(result.discoveries, 3)

    def test_rejects_invalid_alpha_zero(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([0.1, 0.2]), alpha=0.0)

    def test_rejects_invalid_alpha_above_one(self) -> None:
        with self.assertRaises(ValueError):
            compute_bonferroni(np.array([0.1, 0.2]), alpha=1.5)

    def test_m_equals_number_of_pvalues(self) -> None:
        pvalues = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        result = compute_bonferroni(pvalues)
        self.assertEqual(result.m, len(pvalues))


class ComputeBhTests(unittest.TestCase):
    def test_known_reference_vector_with_monotonicity_correction(self) -> None:
        # p_sorted ascending = [0.01, 0.015, 0.02, 0.021, 0.09], m=5.
        # raw*m/rank (ascending order) = [.05, .0375, .033333, .02625, .09] --
        # NOT monotone, so BH's cumulative-minimum-from-the-top step must
        # collapse ranks 1-4 to .02625 and leave rank 5 at .09. Computed by
        # hand from the standard step-up formula, not from a second SciPy call.
        # Fed in shuffled (non-ascending) order to also exercise order restoration.
        pvalues = np.array([0.09, 0.01, 0.021, 0.02, 0.015])
        expected_adjusted = np.array([0.09, 0.02625, 0.02625, 0.02625, 0.02625])
        result = compute_bh(pvalues, fdr_level=0.05)
        np.testing.assert_allclose(result.adjusted, expected_adjusted, rtol=1e-9)

    def test_tie_is_deterministic(self) -> None:
        # Duplicate p-values at consecutive ranks collapse to the same adjusted
        # value regardless of which tied element is assigned which rank:
        # sorted=[0.02,0.02,0.04,0.04], raw*m/rank=[.08,.04,.05333,.04];
        # cumulative min from the top -> [.04,.04,.04,.04] for every position.
        pvalues = np.array([0.02, 0.04, 0.02, 0.04])
        result = compute_bh(pvalues)
        np.testing.assert_allclose(result.adjusted, np.full(4, 0.04), rtol=1e-9)

    def test_adjusted_p_le_fdr_level_matches_significant_mask(self) -> None:
        pvalues = np.array([0.001, 0.2, 0.03, 0.5, 0.0005])
        result = compute_bh(pvalues, fdr_level=0.05)
        np.testing.assert_array_equal(result.significant, result.adjusted <= 0.05)

    def test_zero_discoveries_yields_none_cutoff(self) -> None:
        pvalues = np.array([0.9, 0.8, 0.7])
        result = compute_bh(pvalues, fdr_level=0.05)
        self.assertEqual(result.discoveries, 0)
        self.assertIsNone(result.raw_p_cutoff)

    def test_all_discoveries(self) -> None:
        pvalues = np.array([1e-10, 1e-12, 1e-15])
        result = compute_bh(pvalues, fdr_level=0.05)
        self.assertEqual(result.discoveries, 3)
        self.assertIsNotNone(result.raw_p_cutoff)

    def test_raw_p_cutoff_is_max_raw_p_among_rejected(self) -> None:
        pvalues = np.array([0.001, 0.01, 0.5, 0.02])
        result = compute_bh(pvalues, fdr_level=0.05)
        if result.discoveries > 0:
            rejected_raw = pvalues[result.significant]
            self.assertAlmostEqual(result.raw_p_cutoff, float(rejected_raw.max()))

    def test_adjusted_capped_at_one(self) -> None:
        pvalues = np.array([0.9, 0.95, 1.0])
        result = compute_bh(pvalues)
        self.assertTrue((result.adjusted <= 1.0).all())

    def test_restores_original_input_order(self) -> None:
        ascending = np.array([0.001, 0.01, 0.02, 0.03, 0.5])
        shuffled = ascending[::-1].copy()
        result_ascending = compute_bh(ascending)
        result_shuffled = compute_bh(shuffled)
        np.testing.assert_allclose(
            result_shuffled.adjusted, result_ascending.adjusted[::-1], rtol=1e-9
        )

    def test_rejects_invalid_fdr_level(self) -> None:
        with self.assertRaises(ValueError):
            compute_bh(np.array([0.1, 0.2]), fdr_level=0.0)


class ComputeLambdaGcTests(unittest.TestCase):
    def test_matches_known_chi2_1df_reference_values(self) -> None:
        # median p-value is exactly 0.05 -> chi2 median must equal the
        # well-known chi-square(1) critical value at p=0.05 (~3.84146),
        # not -2*log(0.05) (~5.99146), which would indicate an incorrect
        # -2*log(p) substitution for the survival-function inverse.
        pvalues = np.array([0.01, 0.05, 0.5])
        result = compute_lambda_gc(pvalues, df=1)
        self.assertAlmostEqual(result.expected_median, _CHI2_1DF_MEDIAN, places=6)
        expected_lambda_gc = _CHI2_1DF_CRITICAL_AT_P05 / _CHI2_1DF_MEDIAN
        self.assertAlmostEqual(result.lambda_gc, expected_lambda_gc, places=6)

    def test_p_equal_one_maps_to_chi2_zero(self) -> None:
        pvalues = np.array([1.0, 1.0, 1.0])
        result = compute_lambda_gc(pvalues, df=1)
        self.assertAlmostEqual(result.lambda_gc, 0.0)

    def test_does_not_use_minus_two_log_p_substitute(self) -> None:
        # If the implementation used -2*log(p) instead of chi2.isf, lambda_gc
        # for this all-p=0.05 vector would come out close to
        # (-2*log(0.05)) / chi2_median, not 1.0. The correct chi2.isf(0.05,1)
        # equals the chi2 median's *scaling factor* such that lambda_gc == 1
        # exactly when every p equals chi2.ppf(0.5, df) inverted back to a
        # p-value -- more directly: for p=0.05 repeated, lambda_gc must equal
        # the known ratio below, not the -2log(p)-based ratio.
        pvalues = np.full(5, 0.05)
        result = compute_lambda_gc(pvalues, df=1)
        correct_ratio = _CHI2_1DF_CRITICAL_AT_P05 / _CHI2_1DF_MEDIAN
        wrong_ratio = (-2.0 * np.log(0.05)) / _CHI2_1DF_MEDIAN
        self.assertAlmostEqual(result.lambda_gc, correct_ratio, places=6)
        self.assertGreater(abs(result.lambda_gc - wrong_ratio), 0.5)

    def test_expected_median_not_hardcoded_for_other_df(self) -> None:
        pvalues = np.array([0.1, 0.5, 0.9])
        result_df1 = compute_lambda_gc(pvalues, df=1)
        result_df2 = compute_lambda_gc(pvalues, df=2)
        self.assertNotAlmostEqual(result_df1.expected_median, result_df2.expected_median)

    def test_rejects_non_positive_df(self) -> None:
        with self.assertRaises(ValueError):
            compute_lambda_gc(np.array([0.1, 0.5]), df=0)

    def test_result_is_finite(self) -> None:
        pvalues = np.array([0.001, 0.5, 0.99])
        result = compute_lambda_gc(pvalues, df=1)
        self.assertTrue(np.isfinite(result.lambda_gc))

    def test_independent_of_input_order(self) -> None:
        pvalues = np.array([0.01, 0.5, 0.2, 0.05, 0.8])
        forward = compute_lambda_gc(pvalues, df=1)
        shuffled = compute_lambda_gc(pvalues[::-1].copy(), df=1)
        self.assertAlmostEqual(forward.lambda_gc, shuffled.lambda_gc)


if __name__ == "__main__":
    unittest.main()
