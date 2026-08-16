"""Unit tests for adzuki_gwas_analysis.analysis.qq.

Pure numeric tests -- no file I/O.
"""

from __future__ import annotations

import unittest

import numpy as np

from adzuki_gwas_analysis.analysis.qq import compute_qq_points


class ComputeQqPointsTests(unittest.TestCase):
    def test_observed_is_neg_log10_of_ascending_sorted_pvals(self) -> None:
        # p is sorted ascending, and -log10 is decreasing, so observed is
        # monotonically non-increasing (largest -log10(p) first).
        pvals = np.array([0.5, 1e-4, 0.01])
        points = compute_qq_points(pvals)
        self.assertEqual(points.variant_count, 3)
        np.testing.assert_allclose(points.observed, -np.log10(np.sort(pvals)))
        self.assertTrue(np.all(np.diff(points.observed) <= 0))

    def test_expected_matches_uniform_quantile_formula(self) -> None:
        pvals = np.array([0.9, 0.5, 0.1, 0.01])
        points = compute_qq_points(pvals)
        n = len(pvals)
        expected = -np.log10(np.arange(1, n + 1) / (n + 1))
        np.testing.assert_allclose(points.expected, expected)

    def test_rejects_pvalue_of_zero(self) -> None:
        with self.assertRaises(ValueError):
            compute_qq_points(np.array([0.0, 0.5]))

    def test_rejects_pvalue_above_one(self) -> None:
        with self.assertRaises(ValueError):
            compute_qq_points(np.array([1.5, 0.5]))

    def test_rejects_nan(self) -> None:
        with self.assertRaises(ValueError):
            compute_qq_points(np.array([np.nan, 0.5]))

    def test_accepts_boundary_pvalue_of_one(self) -> None:
        points = compute_qq_points(np.array([1.0, 0.5]))
        self.assertEqual(points.variant_count, 2)


if __name__ == "__main__":
    unittest.main()
