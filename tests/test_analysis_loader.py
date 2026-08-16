"""Unit tests for adzuki_gwas_analysis.analysis.loader."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.analysis.loader import ANALYSIS_COLUMNS, load_analysis_frame
from tests.analysis_support import write_dataset_file

_ROWS = [
    "Chr01\t.\t1000\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-4\t1e-5",
    "Chr01\t.\t2000\t0\tC\tT\t0.45\t-0.02\t0.02\t71.0\t31.0\t29.0\t5e-2\t4e-2\t3e-2",
    "Chr02\t.\t1500\t1\tT\tA\t0.10\t0.11\t0.03\t72.0\t32.0\t30.0\t9e-3\t8e-3\t7e-3",
]


class LoadAnalysisFrameTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.data_dir = Path(self._tmpdir.name)
        self.dataset_path = write_dataset_file(self.data_dir, _ROWS)

    def test_reads_only_the_expected_columns(self) -> None:
        df = load_analysis_frame(self.dataset_path)
        self.assertEqual(tuple(df.columns), ANALYSIS_COLUMNS)

    def test_dtypes_match_contract(self) -> None:
        df = load_analysis_frame(self.dataset_path)
        self.assertEqual(str(df["chr"].dtype), "category")
        self.assertEqual(str(df["allele1"].dtype), "category")
        self.assertEqual(str(df["allele0"].dtype), "category")
        self.assertEqual(str(df["pos"].dtype), "int32")
        self.assertEqual(str(df["af"].dtype), "float64")
        self.assertEqual(str(df["beta"].dtype), "float64")
        self.assertEqual(str(df["pval"].dtype), "float64")

    def test_row_count_and_values(self) -> None:
        df = load_analysis_frame(self.dataset_path)
        self.assertEqual(len(df), 3)
        self.assertEqual(list(df["pos"]), [1000, 2000, 1500])
        self.assertAlmostEqual(df["pval"].iloc[0], 1e-4)


if __name__ == "__main__":
    unittest.main()
