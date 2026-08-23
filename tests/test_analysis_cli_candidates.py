"""Unit tests for the ``candidates`` subcommand of adzuki_gwas_analysis.analysis.cli.

A separate file from tests/test_analysis_cli.py / tests/test_analysis_cli_diagnostics.py
(left unmodified by this Issue) so those pre-existing modules' diffs stay empty.
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.analysis.cli import main
from tests.analysis_support import write_dataset_file, write_manifest

_VALID_ROWS = [
    "Chr01\t.\t1000000\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-4\t1e-5",
    "Chr01\t.\t2000000\t0\tC\tT\t0.45\t-0.02\t0.02\t71.0\t31.0\t29.0\t5e-2\t4e-2\t3e-2",
    "Chr01\t.\t6000000\t0\tG\tA\t0.20\t0.03\t0.01\t70.5\t30.5\t28.5\t2e-3\t2e-4\t2e-5",
    "Chr02\t.\t1500000\t1\tT\tA\t0.10\t0.11\t0.03\t72.0\t32.0\t30.0\t9e-3\t8e-3\t7e-3",
]
_INVALID_ROWS = _VALID_ROWS + [
    "Chr02\t.\t9000000\t0\tG\tC\t0.5\t0.01\t0.01\t70.0\t30.0\t28.0\t1.5\t1.5\t1.5",
]


class CandidatesCliTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.data_dir = self.tmp_path / "data"
        self.data_dir.mkdir()
        self.manifest_path = self.tmp_path / "manifest.toml"
        self.output_dir = self.tmp_path / "out"

    def _write_valid_dataset(self) -> None:
        dataset_path = write_dataset_file(self.data_dir, _VALID_ROWS)
        write_manifest(self.manifest_path, dataset_path)

    def _write_invalid_dataset(self) -> None:
        dataset_path = write_dataset_file(self.data_dir, _INVALID_ROWS)
        write_manifest(self.manifest_path, dataset_path)

    def _run(self, args: list[str]) -> tuple[int, str]:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            code = main(args)
        return code, buf.getvalue()


class CandidatesCommandTests(CandidatesCliTestCase):
    def test_exit_code_zero_and_five_files_on_success(self) -> None:
        self._write_valid_dataset()
        code, output = self._run(
            [
                "candidates",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--dataset-id",
                "miyagi_water_permeability",
                "--output-dir",
                str(self.output_dir),
                "--clustering-distance",
                "1000",
            ]
        )
        self.assertEqual(code, 0)
        produced = {p.name for p in self.output_dir.iterdir()}
        self.assertEqual(
            produced,
            {
                "statistical_diagnostics.tsv",
                "significant_variants.tsv",
                "association_peaks.tsv",
                "candidate_snps.tsv",
                "candidate_ranking.tsv",
            },
        )
        self.assertIn("n_signals=", output)
        self.assertIn("not LD blocks or QTL intervals", output)

    def test_exit_code_nonzero_on_validation_failure_and_no_output(self) -> None:
        self._write_invalid_dataset()
        code, output = self._run(
            [
                "candidates",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--dataset-id",
                "miyagi_water_permeability",
                "--output-dir",
                str(self.output_dir),
                "--clustering-distance",
                "1000",
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertFalse(self.output_dir.exists())
        self.assertIn("ERROR", output)

    def test_missing_clustering_distance_is_a_usage_error(self) -> None:
        self._write_valid_dataset()
        with self.assertRaises(SystemExit):
            self._run(
                [
                    "candidates",
                    "--manifest",
                    str(self.manifest_path),
                    "--data-dir",
                    str(self.data_dir),
                    "--dataset-id",
                    "miyagi_water_permeability",
                    "--output-dir",
                    str(self.output_dir),
                ]
            )

    def test_missing_output_dir_is_a_usage_error(self) -> None:
        self._write_valid_dataset()
        with self.assertRaises(SystemExit):
            self._run(
                [
                    "candidates",
                    "--manifest",
                    str(self.manifest_path),
                    "--data-dir",
                    str(self.data_dir),
                    "--dataset-id",
                    "miyagi_water_permeability",
                    "--clustering-distance",
                    "1000",
                ]
            )

    def test_negative_clustering_distance_exits_nonzero_with_no_output(self) -> None:
        self._write_valid_dataset()
        code, output = self._run(
            [
                "candidates",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--dataset-id",
                "miyagi_water_permeability",
                "--output-dir",
                str(self.output_dir),
                "--clustering-distance",
                "-5",
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertFalse(self.output_dir.exists())


if __name__ == "__main__":
    unittest.main()
