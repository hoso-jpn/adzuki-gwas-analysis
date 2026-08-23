"""Unit tests for the ``report`` subcommand of adzuki_gwas_analysis.analysis.cli.

A separate file from tests/test_analysis_cli.py / tests/test_analysis_cli_batch.py (left
otherwise unmodified by this Issue, aside from the pre-existing subcommand-count pin) so
those modules' diffs stay minimal.
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.analysis.cli import main
from tests.analysis_support import (
    build_candidate_enabled_batch_fixture,
    write_full_manifest,
    write_six_dataset_files,
)


class ReportCliTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.analysis_dir = build_candidate_enabled_batch_fixture(self.tmp_path)
        self.output_dir = self.tmp_path / "delivery"

    def _run(self, args: list[str]) -> tuple[int, str]:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            code = main(args)
        return code, buf.getvalue()


class ReportCommandTests(ReportCliTestCase):
    def test_exit_code_zero_and_delivery_package_on_success(self) -> None:
        code, output = self._run(
            [
                "report",
                "--analysis-dir",
                str(self.analysis_dir),
                "--output-dir",
                str(self.output_dir),
            ]
        )
        self.assertEqual(code, 0)
        self.assertTrue((self.output_dir / "executive_summary.md").is_file())
        self.assertTrue((self.output_dir / "analysis_report.md").is_file())
        self.assertTrue((self.output_dir / "reproducibility" / "run_manifest.json").is_file())
        self.assertIn("datasets=6", output)
        self.assertIn("not a shared", output)

    def test_missing_analysis_dir_flag_is_a_usage_error(self) -> None:
        with self.assertRaises(SystemExit):
            self._run(["report", "--output-dir", str(self.output_dir)])

    def test_missing_output_dir_flag_is_a_usage_error(self) -> None:
        with self.assertRaises(SystemExit):
            self._run(["report", "--analysis-dir", str(self.analysis_dir)])

    def test_schema_v1_analysis_dir_exits_nonzero_with_actionable_message(self) -> None:
        from adzuki_gwas_analysis.analysis import batch

        data_dir = self.tmp_path / "no_candidates_data"
        data_dir.mkdir()
        manifest_path = self.tmp_path / "no_candidates_manifest.toml"
        no_candidates_analysis_dir = self.tmp_path / "no_candidates_analysis_dir"
        rows = ["Chr01\t.\t1000000\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-3\t1e-3"]
        dataset_paths = write_six_dataset_files(
            data_dir,
            {
                did: rows
                for did in (
                    "miyagi_water_permeability",
                    "miyagi_red_seedcoat",
                    "miyagi_mottled_black_seedcoat",
                    "shumari_water_permeability",
                    "shumari_red_seedcoat",
                    "shumari_mottled_black_seedcoat",
                )
            },
        )
        write_full_manifest(manifest_path, dataset_paths)
        batch.run_batch(
            manifest_path=manifest_path, data_dir=data_dir, output_dir=no_candidates_analysis_dir
        )

        code, output = self._run(
            [
                "report",
                "--analysis-dir",
                str(no_candidates_analysis_dir),
                "--output-dir",
                str(self.output_dir),
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertFalse(self.output_dir.exists())
        self.assertIn("--clustering-distance", output)

    def test_nonexistent_analysis_dir_exits_nonzero_with_no_output(self) -> None:
        code, _output = self._run(
            [
                "report",
                "--analysis-dir",
                str(self.tmp_path / "does_not_exist"),
                "--output-dir",
                str(self.output_dir),
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertFalse(self.output_dir.exists())

    def test_existing_nonempty_output_dir_exits_nonzero(self) -> None:
        self.output_dir.mkdir(parents=True)
        (self.output_dir / "stray.txt").write_text("x", encoding="utf-8")
        code, output = self._run(
            [
                "report",
                "--analysis-dir",
                str(self.analysis_dir),
                "--output-dir",
                str(self.output_dir),
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertIn("not safe", output)
        self.assertTrue((self.output_dir / "stray.txt").is_file())


if __name__ == "__main__":
    unittest.main()
