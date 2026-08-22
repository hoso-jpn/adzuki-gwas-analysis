"""Unit tests for the ``diagnostics`` subcommand of adzuki_gwas_analysis.analysis.cli.

A separate file from tests/test_analysis_cli.py (left unmodified by this
Issue) so that pre-existing module's diff stays empty.
"""

from __future__ import annotations

import contextlib
import importlib
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


class DiagnosticsCliTestCase(unittest.TestCase):
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


class DiagnosticsCommandTests(DiagnosticsCliTestCase):
    def test_exit_code_zero_and_two_files_on_success(self) -> None:
        self._write_valid_dataset()
        code, output = self._run(
            [
                "diagnostics",
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
        self.assertEqual(code, 0)
        produced = {p.name for p in self.output_dir.iterdir()}
        self.assertEqual(produced, {"statistical_diagnostics.tsv", "significant_variants.tsv"})
        self.assertIn("n_tests=4", output)

    def test_exit_code_nonzero_on_validation_failure_and_no_output(self) -> None:
        self._write_invalid_dataset()
        code, output = self._run(
            [
                "diagnostics",
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
        self.assertNotEqual(code, 0)
        self.assertIn("validation failed", output)
        self.assertFalse(self.output_dir.exists())

    def test_invalid_alpha_nonzero_exit(self) -> None:
        self._write_valid_dataset()
        code, output = self._run(
            [
                "diagnostics",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--dataset-id",
                "miyagi_water_permeability",
                "--output-dir",
                str(self.output_dir),
                "--alpha",
                "2.0",
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertIn("alpha", output)

    def test_invalid_fdr_level_nonzero_exit(self) -> None:
        self._write_valid_dataset()
        code, output = self._run(
            [
                "diagnostics",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--dataset-id",
                "miyagi_water_permeability",
                "--output-dir",
                str(self.output_dir),
                "--fdr-level",
                "0.0",
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertIn("fdr_level", output)

    def test_unknown_dataset_id_nonzero(self) -> None:
        self._write_valid_dataset()
        code, output = self._run(
            [
                "diagnostics",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--dataset-id",
                "not_a_real_dataset",
                "--output-dir",
                str(self.output_dir),
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertIn("unknown dataset_id", output)

    def test_output_parent_directory_created(self) -> None:
        self._write_valid_dataset()
        nested_output_dir = self.output_dir / "a" / "b"
        code, _ = self._run(
            [
                "diagnostics",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--dataset-id",
                "miyagi_water_permeability",
                "--output-dir",
                str(nested_output_dir),
            ]
        )
        self.assertEqual(code, 0)
        self.assertTrue((nested_output_dir / "statistical_diagnostics.tsv").is_file())

    def test_default_output_dir_is_plots_matching_other_subcommands(self) -> None:
        # Does not run the CLI (would write into ./plots); only checks the
        # parser default, matching --output-dir's default for manhattan/qq/etc.
        module = importlib.import_module("adzuki_gwas_analysis.analysis.cli")
        parser = module._build_parser()
        args = parser.parse_args(["diagnostics"])
        self.assertEqual(args.output_dir, module.DEFAULT_OUTPUT_DIR)
        self.assertEqual(args.alpha, module.DEFAULT_ALPHA)
        self.assertEqual(args.fdr_level, module.DEFAULT_FDR_LEVEL)
        self.assertFalse(hasattr(args, "threshold"))


if __name__ == "__main__":
    unittest.main()
