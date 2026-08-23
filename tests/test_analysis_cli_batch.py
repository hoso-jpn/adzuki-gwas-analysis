"""Unit tests for the ``batch`` subcommand of adzuki_gwas_analysis.analysis.cli.

A separate file from tests/test_analysis_cli.py and tests/test_analysis_cli_diagnostics.py
(both left unmodified by this Issue) so those pre-existing modules' diffs stay empty.
"""

from __future__ import annotations

import contextlib
import importlib
import io
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.analysis.cli import main
from tests.analysis_support import write_full_manifest, write_six_dataset_files

_DATASET_IDS_IN_MANIFEST_ORDER: tuple[str, ...] = (
    "miyagi_water_permeability",
    "miyagi_red_seedcoat",
    "miyagi_mottled_black_seedcoat",
    "shumari_water_permeability",
    "shumari_red_seedcoat",
    "shumari_mottled_black_seedcoat",
)

_ROW_TEMPLATE = (
    "Chr01\t.\t{pos1}\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t{p1}\t1e-3",
    "Chr01\t.\t{pos2}\t0\tC\tT\t0.45\t-0.02\t0.02\t71.0\t31.0\t29.0\t5e-2\t{p2}\t3e-2",
    "Chr02\t.\t{pos3}\t1\tT\tA\t0.10\t0.11\t0.03\t72.0\t32.0\t30.0\t9e-3\t{p3}\t7e-3",
)


def _six_dataset_rows() -> dict[str, list[str]]:
    rows_by_id: dict[str, list[str]] = {}
    for i, dataset_id in enumerate(_DATASET_IDS_IN_MANIFEST_ORDER):
        rows_by_id[dataset_id] = [
            row.format(
                pos1=1_000_000 + i,
                pos2=2_000_000 + i,
                pos3=1_500_000 + i,
                p1="2e-2",
                p2="4e-2",
                p3="8e-3",
            )
            for row in _ROW_TEMPLATE
        ]
    return rows_by_id


class BatchCliTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.data_dir = self.tmp_path / "data"
        self.data_dir.mkdir()
        self.manifest_path = self.tmp_path / "manifest.toml"
        self.output_dir = self.tmp_path / "out"

    def _write_all_six_valid(self) -> None:
        dataset_paths = write_six_dataset_files(self.data_dir, _six_dataset_rows())
        write_full_manifest(self.manifest_path, dataset_paths)

    def _run(self, args: list[str]) -> tuple[int, str]:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            code = main(args)
        return code, buf.getvalue()


class BatchCommandTests(BatchCliTestCase):
    def test_exit_code_zero_and_twenty_five_files_on_success(self) -> None:
        self._write_all_six_valid()
        code, output = self._run(
            [
                "batch",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--output-dir",
                str(self.output_dir),
            ]
        )
        self.assertEqual(code, 0)
        all_files = [p for p in self.output_dir.rglob("*") if p.is_file()]
        self.assertEqual(len(all_files), 25)
        for dataset_id in _DATASET_IDS_IN_MANIFEST_ORDER:
            self.assertIn(dataset_id, output)

    def test_output_dir_is_required(self) -> None:
        module = importlib.import_module("adzuki_gwas_analysis.analysis.cli")
        parser = module._build_parser()
        buf = io.StringIO()
        with contextlib.redirect_stderr(buf), self.assertRaises(SystemExit) as ctx:
            parser.parse_args(["batch"])
        self.assertEqual(ctx.exception.code, 2)
        self.assertIn("--output-dir", buf.getvalue())

    def test_no_dataset_id_argument_exists(self) -> None:
        module = importlib.import_module("adzuki_gwas_analysis.analysis.cli")
        parser = module._build_parser()
        args = parser.parse_args(["batch", "--output-dir", "/tmp/whatever-batch-dir"])
        self.assertFalse(hasattr(args, "dataset_id"))

    def test_default_alpha_fdr_level_threshold(self) -> None:
        module = importlib.import_module("adzuki_gwas_analysis.analysis.cli")
        parser = module._build_parser()
        args = parser.parse_args(["batch", "--output-dir", "/tmp/whatever-batch-dir-2"])
        self.assertEqual(args.alpha, module.DEFAULT_ALPHA)
        self.assertEqual(args.fdr_level, module.DEFAULT_FDR_LEVEL)
        self.assertEqual(args.threshold, module.DEFAULT_THRESHOLD)

    def test_exit_code_nonzero_on_second_dataset_validation_failure(self) -> None:
        self._write_all_six_valid()
        second_dataset_path = self.data_dir / "mapped_to_Miyagi_red_seedcoat.maf_0.05.assoc.txt"
        with second_dataset_path.open("a", encoding="utf-8") as fh:
            fh.write(
                "Chr01\t.\t9999999\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-3\t1e-3\n"
            )
        code, _ = self._run(
            [
                "batch",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--output-dir",
                str(self.output_dir),
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertFalse(self.output_dir.exists())

    def test_rejects_nonempty_output_dir_with_nonzero_exit(self) -> None:
        self._write_all_six_valid()
        self.output_dir.mkdir(parents=True)
        (self.output_dir / "stray.txt").write_text("pre-existing", encoding="utf-8")
        code, output = self._run(
            [
                "batch",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--output-dir",
                str(self.output_dir),
            ]
        )
        self.assertNotEqual(code, 0)
        self.assertIn("not safe", output)
        self.assertTrue((self.output_dir / "stray.txt").is_file())


class ExistingSevenSubcommandsUnaffectedTests(BatchCliTestCase):
    """Pins that adding `batch` (and later `candidates`) did not change any prior subcommand."""

    def test_all_nine_subcommands_are_registered(self) -> None:
        module = importlib.import_module("adzuki_gwas_analysis.analysis.cli")
        parser = module._build_parser()
        subparsers_action = next(
            action
            for action in parser._subparsers._group_actions  # type: ignore[union-attr]
            if hasattr(action, "choices")
        )
        self.assertEqual(
            set(subparsers_action.choices),
            {
                "manhattan",
                "qq",
                "regional",
                "regions",
                "top-variants",
                "all",
                "diagnostics",
                "batch",
                "candidates",
            },
        )

    def test_existing_subcommands_output_dir_and_output_contract_unchanged(self) -> None:
        module = importlib.import_module("adzuki_gwas_analysis.analysis.cli")
        parser = module._build_parser()

        for command in ("manhattan", "qq", "regions", "all"):
            args = parser.parse_args([command])
            self.assertEqual(args.output_dir, module.DEFAULT_OUTPUT_DIR)
            self.assertEqual(args.dataset_id, module.DEFAULT_DATASET_ID)

        regional_args = parser.parse_args(
            ["regional", "--chrom", "Chr01", "--start", "1", "--end", "2", "--output", "x.png"]
        )
        self.assertEqual(regional_args.dataset_id, module.DEFAULT_DATASET_ID)

        top_variants_args = parser.parse_args(["top-variants", "--output", "x.tsv"])
        self.assertEqual(top_variants_args.dataset_id, module.DEFAULT_DATASET_ID)

        diagnostics_args = parser.parse_args(["diagnostics", "--output-dir", "/tmp/x"])
        self.assertEqual(diagnostics_args.dataset_id, module.DEFAULT_DATASET_ID)
        self.assertEqual(diagnostics_args.alpha, module.DEFAULT_ALPHA)
        self.assertEqual(diagnostics_args.fdr_level, module.DEFAULT_FDR_LEVEL)

    def test_diagnostics_output_dir_still_required(self) -> None:
        module = importlib.import_module("adzuki_gwas_analysis.analysis.cli")
        parser = module._build_parser()
        buf = io.StringIO()
        with contextlib.redirect_stderr(buf), self.assertRaises(SystemExit):
            parser.parse_args(["diagnostics"])

    def test_all_command_output_set_unchanged_by_presence_of_batch(self) -> None:
        from tests.analysis_support import write_dataset_file, write_manifest, write_region_config

        valid_rows = [
            "Chr01\t.\t1000000\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-4\t1e-5",
        ]
        dataset_path = write_dataset_file(self.data_dir, valid_rows)
        write_manifest(self.manifest_path, dataset_path)
        one_region = [
            {
                "region_id": "Chr01_1_3Mb",
                "chrom": "Chr01",
                "start": 1_000_000,
                "end": 3_000_000,
                "title": "Region 1",
                "output_filename": "region1.png",
            }
        ]
        config_path = write_region_config(self.tmp_path / "regions.toml", one_region)
        code, _ = self._run(
            [
                "all",
                "--manifest",
                str(self.manifest_path),
                "--data-dir",
                str(self.data_dir),
                "--dataset-id",
                "miyagi_water_permeability",
                "--output-dir",
                str(self.output_dir),
                "--regions-config",
                str(config_path),
            ]
        )
        self.assertEqual(code, 0)
        produced = {p.name for p in self.output_dir.iterdir()}
        self.assertEqual(
            produced,
            {
                "miyagi_water_permeability_manhattan.png",
                "miyagi_water_permeability_qq.png",
                "region1.png",
                "top_variants_by_region.tsv",
            },
        )

    def test_diagnostics_command_output_set_unchanged_by_presence_of_batch(self) -> None:
        from tests.analysis_support import write_dataset_file, write_manifest

        valid_rows = [
            "Chr01\t.\t1000000\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-4\t1e-5",
        ]
        dataset_path = write_dataset_file(self.data_dir, valid_rows)
        write_manifest(self.manifest_path, dataset_path)
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
                str(self.output_dir),
            ]
        )
        self.assertEqual(code, 0)
        produced = {p.name for p in self.output_dir.iterdir()}
        self.assertEqual(produced, {"statistical_diagnostics.tsv", "significant_variants.tsv"})


if __name__ == "__main__":
    unittest.main()
