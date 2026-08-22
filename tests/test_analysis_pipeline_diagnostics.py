"""Unit tests for adzuki_gwas_analysis.analysis.pipeline.run_diagnostics.

A separate file from tests/test_analysis_pipeline.py (which this Issue leaves
unmodified) so the pre-existing test module's diff stays empty. Every test
uses tiny synthetic fixtures under a TemporaryDirectory -- none reads
data/raw/ (the real, un-tracked Dryad file).
"""

from __future__ import annotations

import dataclasses
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from adzuki_gwas_analysis.analysis import pipeline
from adzuki_gwas_analysis.errors import DatasetValidationFailedError, LoadedPvalueCountMismatchError
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


class DiagnosticsPipelineTestCase(unittest.TestCase):
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


class RunDiagnosticsTests(DiagnosticsPipelineTestCase):
    def test_produces_exactly_two_files(self) -> None:
        self._write_valid_dataset()
        pipeline.run_diagnostics(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
        )
        produced = {p.name for p in self.output_dir.iterdir()}
        self.assertEqual(produced, {"statistical_diagnostics.tsv", "significant_variants.tsv"})

    def test_n_tests_equals_row_count(self) -> None:
        self._write_valid_dataset()
        outcome = pipeline.run_diagnostics(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
        )
        self.assertEqual(outcome.result.n_tests, len(_VALID_ROWS))

    def test_pvalue_column_is_pval(self) -> None:
        self._write_valid_dataset()
        outcome = pipeline.run_diagnostics(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
        )
        self.assertEqual(outcome.result.pvalue_column, "pval")

    def test_produces_no_output_when_validation_fails(self) -> None:
        self._write_invalid_dataset()
        with self.assertRaises(DatasetValidationFailedError):
            pipeline.run_diagnostics(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_dir=self.output_dir,
            )
        self.assertFalse(self.output_dir.exists())

    def test_loader_not_called_when_validation_fails(self) -> None:
        self._write_invalid_dataset()
        with mock.patch(
            "adzuki_gwas_analysis.analysis.pipeline.load_analysis_frame"
        ) as mocked_load:
            with self.assertRaises(DatasetValidationFailedError):
                pipeline.run_diagnostics(
                    manifest_path=self.manifest_path,
                    data_dir=self.data_dir,
                    dataset_id="miyagi_water_permeability",
                    output_dir=self.output_dir,
                )
            mocked_load.assert_not_called()

    def test_zero_discoveries_yields_header_only_significant_variants_tsv(self) -> None:
        self._write_valid_dataset()
        outcome = pipeline.run_diagnostics(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
            alpha=1e-300,
            fdr_level=1e-300,
        )
        self.assertEqual(outcome.result.bonferroni.discoveries, 0)
        self.assertEqual(outcome.result.bh.discoveries, 0)
        table = pd.read_csv(outcome.significant_variants_path, sep="\t")
        self.assertEqual(len(table), 0)
        self.assertEqual(
            list(table.columns),
            [
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
            ],
        )

    def test_summary_tsv_has_one_row(self) -> None:
        self._write_valid_dataset()
        outcome = pipeline.run_diagnostics(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
        )
        table = pd.read_csv(outcome.summary_path, sep="\t")
        self.assertEqual(len(table), 1)
        self.assertEqual(table["dataset_id"].iloc[0], "miyagi_water_permeability")

    def test_output_parent_directory_created(self) -> None:
        self._write_valid_dataset()
        nested_output_dir = self.output_dir / "a" / "b"
        outcome = pipeline.run_diagnostics(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=nested_output_dir,
        )
        self.assertTrue(outcome.summary_path.is_file())
        self.assertTrue(outcome.significant_variants_path.is_file())

    def test_rejects_invalid_alpha(self) -> None:
        self._write_valid_dataset()
        with self.assertRaises(ValueError):
            pipeline.run_diagnostics(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_dir=self.output_dir,
                alpha=0.0,
            )

    def test_rejects_invalid_fdr_level(self) -> None:
        self._write_valid_dataset()
        with self.assertRaises(ValueError):
            pipeline.run_diagnostics(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_dir=self.output_dir,
                fdr_level=2.0,
            )

    def test_row_count_mismatch_raises_and_produces_no_output(self) -> None:
        self._write_valid_dataset()
        original_loader = pipeline.load_analysis_frame

        def _truncating_loader(path: Path) -> pd.DataFrame:
            return original_loader(path).iloc[:-1]

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.pipeline.load_analysis_frame",
                side_effect=_truncating_loader,
            ),
            self.assertRaises(LoadedPvalueCountMismatchError) as ctx,
        ):
            pipeline.run_diagnostics(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_dir=self.output_dir,
            )
        self.assertEqual(ctx.exception.validated_row_count, len(_VALID_ROWS))
        self.assertEqual(ctx.exception.loaded_count, len(_VALID_ROWS) - 1)
        self.assertFalse(self.output_dir.exists())

    def test_pvalue_column_comes_from_manifest_not_hardcoded(self) -> None:
        # Point pvalue_columns.primary at a different, already-loaded column
        # ("af") to prove run_diagnostics reads the column name from the
        # Manifest object rather than a literal "pval" string.
        self._write_valid_dataset()
        real_manifest = pipeline.load_manifest(self.manifest_path)
        patched_manifest = dataclasses.replace(
            real_manifest,
            pvalue_columns=dataclasses.replace(real_manifest.pvalue_columns, primary="af"),
        )
        with mock.patch(
            "adzuki_gwas_analysis.analysis.pipeline.load_manifest", return_value=patched_manifest
        ):
            outcome = pipeline.run_diagnostics(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_dir=self.output_dir,
            )
        self.assertEqual(outcome.result.pvalue_column, "af")

    def test_all_output_set_is_unaffected_by_diagnostics_module(self) -> None:
        # run_all must still produce exactly its existing 4-file set; adding
        # diagnostics must not change it.
        self._write_valid_dataset()
        from tests.analysis_support import write_region_config

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
        pipeline.run_all(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            regions_config_path=config_path,
            output_dir=self.output_dir,
        )
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


if __name__ == "__main__":
    unittest.main()
