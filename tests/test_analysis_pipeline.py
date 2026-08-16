"""Unit tests for adzuki_gwas_analysis.analysis.pipeline.

Every test uses tiny synthetic fixtures under a ``TemporaryDirectory`` --
none reads ``data/raw/`` (the real, un-tracked Dryad files).
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from adzuki_gwas_analysis.analysis import pipeline
from adzuki_gwas_analysis.errors import (
    DatasetValidationFailedError,
    EmptyRegionError,
    UnknownDatasetIdError,
)
from tests.analysis_support import write_dataset_file, write_manifest, write_region_config

_VALID_ROWS = [
    "Chr01\t.\t1000000\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-4\t1e-5",
    "Chr01\t.\t2000000\t0\tC\tT\t0.45\t-0.02\t0.02\t71.0\t31.0\t29.0\t5e-2\t4e-2\t3e-2",
    "Chr01\t.\t6000000\t0\tG\tA\t0.20\t0.03\t0.01\t70.5\t30.5\t28.5\t2e-3\t2e-4\t2e-5",
    "Chr02\t.\t1500000\t1\tT\tA\t0.10\t0.11\t0.03\t72.0\t32.0\t30.0\t9e-3\t8e-3\t7e-3",
]

# Row with pval > 1: violates the schema v1 contract.
_INVALID_ROWS = _VALID_ROWS + [
    "Chr02\t.\t9000000\t0\tG\tC\t0.5\t0.01\t0.01\t70.0\t30.0\t28.0\t1.5\t1.5\t1.5",
]

_ONE_REGION = [
    {
        "region_id": "Chr01_1_3Mb",
        "chrom": "Chr01",
        "start": 1_000_000,
        "end": 3_000_000,
        "title": "Region 1",
        "output_filename": "region1.png",
    }
]

_EMPTY_REGION = [
    {
        "region_id": "Chr01_empty",
        "chrom": "Chr01",
        "start": 9_000_000,
        "end": 9_100_000,
        "title": "Empty region",
        "output_filename": "empty.png",
    }
]


class PipelineTestCase(unittest.TestCase):
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


class EnsureValidatedTests(PipelineTestCase):
    def test_succeeds_on_valid_dataset(self) -> None:
        self._write_valid_dataset()
        entry = pipeline.ensure_validated(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
        )
        self.assertEqual(entry.dataset_id, "miyagi_water_permeability")

    def test_raises_on_invalid_dataset(self) -> None:
        self._write_invalid_dataset()
        with self.assertRaises(DatasetValidationFailedError):
            pipeline.ensure_validated(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
            )

    def test_unknown_dataset_id_raises(self) -> None:
        self._write_valid_dataset()
        with self.assertRaises(UnknownDatasetIdError):
            pipeline.ensure_validated(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="not_a_real_dataset",
            )

    def test_never_touches_other_datasets_files(self) -> None:
        # The manifest declares 5 other datasets whose files do not exist on
        # disk at all. Validating miyagi_water_permeability must never
        # attempt to open them.
        self._write_valid_dataset()
        other_files = [
            "mapped_to_Miyagi_red_seedcoat.maf_0.05.assoc.txt",
            "mapped_to_Miyagi_mottled_black_seedcoat.maf_0.05.assoc.txt",
            "mapped_to_Shumari_water_permeability.maf_0.05.assoc.txt",
            "mapped_to_Shumari_red_seedcoat.maf_0.05.assoc.txt",
            "mapped_to_Shumari_mottled_black_seedcoat.maf_0.05.assoc.txt",
        ]
        for name in other_files:
            self.assertFalse((self.data_dir / name).exists())
        # Must not raise despite those files being absent.
        pipeline.ensure_validated(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
        )


class LoadValidatedFrameTests(PipelineTestCase):
    def test_loader_not_called_when_validation_fails(self) -> None:
        self._write_invalid_dataset()
        with mock.patch(
            "adzuki_gwas_analysis.analysis.pipeline.load_analysis_frame"
        ) as mocked_load:
            with self.assertRaises(DatasetValidationFailedError):
                pipeline.load_validated_frame(
                    manifest_path=self.manifest_path,
                    data_dir=self.data_dir,
                    dataset_id="miyagi_water_permeability",
                )
            mocked_load.assert_not_called()


class RunManhattanTests(PipelineTestCase):
    def test_produces_png_for_valid_dataset(self) -> None:
        self._write_valid_dataset()
        output_path = self.output_dir / "manhattan.png"
        result = pipeline.run_manhattan(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_path=output_path,
        )
        self.assertTrue(output_path.is_file())
        self.assertEqual(result.output_path, output_path)

    def test_produces_no_output_when_validation_fails(self) -> None:
        self._write_invalid_dataset()
        output_path = self.output_dir / "manhattan.png"
        with self.assertRaises(DatasetValidationFailedError):
            pipeline.run_manhattan(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_path=output_path,
            )
        self.assertFalse(output_path.exists())

    def test_rejects_out_of_range_threshold(self) -> None:
        self._write_valid_dataset()
        with self.assertRaises(ValueError):
            pipeline.run_manhattan(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_path=self.output_dir / "manhattan.png",
                threshold=1.5,
            )

    def test_creates_output_parent_directory(self) -> None:
        self._write_valid_dataset()
        nested = self.output_dir / "a" / "b" / "c" / "manhattan.png"
        pipeline.run_manhattan(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_path=nested,
        )
        self.assertTrue(nested.is_file())


class RunQqTests(PipelineTestCase):
    def test_produces_png_for_valid_dataset(self) -> None:
        self._write_valid_dataset()
        output_path = self.output_dir / "qq.png"
        pipeline.run_qq(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_path=output_path,
        )
        self.assertTrue(output_path.is_file())

    def test_produces_no_output_when_validation_fails(self) -> None:
        self._write_invalid_dataset()
        output_path = self.output_dir / "qq.png"
        with self.assertRaises(DatasetValidationFailedError):
            pipeline.run_qq(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_path=output_path,
            )
        self.assertFalse(output_path.exists())


class RunSingleRegionalTests(PipelineTestCase):
    def test_top_variant_matches_expected_row(self) -> None:
        self._write_valid_dataset()
        outcome = pipeline.run_single_regional(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            chrom="Chr01",
            start=1_000_000,
            end=3_000_000,
            output_path=self.output_dir / "regional.png",
        )
        self.assertEqual(outcome.top_variant.pos, 1_000_000)
        self.assertEqual(outcome.top_variant.pval, 1e-4)
        self.assertTrue((self.output_dir / "regional.png").is_file())

    def test_empty_region_raises_and_produces_no_output(self) -> None:
        self._write_valid_dataset()
        output_path = self.output_dir / "regional.png"
        with self.assertRaises(EmptyRegionError):
            pipeline.run_single_regional(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                chrom="Chr01",
                start=9_000_000,
                end=9_100_000,
                output_path=output_path,
            )
        self.assertFalse(output_path.exists())


class RunRegionsTests(PipelineTestCase):
    def test_produces_one_plot_per_region(self) -> None:
        self._write_valid_dataset()
        config_path = write_region_config(self.tmp_path / "regions.toml", _ONE_REGION)
        outcomes = pipeline.run_regions(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            regions_config_path=config_path,
            output_dir=self.output_dir,
        )
        self.assertEqual(len(outcomes), 1)
        self.assertTrue((self.output_dir / "region1.png").is_file())

    def test_empty_region_raises_diagnosable_error(self) -> None:
        self._write_valid_dataset()
        config_path = write_region_config(self.tmp_path / "regions.toml", _EMPTY_REGION)
        with self.assertRaises(EmptyRegionError) as ctx:
            pipeline.run_regions(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                regions_config_path=config_path,
                output_dir=self.output_dir,
            )
        self.assertEqual(ctx.exception.region_id, "Chr01_empty")


class RunTopVariantsTests(PipelineTestCase):
    def test_writes_expected_tsv(self) -> None:
        self._write_valid_dataset()
        config_path = write_region_config(self.tmp_path / "regions.toml", _ONE_REGION)
        output_path = self.output_dir / "top_variants.tsv"
        table = pipeline.run_top_variants(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            regions_config_path=config_path,
            output_path=output_path,
        )
        self.assertTrue(output_path.is_file())
        self.assertEqual(list(table["region"]), ["Chr01_1_3Mb"])
        self.assertEqual(table["pos"].iloc[0], 1_000_000)

    def test_no_output_when_validation_fails(self) -> None:
        self._write_invalid_dataset()
        config_path = write_region_config(self.tmp_path / "regions.toml", _ONE_REGION)
        output_path = self.output_dir / "top_variants.tsv"
        with self.assertRaises(DatasetValidationFailedError):
            pipeline.run_top_variants(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                regions_config_path=config_path,
                output_path=output_path,
            )
        self.assertFalse(output_path.exists())


class RunAllTests(PipelineTestCase):
    def test_produces_all_expected_outputs_and_nothing_else(self) -> None:
        self._write_valid_dataset()
        config_path = write_region_config(self.tmp_path / "regions.toml", _ONE_REGION)
        outcome = pipeline.run_all(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            regions_config_path=config_path,
            output_dir=self.output_dir,
        )
        produced = {p.name for p in self.output_dir.iterdir()}
        expected = {
            "miyagi_water_permeability_manhattan.png",
            "miyagi_water_permeability_qq.png",
            "region1.png",
            "top_variants_by_region.tsv",
        }
        self.assertEqual(produced, expected)
        self.assertEqual(len(outcome.regions), 1)


if __name__ == "__main__":
    unittest.main()
