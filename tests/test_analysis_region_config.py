"""Unit tests for adzuki_gwas_analysis.analysis.region_config."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.analysis.region_config import load_region_config
from adzuki_gwas_analysis.errors import RegionConfigError

REPO_ROOT = Path(__file__).resolve().parents[1]
REAL_REGIONS_CONFIG_PATH = REPO_ROOT / "config" / "water_permeability_regions.toml"

_VALID = """
config_schema_version = 1
dataset_id = "miyagi_water_permeability"

[[regions]]
region_id = "region_a"
chrom = "Chr01"
start = 100
end = 200
title = "Region A"
output_filename = "region_a.png"

[[regions]]
region_id = "region_b"
chrom = "Chr02"
start = 1
end = 2
title = "Region B"
output_filename = "region_b.png"
"""


class LoadRegionConfigTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.config_path = Path(self._tmpdir.name) / "regions.toml"

    def _write(self, content: str) -> Path:
        self.config_path.write_text(content, encoding="utf-8")
        return self.config_path

    def test_real_repository_config_is_valid(self) -> None:
        config = load_region_config(
            REAL_REGIONS_CONFIG_PATH, expected_dataset_id="miyagi_water_permeability"
        )
        self.assertEqual(len(config.regions), 5)
        region_ids = [r.region_id for r in config.regions]
        self.assertEqual(len(region_ids), len(set(region_ids)))

    def test_valid_config_loads(self) -> None:
        path = self._write(_VALID)
        config = load_region_config(path, expected_dataset_id="miyagi_water_permeability")
        self.assertEqual(config.config_schema_version, 1)
        self.assertEqual(len(config.regions), 2)
        self.assertEqual(config.regions[0].region_id, "region_a")

    def test_unsupported_schema_version(self) -> None:
        path = self._write(_VALID.replace("config_schema_version = 1", "config_schema_version = 2"))
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_dataset_id_mismatch(self) -> None:
        path = self._write(_VALID)
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="shumari_water_permeability")

    def test_duplicate_region_id(self) -> None:
        path = self._write(_VALID.replace('region_id = "region_b"', 'region_id = "region_a"'))
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_unsafe_region_id_rejected(self) -> None:
        path = self._write(_VALID.replace('region_id = "region_a"', 'region_id = "../etc"'))
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_duplicate_output_filename(self) -> None:
        path = self._write(
            _VALID.replace('output_filename = "region_b.png"', 'output_filename = "region_a.png"')
        )
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_unsafe_output_filename_path_traversal(self) -> None:
        path = self._write(
            _VALID.replace('output_filename = "region_a.png"', 'output_filename = "../a.png"')
        )
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_unsafe_output_filename_path_separator(self) -> None:
        path = self._write(
            _VALID.replace('output_filename = "region_a.png"', 'output_filename = "sub/a.png"')
        )
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_empty_chrom_rejected(self) -> None:
        path = self._write(_VALID.replace('chrom = "Chr01"', 'chrom = ""'))
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_start_not_positive_rejected(self) -> None:
        path = self._write(_VALID.replace("start = 100", "start = 0"))
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_end_before_start_rejected(self) -> None:
        path = self._write(_VALID.replace("start = 100\nend = 200", "start = 200\nend = 100"))
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_zero_regions_rejected(self) -> None:
        path = self._write(
            'config_schema_version = 1\ndataset_id = "miyagi_water_permeability"\nregions = []\n'
        )
        with self.assertRaises(RegionConfigError):
            load_region_config(path, expected_dataset_id="miyagi_water_permeability")

    def test_missing_file(self) -> None:
        with self.assertRaises(RegionConfigError):
            load_region_config(
                Path(self._tmpdir.name) / "nope.toml",
                expected_dataset_id="miyagi_water_permeability",
            )


if __name__ == "__main__":
    unittest.main()
