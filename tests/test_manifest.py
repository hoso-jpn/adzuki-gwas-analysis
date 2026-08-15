"""Unit tests for adzuki_gwas_analysis.manifest."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.errors import ManifestError
from adzuki_gwas_analysis.manifest import load_manifest

REPO_ROOT = Path(__file__).resolve().parents[1]
REAL_MANIFEST_PATH = REPO_ROOT / "manifest.toml"

_VALID_MANIFEST_TOML = """
schema_version = 1

[dryad]
doi = "10.5061/dryad.8w9ghx3xv"
dataset_id = 149675
version_id = 356599
version_number = 6
publication_doi = "10.1126/science.ads2871"

[archive]
filename = "adzuki_GWAS_data.zip"
size_bytes = 257167472
sha256 = "92f0bc6a235c3ee863e96736a831f8652569269a8f61c8699cc9f562743b0ed5"

[pvalue_columns]
p_wald = "Wald test p-value"
pval = "Likelihood ratio test (LRT) p-value"
p_score = "Score test p-value"
primary = "pval"

[[datasets]]
dataset_id = "miyagi_water_permeability"
reference = "Miyagi"
trait = "water_permeability"
member_filename = "mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"
member_sha256 = "1b3b3510fd37a9c8f502b7dbb0450433200afd25750e7ff4e22b8582c05f41d6"
row_count = 1741385
"""


def _write(tmp_dir: str, content: str) -> Path:
    path = Path(tmp_dir) / "manifest.toml"
    path.write_text(content, encoding="utf-8")
    return path


class LoadManifestTests(unittest.TestCase):
    def test_valid_manifest_loads_successfully(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = _write(tmp_dir, _VALID_MANIFEST_TOML)
            manifest = load_manifest(path)
        self.assertEqual(manifest.schema_version, 1)
        self.assertEqual(manifest.archive.filename, "adzuki_GWAS_data.zip")
        self.assertEqual(len(manifest.datasets), 1)
        self.assertEqual(manifest.pvalue_columns.primary, "pval")

    def test_missing_manifest_file_raises_manifest_error(self) -> None:
        with self.assertRaises(ManifestError):
            load_manifest("/nonexistent/path/manifest.toml")

    def test_unsupported_schema_version_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace("schema_version = 1", "schema_version = 99")
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = _write(tmp_dir, broken)
            with self.assertRaises(ManifestError) as ctx:
                load_manifest(path)
        self.assertIn("schema_version", str(ctx.exception))

    def test_duplicate_dataset_id_raises_manifest_error(self) -> None:
        # Append the single [[datasets]] block a second time, verbatim, so
        # the same dataset_id ("miyagi_water_permeability") appears twice.
        datasets_block = _VALID_MANIFEST_TOML[_VALID_MANIFEST_TOML.index("[[datasets]]") :]
        duplicated = _VALID_MANIFEST_TOML + "\n" + datasets_block
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = _write(tmp_dir, duplicated)
            with self.assertRaises(ManifestError) as ctx:
                load_manifest(path)
        self.assertIn("duplicate", str(ctx.exception).lower())

    def test_missing_required_field_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace('doi = "10.5061/dryad.8w9ghx3xv"\n', "")
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = _write(tmp_dir, broken)
            with self.assertRaises(ManifestError) as ctx:
                load_manifest(path)
        self.assertIn("doi", str(ctx.exception))

    def test_invalid_reference_value_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace('reference = "Miyagi"', 'reference = "Nagano"')
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = _write(tmp_dir, broken)
            with self.assertRaises(ManifestError):
                load_manifest(path)


class RealManifestTests(unittest.TestCase):
    """These tests load only the small, committed manifest.toml -- never the
    multi-hundred-megabyte raw data files -- so they run safely in CI."""

    def test_real_manifest_has_exactly_six_datasets(self) -> None:
        manifest = load_manifest(REAL_MANIFEST_PATH)
        self.assertEqual(len(manifest.datasets), 6)

    def test_real_manifest_has_all_six_reference_trait_combinations(self) -> None:
        manifest = load_manifest(REAL_MANIFEST_PATH)
        observed = {(entry.reference, entry.trait) for entry in manifest.datasets}
        expected = {
            ("Miyagi", "water_permeability"),
            ("Miyagi", "red_seedcoat"),
            ("Miyagi", "mottled_black_seedcoat"),
            ("Shumari", "water_permeability"),
            ("Shumari", "red_seedcoat"),
            ("Shumari", "mottled_black_seedcoat"),
        }
        self.assertEqual(observed, expected)

    def test_real_manifest_dataset_ids_are_unique(self) -> None:
        manifest = load_manifest(REAL_MANIFEST_PATH)
        dataset_ids = [entry.dataset_id for entry in manifest.datasets]
        self.assertEqual(len(dataset_ids), len(set(dataset_ids)))

    def test_real_manifest_primary_pvalue_column_is_pval(self) -> None:
        manifest = load_manifest(REAL_MANIFEST_PATH)
        self.assertEqual(manifest.pvalue_columns.primary, "pval")


if __name__ == "__main__":
    unittest.main()
