"""Unit tests for adzuki_gwas_analysis.manifest.

``_VALID_MANIFEST_TOML`` is a fully synthetic, schema-v1-compliant 6-dataset
manifest (fake but correctly-shaped SHA-256 digests -- ``"a" * 64`` etc --
so it is never mistaken for real Dryad provenance). Each negative test
derives a broken variant from it, so every test isolates exactly one
contract violation.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.errors import ManifestError
from adzuki_gwas_analysis.manifest import load_manifest

REPO_ROOT = Path(__file__).resolve().parents[1]
REAL_MANIFEST_PATH = REPO_ROOT / "manifest.toml"


def _fake_sha256(marker: str) -> str:
    """A syntactically-valid (64 lowercase hex chars) but obviously-fake digest.

    Computed, not hand-typed, so its length is never in question.
    """
    return (marker * 64)[:64]


_ARCHIVE_SHA256 = _fake_sha256("0")

_PREAMBLE = f"""
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
sha256 = "{_ARCHIVE_SHA256}"

[pvalue_columns]
p_wald = "Wald test p-value"
pval = "Likelihood ratio test (LRT) p-value"
p_score = "Score test p-value"
primary = "pval"
"""


# One [[datasets]] block per canonical (reference, trait) pair, in the same
# order as adzuki_gwas_analysis.manifest.SCHEMA_V1_DATASETS. member_sha256
# values are computed via _fake_sha256(), not hand-typed, so their length is
# never in question.
_DATASET_SPECS: tuple[tuple[str, str, str, str, int], ...] = (
    (
        "miyagi_water_permeability",
        "Miyagi",
        "water_permeability",
        "mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt",
        1741385,
    ),
    (
        "miyagi_red_seedcoat",
        "Miyagi",
        "red_seedcoat",
        "mapped_to_Miyagi_red_seedcoat.maf_0.05.assoc.txt",
        1255203,
    ),
    (
        "miyagi_mottled_black_seedcoat",
        "Miyagi",
        "mottled_black_seedcoat",
        "mapped_to_Miyagi_mottled_black_seedcoat.maf_0.05.assoc.txt",
        1255203,
    ),
    (
        "shumari_water_permeability",
        "Shumari",
        "water_permeability",
        "mapped_to_Shumari_water_permeability.maf_0.05.assoc.txt",
        1471837,
    ),
    (
        "shumari_red_seedcoat",
        "Shumari",
        "red_seedcoat",
        "mapped_to_Shumari_red_seedcoat.maf_0.05.assoc.txt",
        1232183,
    ),
    (
        "shumari_mottled_black_seedcoat",
        "Shumari",
        "mottled_black_seedcoat",
        "mapped_to_Shumari_mottled_black_seedcoat.maf_0.05.assoc.txt",
        1232183,
    ),
)


def _dataset_block(
    dataset_id: str,
    reference: str,
    trait: str,
    member_filename: str,
    row_count: int,
    *,
    sha_marker: str,
) -> str:
    return f"""
[[datasets]]
dataset_id = "{dataset_id}"
reference = "{reference}"
trait = "{trait}"
member_filename = "{member_filename}"
member_sha256 = "{_fake_sha256(sha_marker)}"
row_count = {row_count}
"""


# sha_marker is a single lowercase hex digit ("a".."f") per dataset, distinct
# per entry and from the archive's own "0" marker -- dataset_id itself isn't
# usable here since e.g. "m" (from "miyagi_...") isn't a valid hex character,
# and a marker needs to include a letter (not just a digit) for the
# uppercase-rejection test below to mean anything (digits have no case).
_DATASET_BLOCKS: dict[str, str] = {
    spec[0]: _dataset_block(*spec, sha_marker=format(index + 10, "x"))
    for index, spec in enumerate(_DATASET_SPECS)
}


def _manifest_toml(*dataset_ids: str) -> str:
    return _PREAMBLE + "".join(_DATASET_BLOCKS[dataset_id] for dataset_id in dataset_ids)


_ALL_SIX_IDS = tuple(_DATASET_BLOCKS.keys())
_VALID_MANIFEST_TOML = _manifest_toml(*_ALL_SIX_IDS)


def _write(tmp_dir: str, content: str) -> Path:
    path = Path(tmp_dir) / "manifest.toml"
    path.write_text(content, encoding="utf-8")
    return path


def _load(content: str) -> object:
    with tempfile.TemporaryDirectory() as tmp_dir:
        return load_manifest(_write(tmp_dir, content))


class LoadManifestTests(unittest.TestCase):
    def test_valid_six_dataset_manifest_loads_successfully(self) -> None:
        manifest = _load(_VALID_MANIFEST_TOML)
        self.assertEqual(manifest.schema_version, 1)
        self.assertEqual(manifest.archive.filename, "adzuki_GWAS_data.zip")
        self.assertEqual(len(manifest.datasets), 6)
        self.assertEqual(manifest.pvalue_columns.primary, "pval")

    def test_missing_manifest_file_raises_manifest_error(self) -> None:
        with self.assertRaises(ManifestError):
            load_manifest("/nonexistent/path/manifest.toml")

    def test_invalid_toml_raises_manifest_error_not_a_traceback(self) -> None:
        broken = _VALID_MANIFEST_TOML + "\nthis is not valid toml [[["
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("toml", str(ctx.exception).lower())

    def test_unsupported_schema_version_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace("schema_version = 1", "schema_version = 99")
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("schema_version", str(ctx.exception))

    def test_missing_required_field_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace('doi = "10.5061/dryad.8w9ghx3xv"\n', "")
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("doi", str(ctx.exception))

    def test_invalid_reference_value_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace('reference = "Miyagi"', 'reference = "Nagano"', 1)
        with self.assertRaises(ManifestError):
            _load(broken)

    # -- P1-1: schema v1's exactly-6-datasets contract --

    def test_one_dataset_manifest_raises_manifest_error(self) -> None:
        with self.assertRaises(ManifestError) as ctx:
            _load(_manifest_toml("miyagi_water_permeability"))
        self.assertIn("6", str(ctx.exception))

    def test_five_dataset_manifest_raises_manifest_error(self) -> None:
        with self.assertRaises(ManifestError) as ctx:
            _load(_manifest_toml(*_ALL_SIX_IDS[:5]))
        self.assertIn("6", str(ctx.exception))

    def test_seven_dataset_manifest_raises_manifest_error(self) -> None:
        # 6 valid entries plus a verbatim repeat of the first.
        seven = _VALID_MANIFEST_TOML + _DATASET_BLOCKS["miyagi_water_permeability"]
        with self.assertRaises(ManifestError):
            _load(seven)

    def test_duplicate_reference_trait_combination_raises_manifest_error(self) -> None:
        # A 7th entry with a brand-new dataset_id/member_filename (so it
        # does NOT trip the dataset_id/member_filename duplicate checks)
        # but the same (reference, trait) as an existing entry -- isolates
        # the (reference, trait) uniqueness check specifically.
        extra_block = """
[[datasets]]
dataset_id = "miyagi_water_permeability_dup"
reference = "Miyagi"
trait = "water_permeability"
member_filename = "mapped_to_Miyagi_water_permeability_dup.maf_0.05.assoc.txt"
member_sha256 = "{sha}"
row_count = 1741385
""".format(sha=_fake_sha256("f"))
        with self.assertRaises(ManifestError) as ctx:
            _load(_VALID_MANIFEST_TOML + extra_block)
        self.assertIn("duplicate", str(ctx.exception).lower())
        self.assertIn("reference", str(ctx.exception).lower())

    def test_unknown_trait_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace(
            'trait = "water_permeability"', 'trait = "unknown_trait"', 1
        )
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("trait", str(ctx.exception).lower())

    def test_dataset_id_mismatched_with_reference_trait_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace(
            'dataset_id = "miyagi_water_permeability"',
            'dataset_id = "miyagi_water_permeability_wrong"',
            1,
        )
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("dataset_id", str(ctx.exception))

    def test_member_filename_mismatched_with_canonical_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace(
            'member_filename = "mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"',
            'member_filename = "wrong_name.maf_0.05.assoc.txt"',
            1,
        )
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("member_filename", str(ctx.exception))

    def test_primary_pvalue_column_must_be_pval(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace('primary = "pval"', 'primary = "p_wald"')
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("pval", str(ctx.exception))

    # -- P2-2: filename safety --

    def test_absolute_path_member_filename_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace(
            'member_filename = "mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"',
            'member_filename = "/etc/passwd.assoc.txt"',
            1,
        )
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("absolute", str(ctx.exception).lower())

    def test_path_traversal_member_filename_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace(
            'member_filename = "mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"',
            'member_filename = "../outside.assoc.txt"',
            1,
        )
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("..", str(ctx.exception))

    def test_wrong_suffix_member_filename_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace(
            'member_filename = "mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"',
            'member_filename = "mapped_to_Miyagi_water_permeability.maf_0.05.txt"',
            1,
        )
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("assoc.txt", str(ctx.exception))

    # -- P2-2: numeric / checksum format --

    def test_non_positive_archive_size_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace("size_bytes = 257167472", "size_bytes = 0")
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("size_bytes", str(ctx.exception))

    def test_non_positive_row_count_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace("row_count = 1741385", "row_count = 0", 1)
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("row_count", str(ctx.exception))

    def test_malformed_archive_sha256_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace(
            f'sha256 = "{_ARCHIVE_SHA256}"',
            'sha256 = "not-a-checksum"',
        )
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("sha256", str(ctx.exception).lower())

    def test_uppercase_member_sha256_raises_manifest_error(self) -> None:
        lowercase_digest = _fake_sha256("a")
        broken = _VALID_MANIFEST_TOML.replace(
            f'member_sha256 = "{lowercase_digest}"',
            f'member_sha256 = "{lowercase_digest.upper()}"',
            1,
        )
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("sha256", str(ctx.exception).lower())

    def test_non_positive_dryad_dataset_id_raises_manifest_error(self) -> None:
        broken = _VALID_MANIFEST_TOML.replace("dataset_id = 149675", "dataset_id = -1")
        with self.assertRaises(ManifestError) as ctx:
            _load(broken)
        self.assertIn("dataset_id", str(ctx.exception))


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
