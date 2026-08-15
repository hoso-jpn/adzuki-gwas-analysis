"""Unit tests for adzuki_gwas_analysis.validate, using small deterministic fixtures.

Every test here builds its own :class:`~adzuki_gwas_analysis.manifest.DatasetEntry`
directly (rather than round-tripping through a manifest TOML file), pointing
at one of the small fixtures under ``tests/fixtures/``. Checksums and row
counts are computed from the fixture itself at test time via the same
loader functions the library uses, so nothing here is a hand-copied value
that could silently drift from the fixture's real content.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from adzuki_gwas_analysis.loader import compute_sha256, count_data_rows
from adzuki_gwas_analysis.manifest import DatasetEntry
from adzuki_gwas_analysis.validate import validate_dataset

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def _entry_for(
    filename: str, *, sha256: str | None = None, row_count: int | None = None
) -> DatasetEntry:
    path = FIXTURES_DIR / filename
    return DatasetEntry(
        dataset_id="test_dataset",
        reference="Miyagi",
        trait="water_permeability",
        member_filename=filename,
        member_sha256=sha256 if sha256 is not None else compute_sha256(path),
        row_count=row_count if row_count is not None else count_data_rows(path),
    )


class ValidateDatasetTests(unittest.TestCase):
    def test_valid_fixture_passes(self) -> None:
        entry = _entry_for("valid.assoc.txt")
        result = validate_dataset(entry, FIXTURES_DIR)
        self.assertTrue(result.success, msg=result.error)
        self.assertEqual(result.row_count, 3)
        self.assertIsNone(result.error)

    def test_missing_file_fails_with_diagnosable_error(self) -> None:
        entry = DatasetEntry(
            dataset_id="test_dataset",
            reference="Miyagi",
            trait="water_permeability",
            member_filename="does_not_exist.assoc.txt",
            member_sha256="0" * 64,
            row_count=3,
        )
        result = validate_dataset(entry, FIXTURES_DIR)
        self.assertFalse(result.success)
        assert result.error is not None
        self.assertIn("not found", result.error)
        self.assertIn(entry.dataset_id, result.error)

    def test_checksum_mismatch_fails_with_diagnosable_error(self) -> None:
        entry = _entry_for("valid.assoc.txt", sha256="0" * 64)
        result = validate_dataset(entry, FIXTURES_DIR)
        self.assertFalse(result.success)
        assert result.error is not None
        self.assertIn("checksum mismatch", result.error)

    def test_row_count_mismatch_fails_with_diagnosable_error(self) -> None:
        entry = _entry_for("valid.assoc.txt", row_count=999)
        result = validate_dataset(entry, FIXTURES_DIR)
        self.assertFalse(result.success)
        assert result.error is not None
        self.assertIn("row count mismatch", result.error)
        self.assertIn("999", result.error)

    def test_duplicate_variant_fails_end_to_end(self) -> None:
        entry = _entry_for("duplicate_variant.assoc.txt")
        result = validate_dataset(entry, FIXTURES_DIR)
        self.assertFalse(result.success)
        assert result.error is not None
        self.assertIn("duplicate variant", result.error.lower())

    def test_bad_row_fails_end_to_end(self) -> None:
        entry = _entry_for("bad_numeric_pos.assoc.txt")
        result = validate_dataset(entry, FIXTURES_DIR)
        self.assertFalse(result.success)
        assert result.error is not None
        self.assertIn("pos", result.error)


if __name__ == "__main__":
    unittest.main()
