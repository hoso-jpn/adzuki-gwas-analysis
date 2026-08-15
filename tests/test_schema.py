"""Unit tests for adzuki_gwas_analysis.schema, using small deterministic fixtures."""

from __future__ import annotations

import unittest
from pathlib import Path

from adzuki_gwas_analysis.errors import RowValidationError, SchemaError
from adzuki_gwas_analysis.loader import iter_data_rows, read_header
from adzuki_gwas_analysis.schema import REQUIRED_COLUMNS, validate_header, validate_row

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def _rows(filename: str) -> list[tuple[int, dict[str, str]]]:
    path = FIXTURES_DIR / filename
    header = read_header(path)
    return list(iter_data_rows(dataset_id="test", path=path, header=header))


class ValidateHeaderTests(unittest.TestCase):
    def test_valid_header_matches_required_columns_exactly(self) -> None:
        header = read_header(FIXTURES_DIR / "valid.assoc.txt")
        validate_header(dataset_id="test", path="valid.assoc.txt", header=header)
        self.assertEqual(tuple(header), REQUIRED_COLUMNS)

    def test_missing_required_column_raises_schema_error(self) -> None:
        header = read_header(FIXTURES_DIR / "missing_column.assoc.txt")
        with self.assertRaises(SchemaError) as ctx:
            validate_header(dataset_id="test", path="missing_column.assoc.txt", header=header)
        self.assertIn("af", str(ctx.exception))

    def test_duplicate_header_column_raises_schema_error(self) -> None:
        header = [*REQUIRED_COLUMNS, "pval"]
        with self.assertRaises(SchemaError) as ctx:
            validate_header(dataset_id="test", path="synthetic", header=header)
        self.assertIn("duplicate", str(ctx.exception).lower())


class ValidateRowTests(unittest.TestCase):
    def test_valid_rows_pass_and_rs_dot_is_permitted(self) -> None:
        rows = _rows("valid.assoc.txt")
        self.assertEqual(len(rows), 3)
        for row_number, row in rows:
            self.assertEqual(row["rs"], ".")
            key = validate_row(
                dataset_id="test", path="valid.assoc.txt", row_number=row_number, row=row
            )
            self.assertEqual(len(key), 4)

    def test_bad_numeric_pos_raises_row_validation_error(self) -> None:
        row_number, row = _rows("bad_numeric_pos.assoc.txt")[0]
        with self.assertRaises(RowValidationError) as ctx:
            validate_row(dataset_id="test", path="x", row_number=row_number, row=row)
        self.assertEqual(ctx.exception.column, "pos")

    def test_pvalue_zero_raises(self) -> None:
        row_number, row = _rows("pvalue_zero.assoc.txt")[0]
        with self.assertRaises(RowValidationError) as ctx:
            validate_row(dataset_id="test", path="x", row_number=row_number, row=row)
        self.assertEqual(ctx.exception.column, "pval")

    def test_pvalue_greater_than_one_raises(self) -> None:
        row_number, row = _rows("pvalue_gt1.assoc.txt")[0]
        with self.assertRaises(RowValidationError) as ctx:
            validate_row(dataset_id="test", path="x", row_number=row_number, row=row)
        self.assertEqual(ctx.exception.column, "pval")

    def test_af_negative_raises(self) -> None:
        row_number, row = _rows("af_negative.assoc.txt")[0]
        with self.assertRaises(RowValidationError) as ctx:
            validate_row(dataset_id="test", path="x", row_number=row_number, row=row)
        self.assertEqual(ctx.exception.column, "af")

    def test_af_greater_than_one_raises(self) -> None:
        row_number, row = _rows("af_gt1.assoc.txt")[0]
        with self.assertRaises(RowValidationError) as ctx:
            validate_row(dataset_id="test", path="x", row_number=row_number, row=row)
        self.assertEqual(ctx.exception.column, "af")

    def test_pos_nonpositive_raises(self) -> None:
        row_number, row = _rows("pos_nonpositive.assoc.txt")[0]
        with self.assertRaises(RowValidationError) as ctx:
            validate_row(dataset_id="test", path="x", row_number=row_number, row=row)
        self.assertEqual(ctx.exception.column, "pos")

    def test_n_miss_negative_raises(self) -> None:
        row_number, row = _rows("n_miss_negative.assoc.txt")[0]
        with self.assertRaises(RowValidationError) as ctx:
            validate_row(dataset_id="test", path="x", row_number=row_number, row=row)
        self.assertEqual(ctx.exception.column, "n_miss")

    def test_n_miss_noninteger_raises(self) -> None:
        row_number, row = _rows("n_miss_noninteger.assoc.txt")[0]
        with self.assertRaises(RowValidationError) as ctx:
            validate_row(dataset_id="test", path="x", row_number=row_number, row=row)
        self.assertEqual(ctx.exception.column, "n_miss")

    def test_duplicate_variant_fixture_rows_share_the_same_key(self) -> None:
        # schema.validate_row only computes the variant key; detecting that
        # two rows share a key is validate.validate_dataset's job (see
        # test_validate.py). This test just proves the fixture's two rows
        # really do produce an identical key, so that end-to-end test is
        # exercising what it claims to.
        rows = _rows("duplicate_variant.assoc.txt")
        keys = [validate_row(dataset_id="test", path="x", row_number=n, row=r) for n, r in rows]
        self.assertEqual(keys[0], keys[1])


if __name__ == "__main__":
    unittest.main()
