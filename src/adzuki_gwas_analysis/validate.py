"""Orchestrates end-to-end validation of one manifest dataset entry against its real file.

:func:`validate_dataset` is deliberately the only entry point most callers
need: it checks file existence, checksum, header/schema, every data row, and
duplicate-variant detection, in that order, stopping at the first
contract violation. It processes exactly one file at a time and never holds
another dataset's rows in memory.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from adzuki_gwas_analysis.errors import (
    ChecksumMismatchError,
    DatasetFileMissingError,
    DuplicateVariantError,
    GwasContractError,
    RowCountMismatchError,
)
from adzuki_gwas_analysis.loader import compute_sha256, iter_data_rows, read_header
from adzuki_gwas_analysis.manifest import DatasetEntry
from adzuki_gwas_analysis.schema import VariantKey, validate_header, validate_row


@dataclass(frozen=True, slots=True)
class ValidationResult:
    """The outcome of validating one dataset entry against its real file."""

    dataset_id: str
    reference: str
    trait: str
    filename: str
    row_count: int
    sha256: str
    success: bool
    error: str | None = None


def validate_dataset(entry: DatasetEntry, data_dir: Path) -> ValidationResult:
    """Validate ``entry`` against ``data_dir / entry.member_filename``.

    Returns a :class:`ValidationResult` describing success or the first
    contract violation found; does not raise for expected contract
    violations (missing file, checksum/row-count mismatch, schema error, a
    bad row, a duplicate variant) so a caller can validate multiple
    datasets and report every failure rather than stopping at the first one.
    """
    path = data_dir / entry.member_filename
    try:
        if not path.is_file():
            raise DatasetFileMissingError(dataset_id=entry.dataset_id, path=str(path))

        actual_sha256 = compute_sha256(path)
        if actual_sha256 != entry.member_sha256:
            raise ChecksumMismatchError(
                dataset_id=entry.dataset_id,
                path=str(path),
                expected_sha256=entry.member_sha256,
                actual_sha256=actual_sha256,
            )

        header = read_header(path)
        validate_header(dataset_id=entry.dataset_id, path=str(path), header=header)

        seen_keys: dict[VariantKey, int] = {}
        row_count = 0
        for row_number, row in iter_data_rows(
            dataset_id=entry.dataset_id, path=path, header=header
        ):
            key = validate_row(
                dataset_id=entry.dataset_id, path=str(path), row_number=row_number, row=row
            )
            first_seen_at = seen_keys.get(key)
            if first_seen_at is not None:
                raise DuplicateVariantError(
                    dataset_id=entry.dataset_id,
                    path=str(path),
                    key=key,
                    first_row_number=first_seen_at,
                    duplicate_row_number=row_number,
                )
            seen_keys[key] = row_number
            row_count = row_number

        if row_count != entry.row_count:
            raise RowCountMismatchError(
                dataset_id=entry.dataset_id,
                path=str(path),
                expected_row_count=entry.row_count,
                actual_row_count=row_count,
            )
    except GwasContractError as exc:
        return ValidationResult(
            dataset_id=entry.dataset_id,
            reference=entry.reference,
            trait=entry.trait,
            filename=entry.member_filename,
            row_count=0,
            sha256="",
            success=False,
            error=str(exc),
        )

    return ValidationResult(
        dataset_id=entry.dataset_id,
        reference=entry.reference,
        trait=entry.trait,
        filename=entry.member_filename,
        row_count=row_count,
        sha256=actual_sha256,
        success=True,
        error=None,
    )


def validate_all(entries: Iterable[DatasetEntry], data_dir: Path) -> list[ValidationResult]:
    """Validate every entry in ``entries`` against ``data_dir``, one file at a time."""
    return [validate_dataset(entry, data_dir) for entry in entries]
