"""Streaming I/O for GWAS summary-statistics files.

Every function here processes one file at a time and never loads a whole
file into memory: :func:`compute_sha256` streams in fixed-size chunks, and
:func:`iter_data_rows` yields one row at a time from an open file handle.
This matters because real files in this dataset are up to ~1.7M data rows
(~230 MB) each, and 6 of them must never be held in memory simultaneously.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Iterator
from pathlib import Path

from adzuki_gwas_analysis.errors import RowValidationError

DEFAULT_CHUNK_SIZE = 1024 * 1024


def compute_sha256(path: Path, *, chunk_size: int = DEFAULT_CHUNK_SIZE) -> str:
    """Compute the SHA-256 digest of ``path``, reading it in fixed-size chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def read_header(path: Path) -> list[str]:
    """Read and tab-split only the first line of ``path``."""
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        return next(reader, [])


def iter_data_rows(
    *, dataset_id: str, path: Path, header: list[str]
) -> Iterator[tuple[int, dict[str, str]]]:
    """Yield ``(row_number, row)`` for every data row in ``path`` after its header.

    ``row_number`` is 1-indexed and counts data rows only. Raises
    :class:`~adzuki_gwas_analysis.errors.RowValidationError` if a row has a
    different number of fields than ``header`` (e.g. a truncated line from
    an interrupted download), rather than silently zipping a short row.
    """
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # header is read/validated separately; skip it here
        for row_number, fields in enumerate(reader, start=1):
            if len(fields) != len(header):
                raise RowValidationError(
                    dataset_id=dataset_id,
                    path=str(path),
                    row_number=row_number,
                    column="<row>",
                    value="\t".join(fields),
                    reason=(
                        f"row has {len(fields)} field(s), expected {len(header)} "
                        f"matching the header"
                    ),
                )
            yield row_number, dict(zip(header, fields, strict=True))


def count_data_rows(path: Path) -> int:
    """Count data rows in ``path`` (total lines minus the header), streaming."""
    with path.open("rb") as fh:
        total_lines = sum(1 for _ in fh)
    return max(0, total_lines - 1)
