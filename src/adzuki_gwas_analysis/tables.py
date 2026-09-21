"""Strict TSV transport preserving identifiers and explicit missing-value strings."""

from __future__ import annotations

import csv
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any


def read_tsv(path: Path, required: Iterable[str] = ()) -> Iterator[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if not header or len(set(header)) != len(header) or not set(required) <= set(header):
            raise ValueError("TSV header is empty, duplicated or missing required columns")
        for number, fields in enumerate(reader, 2):
            if len(fields) != len(header) or any("\n" in x or "\r" in x for x in fields):
                raise ValueError(f"invalid TSV record at line {number}")
            yield dict(zip(header, fields, strict=True))


def write_tsv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
