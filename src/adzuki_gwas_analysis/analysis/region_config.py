"""Machine-readable config for post-hoc GWAS visualization regions.

This is deliberately a separate file and schema from ``manifest.toml``:
``manifest.toml`` is the Dryad input contract (schema v1, #1/#2) -- what a
correct *input file* looks like. A visualization region is a completely
different concern -- an arbitrary window a human picked by eye after looking
at a Manhattan plot -- and mixing the two would make it look like regions are
part of the data's own contract. They are not: the regions here are **not**
independently defined QTL intervals or LD blocks, and nothing in this module
computes or implies genome-wide significance, linkage, or a claim of a novel
association.
"""

from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from adzuki_gwas_analysis.errors import RegionConfigError

SUPPORTED_CONFIG_SCHEMA_VERSION = 1

_REGION_ID_RE = re.compile(r"^[A-Za-z0-9_]+$")


@dataclass(frozen=True, slots=True)
class RegionSpec:
    """One post-hoc visualization window."""

    region_id: str
    chrom: str
    start: int
    end: int
    title: str
    output_filename: str


@dataclass(frozen=True, slots=True)
class RegionConfig:
    """A fully validated set of visualization regions for one dataset."""

    config_schema_version: int
    dataset_id: str
    regions: tuple[RegionSpec, ...]


def _require(table: dict[str, object], key: str, *, region_id: str | None = None) -> object:
    if key not in table:
        raise RegionConfigError(f"missing required field {key!r}", region_id=region_id)
    return table[key]


def _require_str(table: dict[str, object], key: str, *, region_id: str | None = None) -> str:
    value = _require(table, key, region_id=region_id)
    if not isinstance(value, str) or not value:
        raise RegionConfigError(f"field {key!r} must be a non-empty string", region_id=region_id)
    return value


def _require_int(table: dict[str, object], key: str, *, region_id: str | None = None) -> int:
    value = _require(table, key, region_id=region_id)
    if not isinstance(value, int) or isinstance(value, bool):
        raise RegionConfigError(f"field {key!r} must be an integer", region_id=region_id)
    return value


def _require_safe_output_filename(table: dict[str, object], region_id: str) -> str:
    value = _require_str(table, "output_filename", region_id=region_id)
    if PurePosixPath(value).is_absolute() or value.startswith(("/", "\\")):
        raise RegionConfigError(
            f"field 'output_filename' must not be an absolute path: {value!r}",
            region_id=region_id,
        )
    if "/" in value or "\\" in value:
        raise RegionConfigError(
            f"field 'output_filename' must be a plain filename with no path separators: {value!r}",
            region_id=region_id,
        )
    if ".." in value:
        raise RegionConfigError(
            f"field 'output_filename' must not contain '..': {value!r}", region_id=region_id
        )
    if not value.endswith(".png"):
        raise RegionConfigError(
            f"field 'output_filename' must end with '.png': {value!r}", region_id=region_id
        )
    return value


def load_region_config(path: str | Path, *, expected_dataset_id: str) -> RegionConfig:
    """Load and validate the region config at ``path``.

    ``expected_dataset_id`` must match the config's own ``dataset_id`` field,
    so a config file for one dataset can never be silently applied to
    another. Raises :class:`~adzuki_gwas_analysis.errors.RegionConfigError`
    for: a missing file; invalid TOML; an unsupported schema version; a
    ``dataset_id`` mismatch; zero regions; a duplicate or unsafe
    ``region_id``; a non-positive ``start``; an ``end`` before ``start``; an
    empty ``chrom``; or a duplicate/unsafe ``output_filename``.
    """
    config_path = Path(path)
    if not config_path.is_file():
        raise RegionConfigError(f"region config file not found: {config_path}")

    try:
        with config_path.open("rb") as fh:
            raw = tomllib.load(fh)
    except tomllib.TOMLDecodeError as exc:
        raise RegionConfigError(f"region config is not valid TOML: {exc}") from exc
    except OSError as exc:
        raise RegionConfigError(f"could not read region config file: {exc}") from exc

    schema_version = _require_int(raw, "config_schema_version")
    if schema_version != SUPPORTED_CONFIG_SCHEMA_VERSION:
        raise RegionConfigError(
            f"unsupported config_schema_version {schema_version}; "
            f"this package supports version {SUPPORTED_CONFIG_SCHEMA_VERSION}"
        )

    dataset_id = _require_str(raw, "dataset_id")
    if dataset_id != expected_dataset_id:
        raise RegionConfigError(
            f"region config dataset_id {dataset_id!r} does not match "
            f"requested dataset_id {expected_dataset_id!r}"
        )

    regions_raw = _require(raw, "regions")
    if not isinstance(regions_raw, list) or not regions_raw:
        raise RegionConfigError("field 'regions' must be a non-empty array of tables")

    seen_region_ids: set[str] = set()
    seen_output_filenames: set[str] = set()
    regions: list[RegionSpec] = []
    for entry_raw in regions_raw:
        if not isinstance(entry_raw, dict):
            raise RegionConfigError("each entry in 'regions' must be a table")

        region_id = _require_str(entry_raw, "region_id")
        if not _REGION_ID_RE.fullmatch(region_id):
            raise RegionConfigError(
                f"region_id must match {_REGION_ID_RE.pattern!r}: {region_id!r}",
                region_id=region_id,
            )
        if region_id in seen_region_ids:
            raise RegionConfigError(f"duplicate region_id {region_id!r}", region_id=region_id)
        seen_region_ids.add(region_id)

        chrom = _require_str(entry_raw, "chrom", region_id=region_id)

        start = _require_int(entry_raw, "start", region_id=region_id)
        if start <= 0:
            raise RegionConfigError(f"field 'start' must be > 0, got {start}", region_id=region_id)

        end = _require_int(entry_raw, "end", region_id=region_id)
        if end < start:
            raise RegionConfigError(
                f"field 'end' ({end}) must be >= field 'start' ({start})", region_id=region_id
            )

        title = _require_str(entry_raw, "title", region_id=region_id)

        output_filename = _require_safe_output_filename(entry_raw, region_id)
        if output_filename in seen_output_filenames:
            raise RegionConfigError(
                f"duplicate output_filename {output_filename!r}", region_id=region_id
            )
        seen_output_filenames.add(output_filename)

        regions.append(
            RegionSpec(
                region_id=region_id,
                chrom=chrom,
                start=start,
                end=end,
                title=title,
                output_filename=output_filename,
            )
        )

    return RegionConfig(
        config_schema_version=schema_version,
        dataset_id=dataset_id,
        regions=tuple(regions),
    )
