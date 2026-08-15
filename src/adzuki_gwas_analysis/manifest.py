"""Machine-readable manifest describing the 6 Dryad GWAS summary-statistics datasets.

The manifest is the single source of truth for what a "correct" input looks
like: which archive it came from, what each of the 6 files is (reference,
trait, expected checksum, expected row count), and what the three p-value
columns mean. Nothing in this module downloads or reads GWAS data files
themselves -- see :mod:`adzuki_gwas_analysis.loader` and
:mod:`adzuki_gwas_analysis.validate` for that.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

from adzuki_gwas_analysis.errors import ManifestError

SUPPORTED_SCHEMA_VERSION = 1


@dataclass(frozen=True, slots=True)
class DryadInfo:
    """Identifies the exact Dryad dataset version this manifest describes."""

    doi: str
    dataset_id: int
    version_id: int
    version_number: int
    publication_doi: str


@dataclass(frozen=True, slots=True)
class ArchiveInfo:
    """The single ZIP archive all 6 dataset files are extracted from."""

    filename: str
    size_bytes: int
    sha256: str


@dataclass(frozen=True, slots=True)
class PValueColumns:
    """What each of the three p-value columns present in every dataset file means."""

    p_wald: str
    pval: str
    p_score: str
    primary: str


@dataclass(frozen=True, slots=True)
class DatasetEntry:
    """One (reference, trait) GWAS summary-statistics file extracted from the archive."""

    dataset_id: str
    reference: str
    trait: str
    member_filename: str
    member_sha256: str
    row_count: int


@dataclass(frozen=True, slots=True)
class Manifest:
    """The full, validated manifest: archive provenance plus all dataset entries."""

    schema_version: int
    dryad: DryadInfo
    archive: ArchiveInfo
    pvalue_columns: PValueColumns
    datasets: tuple[DatasetEntry, ...]

    def get(self, dataset_id: str) -> DatasetEntry:
        """Return the entry with the given ``dataset_id`` or raise ``KeyError``."""
        for entry in self.datasets:
            if entry.dataset_id == dataset_id:
                return entry
        raise KeyError(dataset_id)


def _require(table: dict[str, object], key: str, *, dataset_id: str | None = None) -> object:
    if key not in table:
        raise ManifestError(f"missing required field {key!r}", dataset_id=dataset_id)
    return table[key]


def _require_str(table: dict[str, object], key: str, *, dataset_id: str | None = None) -> str:
    value = _require(table, key, dataset_id=dataset_id)
    if not isinstance(value, str) or not value:
        raise ManifestError(f"field {key!r} must be a non-empty string", dataset_id=dataset_id)
    return value


def _require_int(table: dict[str, object], key: str, *, dataset_id: str | None = None) -> int:
    value = _require(table, key, dataset_id=dataset_id)
    if not isinstance(value, int) or isinstance(value, bool):
        raise ManifestError(f"field {key!r} must be an integer", dataset_id=dataset_id)
    return value


def load_manifest(path: str | Path) -> Manifest:
    """Load and structurally validate the manifest at ``path``.

    Raises :class:`~adzuki_gwas_analysis.errors.ManifestError` if the manifest
    is missing required fields, declares an unsupported schema version, or
    contains duplicate ``dataset_id`` entries.
    """
    manifest_path = Path(path)
    if not manifest_path.is_file():
        raise ManifestError(f"manifest file not found: {manifest_path}")

    with manifest_path.open("rb") as fh:
        raw = tomllib.load(fh)

    schema_version = _require_int(raw, "schema_version")
    if schema_version != SUPPORTED_SCHEMA_VERSION:
        raise ManifestError(
            f"unsupported schema_version {schema_version}; "
            f"this package supports version {SUPPORTED_SCHEMA_VERSION}"
        )

    dryad_raw = _require(raw, "dryad")
    if not isinstance(dryad_raw, dict):
        raise ManifestError("field 'dryad' must be a table")
    dryad = DryadInfo(
        doi=_require_str(dryad_raw, "doi"),
        dataset_id=_require_int(dryad_raw, "dataset_id"),
        version_id=_require_int(dryad_raw, "version_id"),
        version_number=_require_int(dryad_raw, "version_number"),
        publication_doi=_require_str(dryad_raw, "publication_doi"),
    )

    archive_raw = _require(raw, "archive")
    if not isinstance(archive_raw, dict):
        raise ManifestError("field 'archive' must be a table")
    archive = ArchiveInfo(
        filename=_require_str(archive_raw, "filename"),
        size_bytes=_require_int(archive_raw, "size_bytes"),
        sha256=_require_str(archive_raw, "sha256"),
    )

    pvalue_raw = _require(raw, "pvalue_columns")
    if not isinstance(pvalue_raw, dict):
        raise ManifestError("field 'pvalue_columns' must be a table")
    pvalue_columns = PValueColumns(
        p_wald=_require_str(pvalue_raw, "p_wald"),
        pval=_require_str(pvalue_raw, "pval"),
        p_score=_require_str(pvalue_raw, "p_score"),
        primary=_require_str(pvalue_raw, "primary"),
    )
    if pvalue_columns.primary not in {"p_wald", "pval", "p_score"}:
        raise ManifestError(
            f"pvalue_columns.primary must be one of 'p_wald', 'pval', 'p_score', "
            f"got {pvalue_columns.primary!r}"
        )

    datasets_raw = _require(raw, "datasets")
    if not isinstance(datasets_raw, list) or not datasets_raw:
        raise ManifestError("field 'datasets' must be a non-empty array of tables")

    seen_dataset_ids: set[str] = set()
    seen_member_filenames: set[str] = set()
    datasets: list[DatasetEntry] = []
    for entry_raw in datasets_raw:
        if not isinstance(entry_raw, dict):
            raise ManifestError("each entry in 'datasets' must be a table")
        dataset_id = _require_str(entry_raw, "dataset_id")
        if dataset_id in seen_dataset_ids:
            raise ManifestError("duplicate dataset_id in manifest", dataset_id=dataset_id)
        seen_dataset_ids.add(dataset_id)

        member_filename = _require_str(entry_raw, "member_filename", dataset_id=dataset_id)
        if member_filename in seen_member_filenames:
            raise ManifestError(
                f"duplicate member_filename {member_filename!r} in manifest",
                dataset_id=dataset_id,
            )
        seen_member_filenames.add(member_filename)

        reference = _require_str(entry_raw, "reference", dataset_id=dataset_id)
        if reference not in {"Miyagi", "Shumari"}:
            raise ManifestError(
                f"reference must be 'Miyagi' or 'Shumari', got {reference!r}",
                dataset_id=dataset_id,
            )

        datasets.append(
            DatasetEntry(
                dataset_id=dataset_id,
                reference=reference,
                trait=_require_str(entry_raw, "trait", dataset_id=dataset_id),
                member_filename=member_filename,
                member_sha256=_require_str(entry_raw, "member_sha256", dataset_id=dataset_id),
                row_count=_require_int(entry_raw, "row_count", dataset_id=dataset_id),
            )
        )

    return Manifest(
        schema_version=schema_version,
        dryad=dryad,
        archive=archive,
        pvalue_columns=pvalue_columns,
        datasets=tuple(datasets),
    )
