"""Machine-readable manifest describing the 6 Dryad GWAS summary-statistics datasets.

The manifest is the single source of truth for what a "correct" input looks
like: which archive it came from, what each of the 6 files is (reference,
trait, expected checksum, expected row count), and what the three p-value
columns mean. Nothing in this module downloads or reads GWAS data files
themselves -- see :mod:`adzuki_gwas_analysis.loader` and
:mod:`adzuki_gwas_analysis.validate` for that.

Schema v1's dataset contract -- exactly these 6 (reference, trait) pairs,
each with one canonical ``dataset_id`` and ``member_filename`` -- is
enforced here as data (:data:`SCHEMA_V1_DATASETS`), not just documented in
prose. A manifest with 1, 5, or 7 datasets, an unrecognized trait, a
reference/trait pair that doesn't match its own ``dataset_id``, or a
non-canonical ``member_filename`` all fail :func:`load_manifest` (PR #2
review, P1-1).
"""

from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from adzuki_gwas_analysis.errors import ManifestError

SUPPORTED_SCHEMA_VERSION = 1

REQUIRED_REFERENCES: tuple[str, ...] = ("Miyagi", "Shumari")
REQUIRED_TRAITS: tuple[str, ...] = (
    "water_permeability",
    "red_seedcoat",
    "mottled_black_seedcoat",
)

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True, slots=True)
class CanonicalDatasetSpec:
    """The one correct ``dataset_id``/``member_filename`` for a (reference, trait) pair."""

    dataset_id: str
    member_filename: str


# Schema v1's closed dataset contract: exactly these 6 (reference, trait)
# pairs exist, each with exactly one canonical dataset_id and
# member_filename -- "6 datasets, no more, no less" is Issue #1's own
# explicit deliverable, so it is data here, not just prose.
SCHEMA_V1_DATASETS: dict[tuple[str, str], CanonicalDatasetSpec] = {
    ("Miyagi", "water_permeability"): CanonicalDatasetSpec(
        dataset_id="miyagi_water_permeability",
        member_filename="mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt",
    ),
    ("Miyagi", "red_seedcoat"): CanonicalDatasetSpec(
        dataset_id="miyagi_red_seedcoat",
        member_filename="mapped_to_Miyagi_red_seedcoat.maf_0.05.assoc.txt",
    ),
    ("Miyagi", "mottled_black_seedcoat"): CanonicalDatasetSpec(
        dataset_id="miyagi_mottled_black_seedcoat",
        member_filename="mapped_to_Miyagi_mottled_black_seedcoat.maf_0.05.assoc.txt",
    ),
    ("Shumari", "water_permeability"): CanonicalDatasetSpec(
        dataset_id="shumari_water_permeability",
        member_filename="mapped_to_Shumari_water_permeability.maf_0.05.assoc.txt",
    ),
    ("Shumari", "red_seedcoat"): CanonicalDatasetSpec(
        dataset_id="shumari_red_seedcoat",
        member_filename="mapped_to_Shumari_red_seedcoat.maf_0.05.assoc.txt",
    ),
    ("Shumari", "mottled_black_seedcoat"): CanonicalDatasetSpec(
        dataset_id="shumari_mottled_black_seedcoat",
        member_filename="mapped_to_Shumari_mottled_black_seedcoat.maf_0.05.assoc.txt",
    ),
}


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


def _require_positive_int(
    table: dict[str, object], key: str, *, dataset_id: str | None = None
) -> int:
    value = _require_int(table, key, dataset_id=dataset_id)
    if value <= 0:
        raise ManifestError(
            f"field {key!r} must be a positive integer, got {value}", dataset_id=dataset_id
        )
    return value


def _require_sha256(table: dict[str, object], key: str, *, dataset_id: str | None = None) -> str:
    value = _require_str(table, key, dataset_id=dataset_id)
    if not _SHA256_RE.fullmatch(value):
        raise ManifestError(
            f"field {key!r} must be a lowercase 64-character hexadecimal SHA-256 digest, "
            f"got {value!r}",
            dataset_id=dataset_id,
        )
    return value


def _require_simple_filename(
    table: dict[str, object],
    key: str,
    *,
    expected_suffix: str,
    dataset_id: str | None = None,
) -> str:
    """Require ``table[key]`` to be a plain, safe basename ending in ``expected_suffix``.

    Rejects absolute paths, ``..`` traversal, any path separator, and the
    wrong file extension -- a manifest is not allowed to point outside the
    data directory it will be resolved against.
    """
    value = _require_str(table, key, dataset_id=dataset_id)
    if PurePosixPath(value).is_absolute() or value.startswith(("/", "\\")):
        raise ManifestError(
            f"field {key!r} must not be an absolute path: {value!r}", dataset_id=dataset_id
        )
    if "/" in value or "\\" in value:
        raise ManifestError(
            f"field {key!r} must be a plain filename with no path separators: {value!r}",
            dataset_id=dataset_id,
        )
    if ".." in value:
        raise ManifestError(
            f"field {key!r} must not contain '..': {value!r}", dataset_id=dataset_id
        )
    if not value.endswith(expected_suffix):
        raise ManifestError(
            f"field {key!r} must end with {expected_suffix!r}: {value!r}", dataset_id=dataset_id
        )
    return value


def load_manifest(path: str | Path) -> Manifest:
    """Load and structurally validate the manifest at ``path``.

    Raises :class:`~adzuki_gwas_analysis.errors.ManifestError` for: a
    missing file; invalid TOML syntax; an unsupported ``schema_version``;
    any field that is missing, the wrong type, an unsafe/wrongly-suffixed
    filename, a non-positive size/row-count, or a malformed SHA-256; a
    ``pvalue_columns.primary`` other than ``"pval"``; and, per schema v1's
    closed dataset contract (:data:`SCHEMA_V1_DATASETS`): anything other
    than exactly 6 datasets, an unrecognized reference or trait, a
    duplicate (reference, trait) pair, or a ``dataset_id``/
    ``member_filename`` that doesn't match its canonical value for that
    (reference, trait) pair.
    """
    manifest_path = Path(path)
    if not manifest_path.is_file():
        raise ManifestError(f"manifest file not found: {manifest_path}")

    try:
        with manifest_path.open("rb") as fh:
            raw = tomllib.load(fh)
    except tomllib.TOMLDecodeError as exc:
        raise ManifestError(f"manifest is not valid TOML: {exc}") from exc
    except OSError as exc:
        raise ManifestError(f"could not read manifest file: {exc}") from exc

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
        dataset_id=_require_positive_int(dryad_raw, "dataset_id"),
        version_id=_require_positive_int(dryad_raw, "version_id"),
        version_number=_require_positive_int(dryad_raw, "version_number"),
        publication_doi=_require_str(dryad_raw, "publication_doi"),
    )

    archive_raw = _require(raw, "archive")
    if not isinstance(archive_raw, dict):
        raise ManifestError("field 'archive' must be a table")
    archive = ArchiveInfo(
        filename=_require_simple_filename(archive_raw, "filename", expected_suffix=".zip"),
        size_bytes=_require_positive_int(archive_raw, "size_bytes"),
        sha256=_require_sha256(archive_raw, "sha256"),
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
    if pvalue_columns.primary != "pval":
        raise ManifestError(
            f"schema v1 requires pvalue_columns.primary == 'pval' "
            f"(the likelihood ratio test p-value), got {pvalue_columns.primary!r}"
        )

    datasets_raw = _require(raw, "datasets")
    if not isinstance(datasets_raw, list) or not datasets_raw:
        raise ManifestError("field 'datasets' must be a non-empty array of tables")

    seen_dataset_ids: set[str] = set()
    seen_member_filenames: set[str] = set()
    seen_reference_traits: set[tuple[str, str]] = set()
    datasets: list[DatasetEntry] = []
    for entry_raw in datasets_raw:
        if not isinstance(entry_raw, dict):
            raise ManifestError("each entry in 'datasets' must be a table")
        dataset_id = _require_str(entry_raw, "dataset_id")
        if dataset_id in seen_dataset_ids:
            raise ManifestError("duplicate dataset_id in manifest", dataset_id=dataset_id)
        seen_dataset_ids.add(dataset_id)

        member_filename = _require_simple_filename(
            entry_raw, "member_filename", expected_suffix=".assoc.txt", dataset_id=dataset_id
        )
        if member_filename in seen_member_filenames:
            raise ManifestError(
                f"duplicate member_filename {member_filename!r} in manifest",
                dataset_id=dataset_id,
            )
        seen_member_filenames.add(member_filename)

        reference = _require_str(entry_raw, "reference", dataset_id=dataset_id)
        if reference not in REQUIRED_REFERENCES:
            raise ManifestError(
                f"reference must be one of {REQUIRED_REFERENCES}, got {reference!r}",
                dataset_id=dataset_id,
            )

        trait = _require_str(entry_raw, "trait", dataset_id=dataset_id)
        if trait not in REQUIRED_TRAITS:
            raise ManifestError(
                f"trait must be one of {REQUIRED_TRAITS}, got {trait!r}",
                dataset_id=dataset_id,
            )

        reference_trait = (reference, trait)
        if reference_trait in seen_reference_traits:
            raise ManifestError(
                f"duplicate (reference, trait) combination: {reference_trait!r}",
                dataset_id=dataset_id,
            )
        seen_reference_traits.add(reference_trait)

        canonical = SCHEMA_V1_DATASETS[reference_trait]
        if dataset_id != canonical.dataset_id:
            raise ManifestError(
                f"dataset_id must be {canonical.dataset_id!r} for "
                f"reference={reference!r}, trait={trait!r}; got {dataset_id!r}",
                dataset_id=dataset_id,
            )
        if member_filename != canonical.member_filename:
            raise ManifestError(
                f"member_filename must be {canonical.member_filename!r} for "
                f"reference={reference!r}, trait={trait!r}; got {member_filename!r}",
                dataset_id=dataset_id,
            )

        datasets.append(
            DatasetEntry(
                dataset_id=dataset_id,
                reference=reference,
                trait=trait,
                member_filename=member_filename,
                member_sha256=_require_sha256(entry_raw, "member_sha256", dataset_id=dataset_id),
                row_count=_require_positive_int(entry_raw, "row_count", dataset_id=dataset_id),
            )
        )

    if len(datasets) != len(SCHEMA_V1_DATASETS):
        raise ManifestError(
            f"schema v1 requires exactly {len(SCHEMA_V1_DATASETS)} datasets "
            f"({', '.join(REQUIRED_REFERENCES)} x {', '.join(REQUIRED_TRAITS)}), "
            f"got {len(datasets)}"
        )

    return Manifest(
        schema_version=schema_version,
        dryad=dryad,
        archive=archive,
        pvalue_columns=pvalue_columns,
        datasets=tuple(datasets),
    )
