"""Diagnosable exceptions raised while loading manifests or validating GWAS input files.

Every exception in this module carries enough structured context (dataset id,
file name, row number, column name) for a caller to locate the offending
record without re-reading the whole file.
"""

from __future__ import annotations


class GwasContractError(Exception):
    """Base class for every input-contract violation raised by this package."""


class ManifestError(GwasContractError):
    """The manifest itself is structurally invalid (independent of any data file)."""

    def __init__(self, message: str, *, dataset_id: str | None = None) -> None:
        self.dataset_id = dataset_id
        prefix = f"[dataset={dataset_id}] " if dataset_id else ""
        super().__init__(f"{prefix}{message}")


class ArchiveChecksumMismatchError(GwasContractError):
    """The downloaded archive's SHA-256 does not match the manifest's recorded value."""

    def __init__(self, *, path: str, expected_sha256: str, actual_sha256: str) -> None:
        self.path = path
        self.expected_sha256 = expected_sha256
        self.actual_sha256 = actual_sha256
        super().__init__(
            f"archive checksum mismatch for {path!r}: "
            f"expected sha256={expected_sha256}, got sha256={actual_sha256}"
        )


class DatasetFileMissingError(GwasContractError):
    """A manifest-declared dataset file could not be found in the data directory."""

    def __init__(self, *, dataset_id: str, path: str) -> None:
        self.dataset_id = dataset_id
        self.path = path
        super().__init__(f"[dataset={dataset_id}] file not found: {path!r}")


class ChecksumMismatchError(GwasContractError):
    """A dataset file's SHA-256 does not match the manifest's recorded value."""

    def __init__(
        self, *, dataset_id: str, path: str, expected_sha256: str, actual_sha256: str
    ) -> None:
        self.dataset_id = dataset_id
        self.path = path
        self.expected_sha256 = expected_sha256
        self.actual_sha256 = actual_sha256
        super().__init__(
            f"[dataset={dataset_id}] checksum mismatch for {path!r}: "
            f"expected sha256={expected_sha256}, got sha256={actual_sha256}"
        )


class RowCountMismatchError(GwasContractError):
    """A dataset file's data-row count does not match the manifest's recorded value."""

    def __init__(
        self, *, dataset_id: str, path: str, expected_row_count: int, actual_row_count: int
    ) -> None:
        self.dataset_id = dataset_id
        self.path = path
        self.expected_row_count = expected_row_count
        self.actual_row_count = actual_row_count
        super().__init__(
            f"[dataset={dataset_id}] row count mismatch for {path!r}: "
            f"expected {expected_row_count}, got {actual_row_count}"
        )


class SchemaError(GwasContractError):
    """A dataset file's header does not match the required column contract."""

    def __init__(self, *, dataset_id: str, path: str, reason: str) -> None:
        self.dataset_id = dataset_id
        self.path = path
        super().__init__(f"[dataset={dataset_id}] schema error in {path!r}: {reason}")


class RowValidationError(GwasContractError):
    """A single data row violates a column-level constraint.

    ``row_number`` is 1-indexed and counts data rows only (the header is not
    row 1), matching the row a caller would find by skipping the header and
    counting from the first data line.
    """

    def __init__(
        self,
        *,
        dataset_id: str,
        path: str,
        row_number: int,
        column: str,
        value: str,
        reason: str,
    ) -> None:
        self.dataset_id = dataset_id
        self.path = path
        self.row_number = row_number
        self.column = column
        self.value = value
        super().__init__(
            f"[dataset={dataset_id}] {path!r} row {row_number}, column {column!r} "
            f"(value={value!r}): {reason}"
        )


class UnknownDatasetIdError(GwasContractError):
    """A caller asked for a ``dataset_id`` that does not exist in the manifest."""

    def __init__(self, *, dataset_id: str, available: tuple[str, ...]) -> None:
        self.dataset_id = dataset_id
        self.available = available
        super().__init__(
            f"unknown dataset_id {dataset_id!r}; manifest defines: {', '.join(available)}"
        )


class DatasetValidationFailedError(GwasContractError):
    """A dataset failed :func:`~adzuki_gwas_analysis.validate.validate_dataset`.

    Raised by analysis entry points so that a failed pre-analysis validation
    always stops before any plot or TSV is produced, rather than being
    silently ignored by a caller that forgets to check ``ValidationResult.success``.
    """

    def __init__(self, *, dataset_id: str, reason: str) -> None:
        self.dataset_id = dataset_id
        super().__init__(f"[dataset={dataset_id}] validation failed, no output produced: {reason}")


class RegionConfigError(GwasContractError):
    """A post-hoc visualization region config file is structurally invalid."""

    def __init__(self, message: str, *, region_id: str | None = None) -> None:
        self.region_id = region_id
        prefix = f"[region={region_id}] " if region_id else ""
        super().__init__(f"{prefix}{message}")


class EmptyRegionError(GwasContractError):
    """A configured region contains zero variants in the validated dataset."""

    def __init__(
        self, *, dataset_id: str, region_id: str, chrom: str, start: int, end: int
    ) -> None:
        self.dataset_id = dataset_id
        self.region_id = region_id
        self.chrom = chrom
        self.start = start
        self.end = end
        super().__init__(
            f"[dataset={dataset_id}] region {region_id!r} ({chrom}:{start}-{end}) "
            f"contains no variants"
        )


class LoadedPvalueCountMismatchError(GwasContractError):
    """The number of primary p-values loaded for diagnostics does not match validation's count.

    Raised by :func:`adzuki_gwas_analysis.analysis.pipeline.run_diagnostics` as a
    data-integrity guard: it is not expected to trigger on a file that
    :func:`~adzuki_gwas_analysis.validate.validate_dataset` has just approved. If it
    does, the file changed between validation and loading, or a loader bug exists --
    either way, ``manifest.toml``'s declared ``row_count`` must never be substituted
    for the actually-loaded count without this check.
    """

    def __init__(self, *, dataset_id: str, validated_row_count: int, loaded_count: int) -> None:
        self.dataset_id = dataset_id
        self.validated_row_count = validated_row_count
        self.loaded_count = loaded_count
        super().__init__(
            f"[dataset={dataset_id}] loaded {loaded_count} primary p-value(s) but "
            f"validation counted {validated_row_count} row(s); refusing to compute "
            f"diagnostics against a mismatched count"
        )


class DuplicateVariantError(GwasContractError):
    """The same (chr, pos, allele0, allele1) key appears in more than one row."""

    def __init__(
        self,
        *,
        dataset_id: str,
        path: str,
        key: tuple[str, int, str, str],
        first_row_number: int,
        duplicate_row_number: int,
    ) -> None:
        self.dataset_id = dataset_id
        self.path = path
        self.key = key
        self.first_row_number = first_row_number
        self.duplicate_row_number = duplicate_row_number
        super().__init__(
            f"[dataset={dataset_id}] {path!r}: duplicate variant key {key!r} "
            f"at row {duplicate_row_number} (first seen at row {first_row_number})"
        )
