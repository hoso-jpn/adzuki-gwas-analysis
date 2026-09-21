"""Read-only validation of a candidate-enabled ``batch`` output directory.

This is the sole gate :mod:`adzuki_gwas_analysis.analysis.report` runs through before
generating anything: it never re-validates the original ``.assoc.txt`` files, never
recomputes Bonferroni/BH/lambda_GC, and never re-clusters candidates -- it only reads
``batch_summary.tsv`` and the per-dataset artifact files ``run_batch(...,
clustering_distance=...)`` already wrote, and checks that they still agree with each
other. Every file this module opens is opened for reading only; nothing here ever writes
to ``--analysis-dir``.

``--analysis-dir`` must be a **candidate-enabled** batch output
(``batch_summary.tsv``'s ``schema_version == 2``, i.e. ``run_batch`` was called with an
explicit ``clustering_distance``) -- a schema-v1 (no-candidates) batch output is rejected
with :class:`~adzuki_gwas_analysis.errors.ReportRequiresCandidateEnabledBatchError` rather
than silently proceeding with an implicit clustering distance. There is no repository
default for that parameter (see
:mod:`adzuki_gwas_analysis.analysis.candidates`'s module docstring), so a report can never
invent one on the caller's behalf.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

import pandas as pd

from adzuki_gwas_analysis.analysis.batch import BATCH_SUMMARY_SCHEMA_VERSION_WITH_CANDIDATES
from adzuki_gwas_analysis.errors import (
    ReportRequiresCandidateEnabledBatchError,
    ReportSourceInconsistentError,
    ReportSourcePathUnsafeError,
)

#: The 7 per-dataset artifact-path columns batch_summary.tsv must carry for every row, and
#: the attribute name on DatasetArtifactPaths each one resolves to.
_ARTIFACT_PATH_COLUMNS: tuple[str, ...] = (
    "manhattan_path",
    "qq_path",
    "diagnostics_path",
    "significant_variants_path",
    "association_peaks_path",
    "candidate_snps_path",
    "candidate_ranking_path",
)


@dataclass(frozen=True, slots=True)
class DatasetArtifactPaths:
    """Resolved, safety-checked absolute paths for one dataset's 7 artifact files."""

    manhattan: Path
    qq: Path
    diagnostics: Path
    significant_variants: Path
    association_peaks: Path
    candidate_snps: Path
    candidate_ranking: Path


@dataclass(frozen=True, slots=True)
class ValidatedDataset:
    """One ``batch_summary.tsv`` row, its artifact paths, and its cross-checked counts."""

    dataset_id: str
    reference: str
    trait: str
    source_sha256: str
    pvalue_column: str
    n_tests: int
    bonferroni_discoveries: int
    bh_discoveries: int
    lambda_gc: float
    lambda_gc_df: int
    n_signals: int
    n_candidates: int
    relative_paths: dict[str, str]
    paths: DatasetArtifactPaths
    #: Already read and cross-checked above -- report content generation reuses these
    #: in-memory tables rather than re-reading the same 3 small files a second time.
    association_peaks_df: pd.DataFrame
    candidate_snps_df: pd.DataFrame
    candidate_ranking_df: pd.DataFrame


@dataclass(frozen=True, slots=True)
class ValidatedAnalysisDir:
    """A fully cross-checked candidate-enabled batch output, ready for report generation."""

    analysis_dir: Path
    batch_summary_path: Path
    alpha: float
    fdr_level: float
    visualization_threshold: float
    clustering_distance: int
    datasets: tuple[ValidatedDataset, ...]
    analysis_generation: dict[str, object] | None = None


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _require_safe_relative_path(
    analysis_dir: Path, relative_path: object, *, dataset_id: str | None
) -> Path:
    if not isinstance(relative_path, str) or not relative_path:
        raise ReportSourcePathUnsafeError(
            dataset_id=dataset_id,
            relative_path=str(relative_path),
            reason="must be a non-empty relative path string",
        )
    pure = PurePosixPath(relative_path)
    if pure.is_absolute() or relative_path.startswith(("/", "\\")):
        raise ReportSourcePathUnsafeError(
            dataset_id=dataset_id, relative_path=relative_path, reason="must not be absolute"
        )
    if ".." in pure.parts or "\\" in relative_path:
        raise ReportSourcePathUnsafeError(
            dataset_id=dataset_id,
            relative_path=relative_path,
            reason="must not contain '..' traversal or backslash path separators",
        )

    candidate = analysis_dir / relative_path
    if candidate.is_symlink():
        raise ReportSourcePathUnsafeError(
            dataset_id=dataset_id,
            relative_path=relative_path,
            reason="must not be a symlink (broken or not)",
        )

    resolved_root = analysis_dir.resolve()
    resolved = candidate.resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError:
        raise ReportSourcePathUnsafeError(
            dataset_id=dataset_id,
            relative_path=relative_path,
            reason="resolves outside --analysis-dir",
        ) from None

    if not candidate.is_file():
        raise ReportSourcePathUnsafeError(
            dataset_id=dataset_id,
            relative_path=relative_path,
            reason="does not exist as a plain regular file",
        )
    return candidate


def _require_uniform_column(table: pd.DataFrame, column: str) -> float:
    if column not in table.columns:
        raise ReportSourceInconsistentError(
            dataset_id=None, reason=f"batch_summary.tsv is missing required column {column!r}"
        )
    values = table[column].unique().tolist()
    if len(values) != 1:
        raise ReportSourceInconsistentError(
            dataset_id=None,
            reason=(
                f"batch_summary.tsv column {column!r} is not uniform across rows "
                f"(a report describes one batch run with one set of parameters): {values!r}"
            ),
        )
    value = values[0]
    if not isinstance(value, int | float):
        raise ReportSourceInconsistentError(
            dataset_id=None,
            reason=f"batch_summary.tsv column {column!r} must be numeric, got {value!r}",
        )
    return float(value)


def _require_close(
    *, dataset_id: str, field: str, batch_summary_value: float, diagnostics_value: float
) -> None:
    if not math.isclose(batch_summary_value, diagnostics_value, rel_tol=1e-9, abs_tol=1e-12):
        raise ReportSourceInconsistentError(
            dataset_id=dataset_id,
            reason=(
                f"{field} disagrees between batch_summary.tsv ({batch_summary_value!r}) and "
                f"statistical_diagnostics.tsv ({diagnostics_value!r})"
            ),
        )


def _require_equal(*, dataset_id: str, field: str, expected: object, actual: object) -> None:
    if expected != actual:
        raise ReportSourceInconsistentError(
            dataset_id=dataset_id,
            reason=f"{field} disagrees: batch_summary.tsv says {expected!r}, found {actual!r}",
        )


def _validate_one_dataset(analysis_dir: Path, row: pd.Series) -> ValidatedDataset:
    dataset_id = str(row["dataset_id"])
    reference = str(row["reference"])
    trait = str(row["trait"])

    relative_paths = {column: row[column] for column in _ARTIFACT_PATH_COLUMNS}
    for column in _ARTIFACT_PATH_COLUMNS:
        if column not in row.index:
            raise ReportSourceInconsistentError(
                dataset_id=dataset_id,
                reason=f"batch_summary.tsv is missing required column {column!r}",
            )
    resolved = {
        column: _require_safe_relative_path(analysis_dir, row[column], dataset_id=dataset_id)
        for column in _ARTIFACT_PATH_COLUMNS
    }
    paths = DatasetArtifactPaths(
        manhattan=resolved["manhattan_path"],
        qq=resolved["qq_path"],
        diagnostics=resolved["diagnostics_path"],
        significant_variants=resolved["significant_variants_path"],
        association_peaks=resolved["association_peaks_path"],
        candidate_snps=resolved["candidate_snps_path"],
        candidate_ranking=resolved["candidate_ranking_path"],
    )

    n_signals = int(row["n_signals"])
    n_candidates = int(row["n_candidates"])

    diagnostics = _read_tsv(paths.diagnostics)
    if len(diagnostics) != 1:
        raise ReportSourceInconsistentError(
            dataset_id=dataset_id,
            reason=f"statistical_diagnostics.tsv has {len(diagnostics)} rows, expected exactly 1",
        )
    diag_row = diagnostics.iloc[0]
    _require_equal(
        dataset_id=dataset_id,
        field="dataset_id",
        expected=dataset_id,
        actual=str(diag_row["dataset_id"]),
    )
    _require_equal(
        dataset_id=dataset_id,
        field="source_sha256",
        expected=str(row["source_sha256"]),
        actual=str(diag_row["source_sha256"]),
    )
    _require_equal(
        dataset_id=dataset_id,
        field="n_tests",
        expected=int(row["n_tests"]),
        actual=int(diag_row["n_tests"]),
    )
    _require_equal(
        dataset_id=dataset_id,
        field="bonferroni_discoveries",
        expected=int(row["bonferroni_discoveries"]),
        actual=int(diag_row["bonferroni_discoveries"]),
    )
    _require_equal(
        dataset_id=dataset_id,
        field="bh_discoveries",
        expected=int(row["bh_discoveries"]),
        actual=int(diag_row["bh_discoveries"]),
    )
    _require_close(
        dataset_id=dataset_id,
        field="lambda_gc",
        batch_summary_value=float(row["lambda_gc"]),
        diagnostics_value=float(diag_row["lambda_gc"]),
    )

    significant_variants = _read_tsv(paths.significant_variants)
    _require_equal(
        dataset_id=dataset_id,
        field="len(significant_variants.tsv)",
        expected=n_candidates,
        actual=len(significant_variants),
    )

    association_peaks = _read_tsv(paths.association_peaks)
    _require_equal(
        dataset_id=dataset_id,
        field="len(association_peaks.tsv)",
        expected=n_signals,
        actual=len(association_peaks),
    )
    if n_signals > 0:
        for column in ("dataset_id", "reference", "trait"):
            observed = set(association_peaks[column].astype(str).unique())
            expected_value = {"dataset_id": dataset_id, "reference": reference, "trait": trait}[
                column
            ]
            if observed != {expected_value}:
                raise ReportSourceInconsistentError(
                    dataset_id=dataset_id,
                    reason=(
                        f"association_peaks.tsv column {column!r} contains {sorted(observed)}, "
                        f"expected only {expected_value!r}"
                    ),
                )
        summed = int(association_peaks["n_significant_variants"].sum())
        _require_equal(
            dataset_id=dataset_id,
            field="sum(association_peaks.n_significant_variants)",
            expected=n_candidates,
            actual=summed,
        )
    signal_ids = set(association_peaks["signal_id"].astype(str)) if n_signals > 0 else set()

    candidate_snps = _read_tsv(paths.candidate_snps)
    _require_equal(
        dataset_id=dataset_id,
        field="len(candidate_snps.tsv)",
        expected=n_candidates,
        actual=len(candidate_snps),
    )

    candidate_ranking = _read_tsv(paths.candidate_ranking)
    _require_equal(
        dataset_id=dataset_id,
        field="len(candidate_ranking.tsv)",
        expected=n_candidates,
        actual=len(candidate_ranking),
    )

    for table_name, table in (
        ("candidate_snps.tsv", candidate_snps),
        ("candidate_ranking.tsv", candidate_ranking),
    ):
        if n_candidates == 0:
            continue
        for column in ("dataset_id", "reference", "trait"):
            observed = set(table[column].astype(str).unique())
            expected_value = {"dataset_id": dataset_id, "reference": reference, "trait": trait}[
                column
            ]
            if observed != {expected_value}:
                raise ReportSourceInconsistentError(
                    dataset_id=dataset_id,
                    reason=(
                        f"{table_name} column {column!r} contains {sorted(observed)}, "
                        f"expected only {expected_value!r}"
                    ),
                )
        table_signal_ids = set(table["signal_id"].astype(str))
        unknown = table_signal_ids - signal_ids
        if unknown:
            raise ReportSourceInconsistentError(
                dataset_id=dataset_id,
                reason=(
                    f"{table_name} references signal_id(s) {sorted(unknown)} not present in "
                    f"association_peaks.tsv"
                ),
            )

    if n_candidates > 0:
        ranks = sorted(int(r) for r in candidate_ranking["candidate_rank"])
        if ranks != list(range(1, n_candidates + 1)):
            raise ReportSourceInconsistentError(
                dataset_id=dataset_id,
                reason=(
                    f"candidate_ranking.tsv's candidate_rank values are not a dense "
                    f"1..{n_candidates} sequence: {ranks!r}"
                ),
            )

    return ValidatedDataset(
        dataset_id=dataset_id,
        reference=reference,
        trait=trait,
        source_sha256=str(row["source_sha256"]),
        pvalue_column=str(row["pvalue_column"]),
        n_tests=int(row["n_tests"]),
        bonferroni_discoveries=int(row["bonferroni_discoveries"]),
        bh_discoveries=int(row["bh_discoveries"]),
        lambda_gc=float(row["lambda_gc"]),
        lambda_gc_df=int(row["lambda_gc_df"]),
        n_signals=n_signals,
        n_candidates=n_candidates,
        relative_paths={k: str(v) for k, v in relative_paths.items()},
        paths=paths,
        association_peaks_df=association_peaks,
        candidate_snps_df=candidate_snps,
        candidate_ranking_df=candidate_ranking,
    )


def validate_analysis_dir(
    analysis_dir: Path, *, require_provenance: bool = False
) -> ValidatedAnalysisDir:
    """Read and cross-check every artifact under a candidate-enabled ``batch`` output.

    Raises :class:`~adzuki_gwas_analysis.errors.ReportRequiresCandidateEnabledBatchError`
    if ``batch_summary.tsv``'s ``schema_version`` is not the candidate-enabled value (2),
    :class:`~adzuki_gwas_analysis.errors.ReportSourcePathUnsafeError` for any unsafe
    artifact path, and :class:`~adzuki_gwas_analysis.errors.ReportSourceInconsistentError`
    for any cross-artifact disagreement -- including a dataset with zero candidates, which
    is a normal, fully-checked state (header-only files, ``n_candidates=0``), not an error.
    """
    analysis_dir = Path(analysis_dir)
    from adzuki_gwas_analysis.provenance import validate_provenance

    generation = validate_provenance(analysis_dir, required=require_provenance)
    batch_summary_path = _require_safe_relative_path(
        analysis_dir, "batch_summary.tsv", dataset_id=None
    )
    summary = _read_tsv(batch_summary_path)

    if "schema_version" in summary.columns:
        schema_versions = summary["schema_version"].unique().tolist()
    else:
        schema_versions = []
    if schema_versions != [BATCH_SUMMARY_SCHEMA_VERSION_WITH_CANDIDATES]:
        found = schema_versions[0] if len(schema_versions) == 1 else schema_versions
        raise ReportRequiresCandidateEnabledBatchError(
            analysis_dir=str(analysis_dir), found_schema_version=found
        )

    alpha = float(_require_uniform_column(summary, "alpha"))
    fdr_level = float(_require_uniform_column(summary, "fdr_level"))
    visualization_threshold = float(_require_uniform_column(summary, "visualization_threshold"))
    clustering_distance = int(_require_uniform_column(summary, "clustering_distance"))

    datasets = tuple(_validate_one_dataset(analysis_dir, row) for _, row in summary.iterrows())

    return ValidatedAnalysisDir(
        analysis_dir=analysis_dir,
        batch_summary_path=batch_summary_path,
        alpha=alpha,
        fdr_level=fdr_level,
        visualization_threshold=visualization_threshold,
        clustering_distance=clustering_distance,
        datasets=datasets,
        analysis_generation=generation,
    )
