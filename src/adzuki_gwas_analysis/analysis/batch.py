"""Sequential, single-pass batch orchestration across all 6 manifest-declared datasets.

``manifest.toml`` itself is loaded exactly once, at the start of :func:`run_batch` --
never once per dataset. Each of its 6 declared entries is then validated via
:func:`~adzuki_gwas_analysis.analysis.pipeline.ensure_validated_entry`, which takes that
already-loaded :class:`~adzuki_gwas_analysis.manifest.Manifest` object directly and calls
:func:`~adzuki_gwas_analysis.validate.validate_dataset` without re-parsing the manifest
file (unlike :func:`~adzuki_gwas_analysis.analysis.pipeline.ensure_validated_with_result`,
which the standalone single-dataset ``diagnostics`` subcommand uses, and which does load
the manifest itself since it only ever validates one dataset per call).

For each dataset, in ``manifest.toml``'s own declared order (never re-sorted): resolve the
manifest entry, run schema v1 validation exactly once, load its analysis DataFrame exactly
once, and reuse that one DataFrame for its Manhattan plot, QQ plot, and Bonferroni/BH/
lambda_GC diagnostics -- never re-validating or re-loading the same file for a second
artifact. Only :class:`BatchDatasetResult` (scalars and output-relative POSIX path
strings) survives past one dataset's processing; the loaded DataFrame, its p-value array,
and the diagnostics result's adjusted-p-value/significance arrays are all local to
:func:`_process_one_dataset` and are dropped (eligible for garbage collection) before the
next dataset starts, so peak memory never holds more than one dataset's data, and
``batch_summary.tsv`` never becomes a second, larger multiple-testing family: each row's
``family_scope`` describes that row's own dataset in isolation (see
:func:`adzuki_gwas_analysis.analysis.diagnostics.format_family_scope`), and the 6 datasets'
p-values are never combined, compared, or jointly corrected anywhere in this module.

This module does not reuse :func:`adzuki_gwas_analysis.analysis.pipeline.run_all` (which
also renders regional plots and a top-variant-by-region TSV -- both out of scope for a
batch run across all 6 datasets, since post-hoc regions are defined for
``miyagi_water_permeability`` only) or the standalone ``run_manhattan``/``run_qq``/
``run_diagnostics`` (each of which validates and loads independently; calling three of them
per dataset would validate and load the same file three times).
"""

from __future__ import annotations

import os
import shutil
import tempfile
from dataclasses import dataclass, fields
from pathlib import Path, PurePosixPath

import pandas as pd

from adzuki_gwas_analysis.analysis.chromosomes import compute_manhattan_coordinates
from adzuki_gwas_analysis.analysis.diagnostics import (
    PVALUE_SEMANTICS,
    build_significant_variants_table,
    build_summary_table,
    format_family_scope,
)
from adzuki_gwas_analysis.analysis.loader import load_analysis_frame
from adzuki_gwas_analysis.analysis.pipeline import (
    DEFAULT_ALPHA,
    DEFAULT_FDR_LEVEL,
    DEFAULT_THRESHOLD,
    ValidatedDatasetInfo,
    atomic_write_tsv,
    compute_diagnostics_result,
    ensure_validated_entry,
    validate_threshold,
)
from adzuki_gwas_analysis.analysis.plotting import plot_manhattan, plot_qq
from adzuki_gwas_analysis.analysis.qq import compute_qq_points
from adzuki_gwas_analysis.errors import BatchOutputDirectoryUnsafeError
from adzuki_gwas_analysis.manifest import DatasetEntry, Manifest, load_manifest

#: Bumped only if batch_summary.tsv's column set/meaning changes.
BATCH_SUMMARY_SCHEMA_VERSION = 1

#: Every batch dataset directory holds exactly these 4 artifacts.
ARTIFACTS_PER_DATASET = 4

#: schema v1's 3 trait codes, mapped to the fixed customer-facing label Issue #9 requires.
#: Deliberately built from ``DatasetEntry.trait`` (a manifest-declared value), never by
#: string-splitting ``dataset_id``.
TRAIT_LABELS: dict[str, str] = {
    "water_permeability": "Water Permeability",
    "red_seedcoat": "Red Seed Coat Color",
    "mottled_black_seedcoat": "Mottled Black Seed Coat Color",
}


def format_trait_label(trait: str) -> str:
    """Return the fixed customer-facing label for a manifest ``trait`` code."""
    try:
        return TRAIT_LABELS[trait]
    except KeyError:
        raise ValueError(
            f"no display label configured for trait {trait!r}; expected one of "
            f"{sorted(TRAIT_LABELS)}"
        ) from None


def format_manhattan_title(entry: DatasetEntry) -> str:
    """``"<Trait Label> GWAS (<Reference> reference)"``.

    Built from ``entry.trait``/``entry.reference``, never by string-splitting ``dataset_id``.
    """
    return f"{format_trait_label(entry.trait)} GWAS ({entry.reference} reference)"


def format_qq_title(entry: DatasetEntry) -> str:
    """``"QQ Plot: <Trait Label> GWAS (<Reference> reference)"``."""
    return f"QQ Plot: {format_trait_label(entry.trait)} GWAS ({entry.reference} reference)"


@dataclass(frozen=True, slots=True)
class BatchDatasetResult:
    """One dataset's row in ``batch_summary.tsv`` -- scalars and relative paths only.

    Deliberately holds no DataFrame, p-value array, adjusted-p-value array, significance
    mask, or plot-coordinate array, and no
    :class:`~adzuki_gwas_analysis.analysis.diagnostics.DiagnosticsResult` (which itself
    holds exactly those arrays) -- only what a batch summary row needs once that dataset's
    processing is done.
    """

    dataset_id: str
    reference: str
    trait: str
    source_sha256: str
    pvalue_column: str
    n_tests: int
    visualization_threshold: float
    alpha: float
    bonferroni_threshold: float
    bonferroni_discoveries: int
    fdr_level: float
    bh_raw_p_cutoff: float | None
    bh_discoveries: int
    lambda_gc_df: int
    expected_chi2_median: float
    lambda_gc: float
    manhattan_path: str
    qq_path: str
    diagnostics_path: str
    significant_variants_path: str


@dataclass(frozen=True, slots=True)
class BatchOutcome:
    """Everything produced by one :func:`run_batch` call."""

    output_dir: Path
    datasets: tuple[BatchDatasetResult, ...]
    summary_path: Path
    total_n_tests: int


def _relative_posix(path: Path, base: Path) -> str:
    return str(PurePosixPath(path.relative_to(base).as_posix()))


def _process_one_dataset(
    *,
    manifest: Manifest,
    data_dir: Path,
    entry: DatasetEntry,
    dataset_staging_dir: Path,
    threshold: float,
    alpha: float,
    fdr_level: float,
) -> BatchDatasetResult:
    """Validate once, load once, render both plots, compute diagnostics, write all 4 files.

    Takes the already-loaded ``manifest`` (loaded exactly once, by :func:`run_batch`, for
    the whole 6-dataset run) rather than a ``manifest_path`` -- validating via
    :func:`~adzuki_gwas_analysis.analysis.pipeline.ensure_validated_entry` means this
    function never re-parses ``manifest.toml`` itself. Everything that touches the loaded
    DataFrame or a p-value/adjusted-p-value array is local to this function's stack frame;
    only the returned :class:`BatchDatasetResult` (scalars and path strings) escapes it.
    """
    info: ValidatedDatasetInfo = ensure_validated_entry(
        manifest=manifest, entry=entry, data_dir=data_dir
    )
    variant_df = load_analysis_frame(data_dir / info.entry.member_filename)

    dataset_staging_dir.mkdir(parents=True, exist_ok=True)

    coords = compute_manhattan_coordinates(variant_df)
    manhattan_path = dataset_staging_dir / f"{entry.dataset_id}_manhattan.png"
    plot_manhattan(
        coords,
        variant_df["pval"],
        threshold=threshold,
        output_path=manhattan_path,
        title=format_manhattan_title(entry),
    )

    qq_points = compute_qq_points(variant_df["pval"].to_numpy(dtype="float64"))
    qq_path = dataset_staging_dir / f"{entry.dataset_id}_qq.png"
    plot_qq(qq_points, output_path=qq_path, title=format_qq_title(entry))

    diagnostics_result = compute_diagnostics_result(
        info, variant_df, alpha=alpha, fdr_level=fdr_level
    )
    significant_table = build_significant_variants_table(variant_df, diagnostics_result)

    diagnostics_path = dataset_staging_dir / "statistical_diagnostics.tsv"
    significant_variants_path = dataset_staging_dir / "significant_variants.tsv"
    # This per-dataset statistical_diagnostics.tsv uses the same build_summary_table() as
    # the standalone `diagnostics` subcommand -- batch_summary.tsv (built separately, once
    # all 6 datasets are done) is the batch-specific rollup with its own schema.
    atomic_write_tsv(build_summary_table(diagnostics_result), diagnostics_path)
    atomic_write_tsv(significant_table, significant_variants_path)

    return BatchDatasetResult(
        dataset_id=entry.dataset_id,
        reference=entry.reference,
        trait=entry.trait,
        source_sha256=diagnostics_result.source_sha256,
        pvalue_column=diagnostics_result.pvalue_column,
        n_tests=diagnostics_result.n_tests,
        visualization_threshold=threshold,
        alpha=diagnostics_result.bonferroni.alpha,
        bonferroni_threshold=diagnostics_result.bonferroni.threshold,
        bonferroni_discoveries=diagnostics_result.bonferroni.discoveries,
        fdr_level=diagnostics_result.bh.fdr_level,
        bh_raw_p_cutoff=diagnostics_result.bh.raw_p_cutoff,
        bh_discoveries=diagnostics_result.bh.discoveries,
        lambda_gc_df=diagnostics_result.lambda_gc.df,
        expected_chi2_median=diagnostics_result.lambda_gc.expected_median,
        lambda_gc=diagnostics_result.lambda_gc.lambda_gc,
        manhattan_path=_relative_posix(manhattan_path, dataset_staging_dir.parent),
        qq_path=_relative_posix(qq_path, dataset_staging_dir.parent),
        diagnostics_path=_relative_posix(diagnostics_path, dataset_staging_dir.parent),
        significant_variants_path=_relative_posix(
            significant_variants_path, dataset_staging_dir.parent
        ),
    )


def build_batch_summary_table(results: tuple[BatchDatasetResult, ...]) -> pd.DataFrame:
    """Build ``batch_summary.tsv``: one row per dataset, in the order ``results`` was given.

    ``bh_raw_p_cutoff`` renders as an explicit missing value (``float("nan")``, never a
    Python ``None`` sentinel) when that dataset had zero BH discoveries, matching the
    standalone ``diagnostics`` subcommand's ``statistical_diagnostics.tsv`` convention.
    """
    rows = []
    for result in results:
        raw_p_cutoff = result.bh_raw_p_cutoff
        rows.append(
            {
                "schema_version": BATCH_SUMMARY_SCHEMA_VERSION,
                "dataset_id": result.dataset_id,
                "reference": result.reference,
                "trait": result.trait,
                "source_sha256": result.source_sha256,
                "pvalue_column": result.pvalue_column,
                "pvalue_semantics": PVALUE_SEMANTICS,
                "family_scope": format_family_scope(result.dataset_id),
                "n_tests": result.n_tests,
                "visualization_threshold": result.visualization_threshold,
                "alpha": result.alpha,
                "bonferroni_threshold": result.bonferroni_threshold,
                "bonferroni_discoveries": result.bonferroni_discoveries,
                "fdr_level": result.fdr_level,
                "bh_raw_p_cutoff": raw_p_cutoff if raw_p_cutoff is not None else float("nan"),
                "bh_discoveries": result.bh_discoveries,
                "lambda_gc_df": result.lambda_gc_df,
                "expected_chi2_median": result.expected_chi2_median,
                "lambda_gc": result.lambda_gc,
                "manhattan_path": result.manhattan_path,
                "qq_path": result.qq_path,
                "diagnostics_path": result.diagnostics_path,
                "significant_variants_path": result.significant_variants_path,
            }
        )
    return pd.DataFrame(rows)


def _check_output_dir_is_safe(output_dir: Path) -> None:
    """Fail fast unless ``output_dir`` is a plain, empty (or not-yet-existing) directory.

    Never silently reused, merged into, or replaced: a batch run's 25 files are published
    as one all-or-nothing unit, so anything already there -- or anything this check cannot
    positively confirm is an ordinary real directory -- is rejected before any dataset is
    touched.
    """
    if output_dir.is_symlink():
        raise BatchOutputDirectoryUnsafeError(
            output_dir=str(output_dir), reason="path is a symlink, not a plain directory"
        )
    if output_dir.exists():
        if not output_dir.is_dir():
            raise BatchOutputDirectoryUnsafeError(
                output_dir=str(output_dir), reason="path exists and is not a directory"
            )
        if any(output_dir.iterdir()):
            raise BatchOutputDirectoryUnsafeError(
                output_dir=str(output_dir),
                reason=(
                    "directory already exists and is not empty -- batch output is "
                    "published as one all-or-nothing unit and never merges into an "
                    "existing directory"
                ),
            )


def run_batch(
    *,
    manifest_path: Path,
    data_dir: Path,
    output_dir: Path,
    alpha: float = DEFAULT_ALPHA,
    fdr_level: float = DEFAULT_FDR_LEVEL,
    threshold: float = DEFAULT_THRESHOLD,
) -> BatchOutcome:
    """Process every manifest-declared dataset, in manifest order, into one published output tree.

    Builds the full 25-file tree (6 dataset directories x 4 artifacts, plus
    ``batch_summary.tsv``) in a staging directory created as a sibling of ``output_dir``
    (so the final publish step is a same-filesystem, effectively-atomic directory rename),
    and only moves it into place at ``output_dir`` after every dataset and the summary have
    succeeded. Any failure -- an invalid ``alpha``/``fdr_level``/``threshold``, a failed
    validation, a row-count mismatch, a plotting or TSV-writing error -- removes the
    staging directory and leaves ``output_dir`` untouched; no partial batch is ever
    published.
    """
    validate_threshold(threshold)
    validate_threshold(alpha)
    validate_threshold(fdr_level)

    output_dir = Path(output_dir)
    _check_output_dir_is_safe(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(manifest_path)

    staging_dir = Path(tempfile.mkdtemp(dir=output_dir.parent, prefix=".batch-staging-"))
    try:
        dataset_results = []
        for entry in manifest.datasets:
            dataset_staging_dir = staging_dir / entry.dataset_id
            dataset_results.append(
                _process_one_dataset(
                    manifest=manifest,
                    data_dir=data_dir,
                    entry=entry,
                    dataset_staging_dir=dataset_staging_dir,
                    threshold=threshold,
                    alpha=alpha,
                    fdr_level=fdr_level,
                )
            )
        results = tuple(dataset_results)

        summary_table = build_batch_summary_table(results)
        summary_path = staging_dir / "batch_summary.tsv"
        atomic_write_tsv(summary_table, summary_path)

        os.replace(staging_dir, output_dir)
    except BaseException:
        # os.replace() is inside this try too: if the final publish itself fails (e.g. a
        # race, a permission error, a cross-device rename), the staging directory must be
        # cleaned up the same as any earlier failure -- it must never be left orphaned
        # next to a still-absent output_dir.
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise

    return BatchOutcome(
        output_dir=output_dir,
        datasets=results,
        summary_path=output_dir / "batch_summary.tsv",
        total_n_tests=sum(r.n_tests for r in results),
    )


# BatchDatasetResult's fields are re-exported so tests can assert, by
# introspection, that every field is a scalar (str/int/float/None) -- never a
# DataFrame or array -- without hardcoding the field list twice.
BATCH_DATASET_RESULT_FIELD_NAMES: tuple[str, ...] = tuple(
    f.name for f in fields(BatchDatasetResult)
)
