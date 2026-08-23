"""High-level orchestration: validate, then load, then analyze/plot.

Every public function here follows the same shape: resolve the manifest
entry, run :func:`adzuki_gwas_analysis.validate.validate_dataset`, and raise
if it fails -- before touching pandas, matplotlib, or the output directory
at all. No function in this module produces a plot or TSV for a dataset that
failed validation.

``threshold`` throughout (``run_manhattan``/``run_single_regional``/``run_regions``/
``run_all``) is the legacy ``1e-5`` visualization line used by the original scripts.
It is not a Bonferroni-corrected or genome-wide significance threshold, and drawing
it does not itself perform any multiple-testing correction.

Multiple-testing correction is a separate responsibility, orchestrated by
``run_diagnostics`` below: it validates and loads a dataset exactly like every
other ``run_*`` function here, then delegates the actual Bonferroni/
Benjamini-Hochberg/lambda_GC computation to
:mod:`adzuki_gwas_analysis.analysis.statistics` (no numeric logic lives in this
module) before writing its own, independent output files. Plotting's
``threshold`` and ``run_diagnostics``'s ``alpha``/``fdr_level`` are unrelated
values for unrelated purposes -- one draws a line on a plot, the other
computes a statistical correction -- and neither substitutes for the other.
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from adzuki_gwas_analysis.analysis.candidates import CandidatesResult, build_candidates_result
from adzuki_gwas_analysis.analysis.chromosomes import compute_manhattan_coordinates
from adzuki_gwas_analysis.analysis.diagnostics import (
    DiagnosticsResult,
    build_significant_variants_table,
    build_summary_table,
)
from adzuki_gwas_analysis.analysis.loader import load_analysis_frame
from adzuki_gwas_analysis.analysis.plotting import (
    PlotResult,
    plot_manhattan,
    plot_qq,
    plot_regional,
)
from adzuki_gwas_analysis.analysis.qq import compute_qq_points
from adzuki_gwas_analysis.analysis.region_config import RegionConfig, load_region_config
from adzuki_gwas_analysis.analysis.regions import TopVariant, select_top_variant, subset_region
from adzuki_gwas_analysis.analysis.statistics import (
    compute_bh,
    compute_bonferroni,
    compute_lambda_gc,
)
from adzuki_gwas_analysis.errors import (
    DatasetValidationFailedError,
    LoadedPvalueCountMismatchError,
    UnknownDatasetIdError,
)
from adzuki_gwas_analysis.manifest import DatasetEntry, Manifest, load_manifest
from adzuki_gwas_analysis.validate import ValidationResult, validate_dataset

DEFAULT_THRESHOLD = 1e-5


def validate_threshold(threshold: float) -> None:
    """Raise ``ValueError`` unless ``0 < threshold <= 1``."""
    if not (0.0 < threshold <= 1.0):
        raise ValueError(f"threshold must satisfy 0 < threshold <= 1, got {threshold}")


def _get_entry(manifest: Manifest, dataset_id: str) -> DatasetEntry:
    try:
        return manifest.get(dataset_id)
    except KeyError:
        raise UnknownDatasetIdError(
            dataset_id=dataset_id,
            available=tuple(e.dataset_id for e in manifest.datasets),
        ) from None


def ensure_validated(*, manifest_path: Path, data_dir: Path, dataset_id: str) -> DatasetEntry:
    """Resolve and validate ``dataset_id`` against ``data_dir``; raise if invalid.

    Returns the manifest's :class:`~adzuki_gwas_analysis.manifest.DatasetEntry`
    on success. This is the single mandatory gate every analysis operation in
    this module runs through -- nothing downstream ever sees an unvalidated
    file.
    """
    manifest = load_manifest(manifest_path)
    entry = _get_entry(manifest, dataset_id)
    result = validate_dataset(entry, data_dir)
    if not result.success:
        raise DatasetValidationFailedError(dataset_id=dataset_id, reason=result.error or "unknown")
    return entry


@dataclass(frozen=True, slots=True)
class ValidatedDatasetInfo:
    """The manifest, dataset entry, and validation result for one validated dataset.

    ``ensure_validated`` (above) discards :class:`~adzuki_gwas_analysis.validate.ValidationResult`
    after checking ``success``, keeping only the manifest's declared
    :class:`~adzuki_gwas_analysis.manifest.DatasetEntry`. Diagnostics needs more than that:
    it must confirm the number of p-values it loads equals the row count
    *validation actually counted* -- not ``manifest.toml``'s declared ``row_count`` taken on
    faith -- so this variant returns the full :class:`Manifest` (for
    ``pvalue_columns.primary``) and the :class:`~adzuki_gwas_analysis.validate.ValidationResult`
    (for its counted ``row_count`` and ``sha256``) as well.
    """

    manifest: Manifest
    entry: DatasetEntry
    validation: ValidationResult


def ensure_validated_entry(
    *, manifest: Manifest, entry: DatasetEntry, data_dir: Path
) -> ValidatedDatasetInfo:
    """Validate one already-resolved ``entry`` against an already-loaded ``manifest``.

    Unlike :func:`ensure_validated_with_result`, this never calls
    :func:`~adzuki_gwas_analysis.manifest.load_manifest` itself -- callers that already hold
    a :class:`Manifest` (e.g. :func:`adzuki_gwas_analysis.analysis.batch.run_batch`, which
    loads it exactly once for the whole batch and then validates each of its 6 entries in
    turn) use this to avoid re-parsing ``manifest.toml`` once per dataset. Calls
    :func:`~adzuki_gwas_analysis.validate.validate_dataset` exactly once, the same single
    pass over the file every other validation path here performs.
    """
    result = validate_dataset(entry, data_dir)
    if not result.success:
        raise DatasetValidationFailedError(
            dataset_id=entry.dataset_id, reason=result.error or "unknown"
        )
    return ValidatedDatasetInfo(manifest=manifest, entry=entry, validation=result)


def ensure_validated_with_result(
    *, manifest_path: Path, data_dir: Path, dataset_id: str
) -> ValidatedDatasetInfo:
    """Like :func:`ensure_validated`, but also returns the manifest and full validation result.

    Loads the manifest and resolves ``dataset_id`` to a
    :class:`~adzuki_gwas_analysis.manifest.DatasetEntry`, then delegates to
    :func:`ensure_validated_entry`. Calls
    :func:`~adzuki_gwas_analysis.validate.validate_dataset` exactly once -- the same single
    pass over the file that :func:`ensure_validated` already performs -- so using this
    instead adds no extra raw-file scan.
    """
    manifest = load_manifest(manifest_path)
    entry = _get_entry(manifest, dataset_id)
    return ensure_validated_entry(manifest=manifest, entry=entry, data_dir=data_dir)


def load_validated_frame(*, manifest_path: Path, data_dir: Path, dataset_id: str) -> pd.DataFrame:
    """Validate ``dataset_id`` then load its analysis columns into a DataFrame."""
    entry = ensure_validated(manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id)
    return load_analysis_frame(data_dir / entry.member_filename)


def run_manhattan(
    *,
    manifest_path: Path,
    data_dir: Path,
    dataset_id: str,
    output_path: Path,
    threshold: float = DEFAULT_THRESHOLD,
    title: str = "Water Permeability GWAS",
) -> PlotResult:
    """Validate, load, and render the genome-wide Manhattan plot for ``dataset_id``."""
    validate_threshold(threshold)
    df = load_validated_frame(manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id)
    coords = compute_manhattan_coordinates(df)
    return plot_manhattan(
        coords, df["pval"], threshold=threshold, output_path=output_path, title=title
    )


def run_qq(
    *,
    manifest_path: Path,
    data_dir: Path,
    dataset_id: str,
    output_path: Path,
    title: str = "QQ Plot: Water Permeability GWAS",
) -> PlotResult:
    """Validate, load, and render the QQ plot for ``dataset_id``."""
    df = load_validated_frame(manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id)
    qq_points = compute_qq_points(df["pval"].to_numpy(dtype="float64"))
    return plot_qq(qq_points, output_path=output_path, title=title)


@dataclass(frozen=True, slots=True)
class RegionalPlotOutcome:
    """A rendered regional plot plus the top variant found in that window."""

    plot: PlotResult
    top_variant: TopVariant


def run_single_regional(
    *,
    manifest_path: Path,
    data_dir: Path,
    dataset_id: str,
    chrom: str,
    start: int,
    end: int,
    output_path: Path,
    threshold: float = DEFAULT_THRESHOLD,
    title: str | None = None,
    region_id: str = "adhoc",
) -> RegionalPlotOutcome:
    """Validate, load, subset one ad-hoc region, and render its regional plot.

    Used both by the ``regional`` CLI subcommand and by the
    ``scripts/03_regional_plot.py`` backward-compatible wrapper, which takes
    ``--chrom``/``--start``/``--end`` directly rather than a region config.
    """
    validate_threshold(threshold)
    df = load_validated_frame(manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id)
    region_df = subset_region(df, chrom=chrom, start=start, end=end)
    top_variant = select_top_variant(
        region_df, dataset_id=dataset_id, region_id=region_id, chrom=chrom, start=start, end=end
    )
    plot_title = title or f"Regional Plot: {chrom}:{start}-{end}"
    plot = plot_regional(
        region_df, top_variant, threshold=threshold, output_path=output_path, title=plot_title
    )
    return RegionalPlotOutcome(plot=plot, top_variant=top_variant)


def _plot_all_regions(
    df: pd.DataFrame,
    region_config: RegionConfig,
    *,
    dataset_id: str,
    output_dir: Path,
    threshold: float,
) -> list[RegionalPlotOutcome]:
    outcomes = []
    for region in region_config.regions:
        region_df = subset_region(df, chrom=region.chrom, start=region.start, end=region.end)
        top_variant = select_top_variant(
            region_df,
            dataset_id=dataset_id,
            region_id=region.region_id,
            chrom=region.chrom,
            start=region.start,
            end=region.end,
        )
        plot = plot_regional(
            region_df,
            top_variant,
            threshold=threshold,
            output_path=output_dir / region.output_filename,
            title=region.title,
        )
        outcomes.append(RegionalPlotOutcome(plot=plot, top_variant=top_variant))
    return outcomes


def run_regions(
    *,
    manifest_path: Path,
    data_dir: Path,
    dataset_id: str,
    regions_config_path: Path,
    output_dir: Path,
    threshold: float = DEFAULT_THRESHOLD,
) -> list[RegionalPlotOutcome]:
    """Validate, load once, and render every region in ``regions_config_path``."""
    validate_threshold(threshold)
    region_config = load_region_config(regions_config_path, expected_dataset_id=dataset_id)
    df = load_validated_frame(manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id)
    return _plot_all_regions(
        df, region_config, dataset_id=dataset_id, output_dir=output_dir, threshold=threshold
    )


def top_variants_table(
    region_config: RegionConfig, df: pd.DataFrame, *, dataset_id: str
) -> pd.DataFrame:
    """Build the top-variant-per-region table, in the region config's own order.

    Column order/names (``region, chr, pos, allele1, allele0, af, beta,
    pval``) match ``results/water_permeability/top_variants_by_region.tsv``
    exactly, for byte-for-byte equivalence with that pre-existing file.
    """
    rows = []
    for region in region_config.regions:
        region_df = subset_region(df, chrom=region.chrom, start=region.start, end=region.end)
        top = select_top_variant(
            region_df,
            dataset_id=dataset_id,
            region_id=region.region_id,
            chrom=region.chrom,
            start=region.start,
            end=region.end,
        )
        rows.append(
            {
                "region": region.region_id,
                "chr": top.chrom,
                "pos": top.pos,
                "allele1": top.allele1,
                "allele0": top.allele0,
                "af": top.af,
                "beta": top.beta,
                "pval": top.pval,
            }
        )
    return pd.DataFrame(rows)


def run_top_variants(
    *,
    manifest_path: Path,
    data_dir: Path,
    dataset_id: str,
    regions_config_path: Path,
    output_path: Path,
) -> pd.DataFrame:
    """Validate, load once, and write the top-variant-per-region TSV."""
    region_config = load_region_config(regions_config_path, expected_dataset_id=dataset_id)
    df = load_validated_frame(manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id)
    table = top_variants_table(region_config, df, dataset_id=dataset_id)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_tsv(table, output_path)
    return table


def atomic_write_tsv(table: pd.DataFrame, output_path: Path) -> None:
    fd, tmp_name = tempfile.mkstemp(dir=output_path.parent, suffix=".tsv.tmp")
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        table.to_csv(tmp_path, sep="\t", index=False)
        os.replace(tmp_path, output_path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


@dataclass(frozen=True, slots=True)
class AllOutcome:
    """Every output produced by a single ``all`` run."""

    manhattan: PlotResult
    qq: PlotResult
    regions: list[RegionalPlotOutcome]
    top_variants: pd.DataFrame
    top_variants_path: Path


def run_all(
    *,
    manifest_path: Path,
    data_dir: Path,
    dataset_id: str,
    regions_config_path: Path,
    output_dir: Path,
    threshold: float = DEFAULT_THRESHOLD,
) -> AllOutcome:
    """Validate once, load once, and produce Manhattan + QQ + all regions + top-variant TSV."""
    validate_threshold(threshold)
    entry = ensure_validated(manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id)
    df = load_analysis_frame(data_dir / entry.member_filename)
    region_config = load_region_config(regions_config_path, expected_dataset_id=dataset_id)

    coords = compute_manhattan_coordinates(df)
    manhattan = plot_manhattan(
        coords,
        df["pval"],
        threshold=threshold,
        output_path=output_dir / f"{dataset_id}_manhattan.png",
    )

    qq_points = compute_qq_points(df["pval"].to_numpy(dtype="float64"))
    qq = plot_qq(qq_points, output_path=output_dir / f"{dataset_id}_qq.png")

    region_outcomes = _plot_all_regions(
        df, region_config, dataset_id=dataset_id, output_dir=output_dir, threshold=threshold
    )

    table = top_variants_table(region_config, df, dataset_id=dataset_id)
    top_variants_path = output_dir / "top_variants_by_region.tsv"
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_tsv(table, top_variants_path)

    return AllOutcome(
        manhattan=manhattan,
        qq=qq,
        regions=region_outcomes,
        top_variants=table,
        top_variants_path=top_variants_path,
    )


DEFAULT_ALPHA = 0.05
DEFAULT_FDR_LEVEL = 0.05
DEFAULT_LAMBDA_GC_DF = 1


@dataclass(frozen=True, slots=True)
class DiagnosticsOutcome:
    """Every output produced by a single ``diagnostics`` run."""

    result: DiagnosticsResult
    summary_path: Path
    significant_variants_path: Path


def extract_primary_pvalues(
    info: ValidatedDatasetInfo, variant_df: pd.DataFrame
) -> tuple[str, pd.Series[float]]:
    """Return ``(pvalue_column, pvalues)`` for an already-validated, already-loaded dataset.

    ``pvalue_column`` comes from ``info.manifest.pvalue_columns.primary`` -- never a
    hardcoded ``"pval"`` literal -- and the returned p-value count is checked against
    ``info.validation.row_count`` (the row count schema v1 validation actually counted,
    not ``manifest.toml``'s declared ``row_count`` taken on faith), raising
    :class:`~adzuki_gwas_analysis.errors.LoadedPvalueCountMismatchError` on a mismatch.
    Shared by :func:`compute_diagnostics_result` and
    :mod:`adzuki_gwas_analysis.analysis.batch`, which both start from a dataset that has
    already been validated and loaded exactly once.
    """
    pvalue_column = info.manifest.pvalue_columns.primary
    pvalues = variant_df[pvalue_column]
    if len(pvalues) != info.validation.row_count:
        raise LoadedPvalueCountMismatchError(
            dataset_id=info.entry.dataset_id,
            validated_row_count=info.validation.row_count,
            loaded_count=len(pvalues),
        )
    return pvalue_column, pvalues


def compute_diagnostics_result(
    info: ValidatedDatasetInfo,
    variant_df: pd.DataFrame,
    *,
    alpha: float = DEFAULT_ALPHA,
    fdr_level: float = DEFAULT_FDR_LEVEL,
) -> DiagnosticsResult:
    """Compute Bonferroni/BH/lambda_GC for an already-validated, already-loaded dataset.

    Pure computation, no file I/O: takes the :class:`ValidatedDatasetInfo` and analysis
    ``variant_df`` a caller has already produced via one
    :func:`ensure_validated_with_result` and one
    :func:`~adzuki_gwas_analysis.analysis.loader.load_analysis_frame` call, and returns the
    same :class:`~adzuki_gwas_analysis.analysis.diagnostics.DiagnosticsResult`
    that :func:`run_diagnostics` writes to disk. This is the internal API
    :mod:`adzuki_gwas_analysis.analysis.batch` uses so that a batch run's Manhattan plot, QQ
    plot, and diagnostics for one dataset all share the exact same validation pass and
    loaded DataFrame, rather than each re-validating and re-loading the same file.

    The multiple-testing family is fixed: this one dataset's manifest-declared primary
    p-value column, over every variant validation counted for this one file -- never the
    other 5 Dryad files, never Miyagi+Shumari, never the 3 traits, never a post-hoc region,
    and never ``p_wald``/``p_score``.
    """
    pvalue_column, pvalues_series = extract_primary_pvalues(info, variant_df)
    pvalues = pvalues_series.to_numpy(dtype="float64")

    bonferroni = compute_bonferroni(pvalues, alpha=alpha)
    bh = compute_bh(pvalues, fdr_level=fdr_level)
    lambda_gc = compute_lambda_gc(pvalues, df=DEFAULT_LAMBDA_GC_DF)

    return DiagnosticsResult(
        dataset_id=info.entry.dataset_id,
        source_sha256=info.validation.sha256,
        pvalue_column=pvalue_column,
        n_tests=len(pvalues),
        bonferroni=bonferroni,
        bh=bh,
        lambda_gc=lambda_gc,
    )


def run_diagnostics(
    *,
    manifest_path: Path,
    data_dir: Path,
    dataset_id: str,
    output_dir: Path,
    alpha: float = DEFAULT_ALPHA,
    fdr_level: float = DEFAULT_FDR_LEVEL,
) -> DiagnosticsOutcome:
    """Validate, load once, compute Bonferroni/BH/lambda_GC, and write both TSVs atomically.

    ``alpha``/``fdr_level`` validation and the family's numeric computation (delegated to
    :func:`compute_diagnostics_result`) all happen in memory before either output file is
    created, so a contract, count-mismatch, or numeric-input error leaves no output behind
    -- matching every other ``run_*`` function in this module.
    """
    info = ensure_validated_with_result(
        manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id
    )
    variant_df = load_analysis_frame(data_dir / info.entry.member_filename)
    result = compute_diagnostics_result(info, variant_df, alpha=alpha, fdr_level=fdr_level)

    summary_table = build_summary_table(result)
    significant_table = build_significant_variants_table(variant_df, result)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "statistical_diagnostics.tsv"
    significant_variants_path = output_dir / "significant_variants.tsv"
    atomic_write_tsv(summary_table, summary_path)
    atomic_write_tsv(significant_table, significant_variants_path)

    return DiagnosticsOutcome(
        result=result,
        summary_path=summary_path,
        significant_variants_path=significant_variants_path,
    )


@dataclass(frozen=True, slots=True)
class CandidatesOutcome:
    """Every output produced by a single ``candidates`` run."""

    diagnostics_result: DiagnosticsResult
    candidates_result: CandidatesResult
    summary_path: Path
    significant_variants_path: Path
    association_peaks_path: Path
    candidate_snps_path: Path
    candidate_ranking_path: Path


def run_candidates(
    *,
    manifest_path: Path,
    data_dir: Path,
    dataset_id: str,
    output_dir: Path,
    clustering_distance: int,
    alpha: float = DEFAULT_ALPHA,
    fdr_level: float = DEFAULT_FDR_LEVEL,
) -> CandidatesOutcome:
    """Validate, load once, compute diagnostics, cluster, and write all 5 output TSVs.

    Self-sufficient like :func:`run_diagnostics`: this does not read a previous
    ``diagnostics`` run's output back from disk. It validates ``dataset_id`` once, loads its
    analysis DataFrame once, computes the same Bonferroni/BH/lambda_GC diagnostics and the
    same significant-variant population :func:`run_diagnostics` would, and then clusters that
    in-memory population (via
    :func:`~adzuki_gwas_analysis.analysis.candidates.build_candidates_result`) into
    ``association_peaks.tsv``/``candidate_snps.tsv``/``candidate_ranking.tsv`` -- so a
    ``candidates`` run's ``statistical_diagnostics.tsv``/``significant_variants.tsv`` are
    always consistent with the ``alpha``/``fdr_level`` that same invocation used, never a
    stale file from an earlier, possibly differently-parameterized ``diagnostics`` run.

    ``clustering_distance`` has no default anywhere in this repository (see
    :mod:`adzuki_gwas_analysis.analysis.candidates`'s module docstring) and must always be
    supplied explicitly.
    """
    info = ensure_validated_with_result(
        manifest_path=manifest_path, data_dir=data_dir, dataset_id=dataset_id
    )
    variant_df = load_analysis_frame(data_dir / info.entry.member_filename)
    diagnostics_result = compute_diagnostics_result(
        info, variant_df, alpha=alpha, fdr_level=fdr_level
    )

    summary_table = build_summary_table(diagnostics_result)
    significant_table = build_significant_variants_table(variant_df, diagnostics_result)
    candidates_result = build_candidates_result(
        significant_table,
        dataset_id=info.entry.dataset_id,
        reference=info.entry.reference,
        trait=info.entry.trait,
        clustering_distance=clustering_distance,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "statistical_diagnostics.tsv"
    significant_variants_path = output_dir / "significant_variants.tsv"
    association_peaks_path = output_dir / "association_peaks.tsv"
    candidate_snps_path = output_dir / "candidate_snps.tsv"
    candidate_ranking_path = output_dir / "candidate_ranking.tsv"
    atomic_write_tsv(summary_table, summary_path)
    atomic_write_tsv(significant_table, significant_variants_path)
    atomic_write_tsv(candidates_result.association_peaks, association_peaks_path)
    atomic_write_tsv(candidates_result.candidate_snps, candidate_snps_path)
    atomic_write_tsv(candidates_result.candidate_ranking, candidate_ranking_path)

    return CandidatesOutcome(
        diagnostics_result=diagnostics_result,
        candidates_result=candidates_result,
        summary_path=summary_path,
        significant_variants_path=significant_variants_path,
        association_peaks_path=association_peaks_path,
        candidate_snps_path=candidate_snps_path,
        candidate_ranking_path=candidate_ranking_path,
    )
