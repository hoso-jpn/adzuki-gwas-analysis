"""High-level orchestration: validate, then load, then analyze/plot.

Every public function here follows the same shape: resolve the manifest
entry, run :func:`adzuki_gwas_analysis.validate.validate_dataset`, and raise
if it fails -- before touching pandas, matplotlib, or the output directory
at all. No function in this module produces a plot or TSV for a dataset that
failed validation.

``threshold`` throughout is the legacy ``1e-5`` visualization line used by
the original scripts. It is not a Bonferroni-corrected or genome-wide
significance threshold -- it is not derived from any multiple-testing
correction, and this module does not compute one.
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from adzuki_gwas_analysis.analysis.chromosomes import compute_manhattan_coordinates
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
from adzuki_gwas_analysis.errors import DatasetValidationFailedError, UnknownDatasetIdError
from adzuki_gwas_analysis.manifest import DatasetEntry, Manifest, load_manifest
from adzuki_gwas_analysis.validate import validate_dataset

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
    _atomic_write_tsv(table, output_path)
    return table


def _atomic_write_tsv(table: pd.DataFrame, output_path: Path) -> None:
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
    _atomic_write_tsv(table, top_variants_path)

    return AllOutcome(
        manhattan=manhattan,
        qq=qq,
        regions=region_outcomes,
        top_variants=table,
        top_variants_path=top_variants_path,
    )
