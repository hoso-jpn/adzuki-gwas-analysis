"""Command-line entry point for validated GWAS re-analysis and plotting.

Usage::

    adzuki-gwas-analyze manhattan --output-dir out/
    adzuki-gwas-analyze qq --output-dir out/
    adzuki-gwas-analyze regional --chrom Chr07 --start 5000000 --end 7000000 \\
        --output out/regional.png
    adzuki-gwas-analyze regions --output-dir out/
    adzuki-gwas-analyze top-variants --output out/top_variants_by_region.tsv
    adzuki-gwas-analyze all --output-dir out/

Every subcommand validates its dataset against the schema v1 manifest
(:mod:`adzuki_gwas_analysis.validate`) before producing any output, and exits
non-zero with no output written if that validation fails. ``pval`` is always
the likelihood-ratio-test (LRT) p-value (the manifest's declared
``pvalue_columns.primary``); the ``--threshold`` line drawn on Manhattan/
regional plots is a legacy visualization threshold (default ``1e-5``), not a
Bonferroni-corrected or genome-wide significance level.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from adzuki_gwas_analysis.analysis.pipeline import (
    DEFAULT_THRESHOLD,
    run_all,
    run_manhattan,
    run_qq,
    run_regions,
    run_single_regional,
    run_top_variants,
)
from adzuki_gwas_analysis.errors import GwasContractError

DEFAULT_DATASET_ID = "miyagi_water_permeability"
DEFAULT_MANIFEST = Path("manifest.toml")
DEFAULT_DATA_DIR = Path("data/raw")
DEFAULT_REGIONS_CONFIG = Path("config/water_permeability_regions.toml")
DEFAULT_OUTPUT_DIR = Path("plots")


def _add_common_arguments(parser: argparse.ArgumentParser, *, needs_output_dir: bool) -> None:
    parser.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help="Path to manifest.toml"
    )
    parser.add_argument(
        "--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Directory with .assoc.txt files"
    )
    parser.add_argument(
        "--dataset-id", default=DEFAULT_DATASET_ID, help="Manifest dataset_id to analyze"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help="Legacy visualization threshold line (not a corrected significance level)",
    )
    if needs_output_dir:
        parser.add_argument(
            "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Output directory"
        )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="adzuki-gwas-analyze",
        description=(
            "Validated re-analysis and visualization of GWAS summary statistics. "
            "pval is the likelihood-ratio-test p-value; --threshold is a legacy "
            "visualization line, not a statistically corrected significance level."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    manhattan = subparsers.add_parser("manhattan", help="Genome-wide Manhattan plot")
    _add_common_arguments(manhattan, needs_output_dir=True)

    qq = subparsers.add_parser("qq", help="QQ plot")
    _add_common_arguments(qq, needs_output_dir=True)

    regional = subparsers.add_parser("regional", help="Single ad-hoc regional plot")
    _add_common_arguments(regional, needs_output_dir=False)
    regional.add_argument("--chrom", required=True, help="Chromosome name, e.g. Chr07")
    regional.add_argument("--start", type=int, required=True, help="Start position")
    regional.add_argument("--end", type=int, required=True, help="End position")
    regional.add_argument("--output", type=Path, required=True, help="Output PNG file")
    regional.add_argument("--title", default=None, help="Plot title")

    regions = subparsers.add_parser("regions", help="All regional plots from --regions-config")
    _add_common_arguments(regions, needs_output_dir=True)
    regions.add_argument("--regions-config", type=Path, default=DEFAULT_REGIONS_CONFIG)

    top_variants = subparsers.add_parser("top-variants", help="Top-variant-per-region TSV")
    _add_common_arguments(top_variants, needs_output_dir=False)
    top_variants.add_argument("--regions-config", type=Path, default=DEFAULT_REGIONS_CONFIG)
    top_variants.add_argument("--output", type=Path, required=True, help="Output TSV file")

    all_cmd = subparsers.add_parser("all", help="Manhattan + QQ + regions + top-variants")
    _add_common_arguments(all_cmd, needs_output_dir=True)
    all_cmd.add_argument("--regions-config", type=Path, default=DEFAULT_REGIONS_CONFIG)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the analysis CLI; returns the process exit code (0 = success)."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "manhattan":
            result = run_manhattan(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                dataset_id=args.dataset_id,
                output_path=args.output_dir / f"{args.dataset_id}_manhattan.png",
                threshold=args.threshold,
            )
            print(f"Saved: {result.output_path} ({result.width_px}x{result.height_px}px)")

        elif args.command == "qq":
            result = run_qq(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                dataset_id=args.dataset_id,
                output_path=args.output_dir / f"{args.dataset_id}_qq.png",
            )
            print(f"Saved: {result.output_path} ({result.width_px}x{result.height_px}px)")

        elif args.command == "regional":
            outcome = run_single_regional(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                dataset_id=args.dataset_id,
                chrom=args.chrom,
                start=args.start,
                end=args.end,
                output_path=args.output,
                threshold=args.threshold,
                title=args.title,
            )
            print(f"Saved: {outcome.plot.output_path}")
            print(f"Top variant in region: {outcome.top_variant}")

        elif args.command == "regions":
            outcomes = run_regions(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                dataset_id=args.dataset_id,
                regions_config_path=args.regions_config,
                output_dir=args.output_dir,
                threshold=args.threshold,
            )
            for outcome in outcomes:
                print(f"Saved: {outcome.plot.output_path}")

        elif args.command == "top-variants":
            table = run_top_variants(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                dataset_id=args.dataset_id,
                regions_config_path=args.regions_config,
                output_path=args.output,
            )
            print(table.to_string(index=False))
            print(f"Saved: {args.output}")

        elif args.command == "all":
            all_outcome = run_all(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                dataset_id=args.dataset_id,
                regions_config_path=args.regions_config,
                output_dir=args.output_dir,
                threshold=args.threshold,
            )
            print(f"Saved: {all_outcome.manhattan.output_path}")
            print(f"Saved: {all_outcome.qq.output_path}")
            for region_outcome in all_outcome.regions:
                print(f"Saved: {region_outcome.plot.output_path}")
            print(f"Saved: {all_outcome.top_variants_path}")

        else:  # pragma: no cover - argparse enforces valid choices
            parser.error(f"unknown command {args.command!r}")
            return 2

    except GwasContractError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
