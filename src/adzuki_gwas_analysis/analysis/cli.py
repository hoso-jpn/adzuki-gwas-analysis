"""Command-line entry point for validated GWAS re-analysis and plotting.

Usage::

    adzuki-gwas-analyze manhattan --output-dir out/
    adzuki-gwas-analyze qq --output-dir out/
    adzuki-gwas-analyze regional --chrom Chr07 --start 5000000 --end 7000000 \\
        --output out/regional.png
    adzuki-gwas-analyze regions --output-dir out/
    adzuki-gwas-analyze top-variants --output out/top_variants_by_region.tsv
    adzuki-gwas-analyze all --output-dir out/
    adzuki-gwas-analyze diagnostics --output-dir out/
    adzuki-gwas-analyze batch --output-dir out/
    adzuki-gwas-analyze candidates --output-dir out/ --clustering-distance 50000
    adzuki-gwas-analyze batch --output-dir out/ --clustering-distance 50000
    adzuki-gwas-analyze report --analysis-dir batch-out/ --output-dir delivery/

Every subcommand validates its dataset against the schema v1 manifest
(:mod:`adzuki_gwas_analysis.validate`) before producing any output, and exits
non-zero with no output written if that validation fails. ``pval`` is always
the likelihood-ratio-test (LRT) p-value (the manifest's declared
``pvalue_columns.primary``); the ``--threshold`` line drawn on Manhattan/
regional plots is a legacy visualization threshold (default ``1e-5``), not a
Bonferroni-corrected or genome-wide significance level.

``diagnostics`` is a separate multiple-testing-correction subcommand (Bonferroni
FWER, Benjamini-Hochberg FDR, and the genomic inflation factor lambda_GC), computed
over one dataset's manifest-declared primary p-value column only -- see
:mod:`adzuki_gwas_analysis.analysis.diagnostics`. It takes its own ``--alpha``/
``--fdr-level`` flags, not ``--threshold`` (a distinct, legacy visualization
concept), and does not change what ``manhattan``/``qq``/``regional``/``regions``/
``top-variants``/``all`` produce.

``batch`` (see :mod:`adzuki_gwas_analysis.analysis.batch`) processes every dataset
declared in the manifest, in the manifest's own order -- there is no ``--dataset-id``.
For each dataset it validates once, loads its analysis DataFrame once, and reuses that
one DataFrame for its Manhattan plot, QQ plot, and diagnostics before moving to the
next dataset; no two datasets' DataFrames are ever held at once, and each dataset's
Bonferroni/BH/lambda_GC family stays independent of the other 5. Output is a
6-directory, 25-file tree (4 files per dataset plus a root ``batch_summary.tsv``),
published to ``--output-dir`` only after every dataset succeeds -- ``--output-dir``
must not already exist as a non-empty directory. ``batch`` does not render regional
plots or a top-variant-by-region TSV (both specific to ``miyagi_water_permeability``'s
post-hoc regions) and does not change ``all``'s own output.

``candidates`` (see :mod:`adzuki_gwas_analysis.analysis.candidates`) clusters one dataset's
own Bonferroni-or-BH-significant variants (the same population ``diagnostics`` writes to
``significant_variants.tsv``) into physical-distance "signals" and writes
``association_peaks.tsv``/``candidate_snps.tsv``/``candidate_ranking.tsv`` alongside its own
``statistical_diagnostics.tsv``/``significant_variants.tsv`` (5 files total; it is
self-sufficient and does not read a previous ``diagnostics`` run's output back from disk).
``--clustering-distance`` (base pairs) is **required and has no default anywhere in this
package** -- there is no LD or independently-defined-QTL estimate in this repository from
which a default window could be justified, so it must always be chosen and stated
explicitly. A physical-distance cluster is not an LD block or a QTL interval, a cluster's
lead variant is not asserted to be causal, and ``priority_tier``/``priority_reasons``
describe downstream validation priority only -- never a biological-importance ranking or a
validated breeding marker. ``batch --clustering-distance N`` additively produces the same 3
files inside every one of the 6 dataset directories (43 files total) when supplied; omitting
it (the default) leaves ``batch``'s original 25-file, 4-files-per-dataset output completely
unchanged. Candidates are always clustered and ranked per dataset -- the 6 datasets' (or,
for ``candidates``, one dataset's) significant variants are never pooled into one shared
ranking population, and Miyagi/Shumari coordinates are never combined.

``report`` (see :mod:`adzuki_gwas_analysis.analysis.report`) is a **consumer** of an
existing ``batch`` output, not a new analysis: it never re-validates ``.assoc.txt`` files,
never recomputes Bonferroni/BH/lambda_GC, and never re-clusters candidates. ``--analysis-dir``
must be a **candidate-enabled** ``batch`` output (``batch_summary.tsv`` with
``schema_version=2``, i.e. that ``batch`` run was given an explicit
``--clustering-distance``) -- a schema-v1 (no-candidates) ``batch`` output is rejected with
an actionable error rather than silently re-clustering with an implicit distance. Publishes
a self-contained delivery package to ``--output-dir``: ``executive_summary.md``,
``analysis_report.md``, an ``artifacts/`` copy of every dataset's derived files (never the
raw ``.assoc.txt`` files), and a ``reproducibility/`` folder
(``input_checksums.tsv``/``software_versions.json``/``run_manifest.json``). Every
cross-artifact count in the report (n_tests, discoveries, lambda_GC, signal/candidate
counts) is checked against its source file before anything is written; a mismatch aborts
with no output produced. No network access or external service is used anywhere in
``report``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from adzuki_gwas_analysis.analysis.batch import run_batch
from adzuki_gwas_analysis.analysis.pipeline import (
    DEFAULT_ALPHA,
    DEFAULT_FDR_LEVEL,
    DEFAULT_THRESHOLD,
    run_all,
    run_candidates,
    run_diagnostics,
    run_manhattan,
    run_qq,
    run_regions,
    run_single_regional,
    run_top_variants,
)
from adzuki_gwas_analysis.analysis.report import run_report
from adzuki_gwas_analysis.errors import GwasContractError
from adzuki_gwas_analysis.provenance import run_audited_analysis

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

    diagnostics = subparsers.add_parser(
        "diagnostics", help="Bonferroni / BH-FDR / genomic inflation factor (lambda_GC)"
    )
    diagnostics.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help="Path to manifest.toml"
    )
    diagnostics.add_argument(
        "--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Directory with .assoc.txt files"
    )
    diagnostics.add_argument(
        "--dataset-id", default=DEFAULT_DATASET_ID, help="Manifest dataset_id to analyze"
    )
    diagnostics.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help=(
            "Output directory (required -- unlike the other subcommands, there is no "
            "default, so a bare `diagnostics` invocation can never write into the "
            f"tracked {DEFAULT_OUTPUT_DIR}/ during real-data smoke testing)"
        ),
    )
    diagnostics.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help="Bonferroni family-wise alpha (not the legacy --threshold)",
    )
    diagnostics.add_argument(
        "--fdr-level",
        type=float,
        default=DEFAULT_FDR_LEVEL,
        help="Benjamini-Hochberg FDR level (not the legacy --threshold)",
    )

    candidates = subparsers.add_parser(
        "candidates",
        help=(
            "Cluster one dataset's significant variants into physical-distance signals "
            "and rank them for downstream validation priority"
        ),
    )
    candidates.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help="Path to manifest.toml"
    )
    candidates.add_argument(
        "--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Directory with .assoc.txt files"
    )
    candidates.add_argument(
        "--dataset-id", default=DEFAULT_DATASET_ID, help="Manifest dataset_id to analyze"
    )
    candidates.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory (required, same convention as diagnostics/batch)",
    )
    candidates.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help="Bonferroni family-wise alpha (not the legacy --threshold)",
    )
    candidates.add_argument(
        "--fdr-level",
        type=float,
        default=DEFAULT_FDR_LEVEL,
        help="Benjamini-Hochberg FDR level (not the legacy --threshold)",
    )
    candidates.add_argument(
        "--clustering-distance",
        type=int,
        required=True,
        help=(
            "Physical-distance clustering window in base pairs (required -- no default "
            "exists anywhere in this package; not an LD or QTL-interval estimate)"
        ),
    )

    batch = subparsers.add_parser(
        "batch",
        help=(
            "Sequentially process all 6 manifest-declared datasets: Manhattan + QQ + "
            "diagnostics for each, plus one batch_summary.tsv"
        ),
    )
    batch.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help="Path to manifest.toml"
    )
    batch.add_argument(
        "--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Directory with .assoc.txt files"
    )
    batch.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help=(
            "Output directory (required, must not already exist non-empty). No "
            "--dataset-id: batch always processes every dataset in manifest.toml, in "
            "its own declared order."
        ),
    )
    batch.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help="Bonferroni family-wise alpha, applied independently per dataset",
    )
    batch.add_argument(
        "--fdr-level",
        type=float,
        default=DEFAULT_FDR_LEVEL,
        help="Benjamini-Hochberg FDR level, applied independently per dataset",
    )
    batch.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help=(
            "Legacy visualization threshold line for each dataset's Manhattan/QQ plots "
            "(not a Bonferroni-corrected or genome-wide-significance level)"
        ),
    )
    batch.add_argument(
        "--clustering-distance",
        type=int,
        default=None,
        help=(
            "Optional physical-distance clustering window in base pairs. Omitted by "
            "default -- batch's original 25-file output is unchanged. When supplied, "
            "adds association_peaks.tsv/candidate_snps.tsv/candidate_ranking.tsv to every "
            "dataset directory, clustered and ranked independently per dataset. No "
            "default value exists anywhere in this package; not an LD or QTL-interval "
            "estimate."
        ),
    )

    for producer in (batch, candidates):
        producer.add_argument(
            "--record-provenance",
            action="store_true",
            help="Atomically record generation environment and artifact checksums",
        )

    report = subparsers.add_parser(
        "report",
        help=(
            "Generate a customer-facing report and audit package from an existing "
            "candidate-enabled batch output (does not re-run any analysis)"
        ),
    )
    report.add_argument(
        "--analysis-dir",
        type=Path,
        required=True,
        help=(
            "A candidate-enabled `batch` output directory (batch_summary.tsv with "
            "schema_version=2, i.e. that batch run used an explicit "
            "--clustering-distance). A schema-v1 batch output is rejected."
        ),
    )
    report.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Delivery package output directory (required, must not already exist non-empty)",
    )

    report.add_argument(
        "--require-provenance",
        action="store_true",
        help="Reject legacy or incomplete generation provenance",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the analysis CLI; returns the process exit code (0 = success)."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command in ("batch", "candidates") and args.record_provenance:
            run_audited_analysis(
                command=args.command,
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                output_dir=args.output_dir,
                dataset_id=getattr(args, "dataset_id", None),
                alpha=args.alpha,
                fdr_level=args.fdr_level,
                clustering_distance=args.clustering_distance,
                threshold=getattr(args, "threshold", DEFAULT_THRESHOLD),
            )
            print(f"Saved audited analysis: {args.output_dir}")
            return 0
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

        elif args.command == "diagnostics":
            diagnostics_outcome = run_diagnostics(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                dataset_id=args.dataset_id,
                output_dir=args.output_dir,
                alpha=args.alpha,
                fdr_level=args.fdr_level,
            )
            diagnostics_result = diagnostics_outcome.result
            print(
                f"dataset_id={diagnostics_result.dataset_id} "
                f"pvalue_column={diagnostics_result.pvalue_column} "
                f"n_tests={diagnostics_result.n_tests}"
            )
            print(
                f"Bonferroni: threshold={diagnostics_result.bonferroni.threshold:.6g} "
                f"discoveries={diagnostics_result.bonferroni.discoveries}"
            )
            print(
                f"BH: fdr_level={diagnostics_result.bh.fdr_level} "
                f"discoveries={diagnostics_result.bh.discoveries} "
                f"raw_p_cutoff={diagnostics_result.bh.raw_p_cutoff}"
            )
            print(
                f"lambda_GC: df={diagnostics_result.lambda_gc.df} "
                f"value={diagnostics_result.lambda_gc.lambda_gc:.6g} "
                f"(diagnostic only; not a genome-wide-significance verdict)"
            )
            print(f"Saved: {diagnostics_outcome.summary_path}")
            print(f"Saved: {diagnostics_outcome.significant_variants_path}")

        elif args.command == "candidates":
            candidates_outcome = run_candidates(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                dataset_id=args.dataset_id,
                output_dir=args.output_dir,
                clustering_distance=args.clustering_distance,
                alpha=args.alpha,
                fdr_level=args.fdr_level,
            )
            candidates_result = candidates_outcome.candidates_result
            print(
                f"dataset_id={args.dataset_id} "
                f"clustering_distance={args.clustering_distance}bp "
                f"n_signals={candidates_result.n_signals} "
                f"n_candidates={candidates_result.n_candidates} "
                f"(physical-distance clusters, not LD blocks or QTL intervals)"
            )
            print(f"Saved: {candidates_outcome.summary_path}")
            print(f"Saved: {candidates_outcome.significant_variants_path}")
            print(f"Saved: {candidates_outcome.association_peaks_path}")
            print(f"Saved: {candidates_outcome.candidate_snps_path}")
            print(f"Saved: {candidates_outcome.candidate_ranking_path}")

        elif args.command == "batch":
            batch_outcome = run_batch(
                manifest_path=args.manifest,
                data_dir=args.data_dir,
                output_dir=args.output_dir,
                alpha=args.alpha,
                fdr_level=args.fdr_level,
                threshold=args.threshold,
                clustering_distance=args.clustering_distance,
            )
            for dataset_result in batch_outcome.datasets:
                line = (
                    f"{dataset_result.dataset_id}: n_tests={dataset_result.n_tests} "
                    f"bonferroni_discoveries={dataset_result.bonferroni_discoveries} "
                    f"bh_discoveries={dataset_result.bh_discoveries} "
                    f"lambda_gc={dataset_result.lambda_gc:.6g}"
                )
                if dataset_result.n_signals is not None:
                    line += (
                        f" n_signals={dataset_result.n_signals} "
                        f"n_candidates={dataset_result.n_candidates}"
                    )
                print(line)
            print(
                f"Total n_tests across {len(batch_outcome.datasets)} datasets: "
                f"{batch_outcome.total_n_tests} (not a shared multiple-testing family)"
            )
            print(f"Saved: {batch_outcome.summary_path}")
            n_datasets = len(batch_outcome.datasets)
            print(f"Saved: {batch_outcome.output_dir} ({n_datasets} dataset directories)")

        elif args.command == "report":
            report_outcome = run_report(
                analysis_dir=args.analysis_dir,
                output_dir=args.output_dir,
                require_provenance=args.require_provenance,
            )
            print(
                f"datasets={report_outcome.n_datasets} "
                f"n_signals={report_outcome.n_signals} "
                f"n_candidates={report_outcome.n_candidates} "
                f"(simple sums across independent per-dataset counts, not a shared "
                f"multiple-testing family or a common biological locus count)"
            )
            print(f"Saved: {report_outcome.executive_summary_path}")
            print(f"Saved: {report_outcome.analysis_report_path}")
            print(f"Saved: {report_outcome.run_manifest_path}")
            print(f"Saved: {report_outcome.output_dir}")

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
