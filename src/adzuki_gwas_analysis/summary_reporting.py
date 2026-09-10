"""Plot and report explicit summary rows without reusing Dryad-specific test assumptions."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.tables import write_tsv


def lambda_diagnostic(rows: list[dict[str, Any]], meta: dict[str, Any]) -> dict[str, Any]:
    if "chi_square_df" not in meta:
        return {
            "value": None,
            "status": "not_assessed",
            "reason": "chi-square df/basis not supplied",
        }
    from scipy.stats import chi2

    values = sorted(row["neg_log10_pvalue"] for row in rows)
    middle = (
        [values[len(values) // 2]]
        if len(values) % 2
        else values[len(values) // 2 - 1 : len(values) // 2 + 1]
    )
    if max(middle) > 300:
        return {
            "value": None,
            "status": "not_assessed",
            "reason": "median exceeds supported chi-square tail",
        }
    df = meta["chi_square_df"]
    observed = sum(float(chi2.isf(10**-value, df)) for value in middle) / len(middle)
    return {
        "value": observed / float(chi2.ppf(0.5, df)),
        "df": df,
        "status": "diagnostic_only",
        "basis": meta["chi_square_basis"],
    }


def plot_summary(
    rows: list[dict[str, Any]], contig_order: list[str], stage: Path, *, alpha: float
) -> None:
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(figsize=(11, 4))
    cursor = 0
    ticks, labels = [], []
    for index, chrom in enumerate(contig_order):
        values = [row for row in rows if row["chr"] == chrom]
        if not values:
            continue
        start, end = min(row["pos"] for row in values), max(row["pos"] for row in values)
        axes.scatter(
            [cursor + row["pos"] - start for row in values],
            [row["neg_log10_pvalue"] for row in values],
            s=4,
            color=("#356a9a", "#89a9bc")[index % 2],
        )
        ticks.append(cursor + (end - start) / 2)
        labels.append(chrom)
        cursor += end - start + 1
    axes.axhline(
        -math.log10(alpha / len(rows)), color="#c45a37", linestyle="--", label="Bonferroni"
    )
    axes.set(xticks=ticks, xticklabels=labels, ylabel="-log10(p)", xlabel="Contig")
    axes.legend()
    figure.tight_layout()
    figure.savefig(stage / "manhattan.png", dpi=150)
    plt.close(figure)
    figure, axes = plt.subplots(figsize=(5, 5))
    observed = sorted((row["neg_log10_pvalue"] for row in rows), reverse=True)
    expected = [-math.log10((rank + 0.5) / len(rows)) for rank in range(len(rows))]
    axes.scatter(expected, observed, s=5)
    end = max(expected)
    axes.plot([0, end], [0, end], color="#c45a37", linestyle="--")
    axes.set(xlabel="Expected -log10(p)", ylabel="Observed -log10(p)")
    figure.tight_layout()
    figure.savefig(stage / "qq.png", dpi=150)
    plt.close(figure)


def write_summary_artifacts(
    stage: Path,
    *,
    rows: list[dict[str, Any]],
    candidates: list[dict[str, Any]],
    excluded: list[dict[str, Any]],
    meta: dict[str, Any],
    contig_order: list[str],
    alpha: float,
    fdr_level: float,
) -> dict[str, Any]:
    diagnostics = {
        "schema_version": 1,
        "n_input": len(rows) + len(excluded),
        "n_tests": len(rows),
        "n_excluded": len(excluded),
        "n_candidates": len(candidates),
        "bonferroni_discoveries": sum(row["bonferroni_significant"] for row in rows),
        "bh_discoveries": sum(row["bh_significant"] for row in rows),
        "alpha": alpha,
        "fdr_level": fdr_level,
        "lambda_gc": lambda_diagnostic(rows, meta),
        "family": {
            key: meta[key] for key in ("dataset_id", "cohort_id", "trait", "analysis_id", "test")
        },
    }
    write_tsv(stage / "normalized_summary.tsv", rows[0].keys(), rows)
    write_tsv(
        stage / "candidate_snps.tsv",
        tuple(rows[0]) + ("signal_id", "candidate_rank", "priority_tier"),
        candidates,
    )
    write_tsv(stage / "excluded_rows.tsv", ("source_row", "reason"), excluded)
    plot_summary(rows, contig_order, stage, alpha=alpha)
    unresolved = sum(row["effect_orientation"] == "unresolved_strand" for row in rows)
    (stage / "analysis_report.md").write_text(
        "# Customer summary-statistics report\n\n"
        f"Rows analyzed: {len(rows)}; rows excluded: {len(excluded)}; "
        f"candidates: {len(candidates)}.\n\n"
        f"Effect orientation unresolved: {unresolved}. Favorable alleles are not inferred.\n\n"
        "The supplied summary statistics were analyzed; individual-level GWAS was not run here. "
        "One dataset/cohort/trait/analysis/test defines the correction family. "
        "Bonferroni controls family-wise error; BH-adjusted p-values are not Storey q-values. "
        "BH guarantees require independence or PRDS, which has not been established here.\n\n"
        "Physical-distance signals are not independently defined QTLs or LD blocks. "
        "Candidate ranking is follow-up priority, not causality. Test and trait-coding metadata "
        "must be reviewed before interpreting effects. Lambda GC is diagnostic only and is "
        "not assessed without an explicit chi-square basis. "
        "No genomic-control correction is applied.\n",
        encoding="utf-8",
    )
    return diagnostics
