"""Pure, deterministic content builders for the customer report delivery package.

Every function here takes an already-validated
:class:`~adzuki_gwas_analysis.analysis.report_validation.ValidatedAnalysisDir` (or plain
scalars) and returns a string or a JSON-serializable ``dict`` -- no file I/O, no network
access, no LLM call, and no wall-clock timestamp (so the same input always produces
byte-identical output; see :mod:`adzuki_gwas_analysis.analysis.report`'s module docstring
for why determinism is preferred over embedding a generation timestamp).

Every reference to a figure or TSV in the generated text is a path relative to the
delivery package root (``artifacts/<dataset_id>/...``) -- never an absolute path, and
never anything from the local machine (username, hostname, ``$HOME``, a temp directory).
"""

from __future__ import annotations

from adzuki_gwas_analysis.analysis.batch import TRAIT_LABELS
from adzuki_gwas_analysis.analysis.report_validation import ValidatedAnalysisDir, ValidatedDataset

#: Bumped only if run_manifest.json's field set/meaning changes.
RUN_MANIFEST_SCHEMA_VERSION = 1

#: Bumped only if software_versions.json's field set/meaning changes.
SOFTWARE_VERSIONS_SCHEMA_VERSION = 1

#: Sentinel used whenever this report cannot observe a fact about the *analysis*-generation
#: environment (as opposed to *this report's own* generation environment) -- see
#: :mod:`adzuki_gwas_analysis.analysis.report`'s provenance-collection code. Never replaced
#: with a guess.
UNAVAILABLE_FROM_SOURCE_ARTIFACTS = "unavailable_from_source_artifacts"

#: How many candidates to preview inline in executive_summary.md per dataset before
#: pointing the reader to candidate_ranking.tsv for the rest.
_EXEC_SUMMARY_PREVIEW_ROWS = 3

#: How many candidates to preview inline in analysis_report.md per dataset.
_ANALYSIS_REPORT_PREVIEW_ROWS = 10

SCIENTIFIC_SCOPE: dict[str, object] = {
    "post_hoc_analysis_of_published_summary_statistics": True,
    "gwas_rerun": False,
    "association_signal_definition": "physical_distance_cluster_only",
    "association_signal_is_ld_block": False,
    "association_signal_is_independent_qtl": False,
    "lead_variant_is_causal_variant_claim": False,
    "candidate_is_causal_variant_claim": False,
    "candidate_is_validated_breeding_marker_claim": False,
    "priority_tier_is_biological_importance_ranking": False,
    "cross_reference_coordinate_integration": False,
    "cross_dataset_shared_candidate_ranking": False,
}


def _trait_label(trait: str) -> str:
    return TRAIT_LABELS.get(trait, trait)


def _preview_candidates(dataset: ValidatedDataset, *, limit: int) -> list[dict[str, object]]:
    """Return up to ``limit`` top-ranked candidates, joined with pval/beta for display.

    ``candidate_ranking.tsv`` alone has no ``pval``/``beta`` column (see
    :mod:`adzuki_gwas_analysis.analysis.candidates`'s ``CANDIDATE_RANKING_COLUMNS``); this
    joins it against the same dataset's already-loaded ``candidate_snps.tsv`` on the
    variant-identity key schema v1 uses everywhere (``chr, pos, allele1, allele0``) purely
    for display -- it recomputes nothing.
    """
    if dataset.n_candidates == 0:
        return []
    ranking = dataset.candidate_ranking_df.sort_values("candidate_rank").head(limit)
    join_keys = ["chr", "pos", "allele1", "allele0"]
    merged = ranking.merge(
        dataset.candidate_snps_df[[*join_keys, "pval", "beta"]],
        on=join_keys,
        how="left",
        validate="one_to_one",
    )
    return [
        {
            "candidate_rank": int(row["candidate_rank"]),
            "priority_tier": int(row["priority_tier"]),
            "signal_id": str(row["signal_id"]),
            "chr": str(row["chr"]),
            "pos": int(row["pos"]),
            "pval": float(row["pval"]),
            "beta": float(row["beta"]),
            "priority_reasons": str(row["priority_reasons"]),
        }
        for _, row in merged.iterrows()
    ]


def _dataset_artifact_links(dataset: ValidatedDataset) -> dict[str, str]:
    prefix = f"artifacts/{dataset.dataset_id}"
    return {
        "manhattan": f"{prefix}/{dataset.relative_paths['manhattan_path'].split('/')[-1]}",
        "qq": f"{prefix}/{dataset.relative_paths['qq_path'].split('/')[-1]}",
        "diagnostics": f"{prefix}/statistical_diagnostics.tsv",
        "significant_variants": f"{prefix}/significant_variants.tsv",
        "association_peaks": f"{prefix}/association_peaks.tsv",
        "candidate_snps": f"{prefix}/candidate_snps.tsv",
        "candidate_ranking": f"{prefix}/candidate_ranking.tsv",
    }


def build_executive_summary_markdown(validated: ValidatedAnalysisDir) -> str:
    """Build ``executive_summary.md``: a short, non-specialist-readable overview."""
    references = sorted({d.reference for d in validated.datasets})
    traits = sorted({d.trait for d in validated.datasets})
    total_candidates = sum(d.n_candidates for d in validated.datasets)
    total_signals = sum(d.n_signals for d in validated.datasets)

    lines = [
        "# Executive Summary",
        "",
        "## Analysis scope",
        "",
        "- Post-hoc re-analysis of already-published GWAS summary statistics. "
        "**The underlying GWAS was not re-run.**",
        f"- {len(validated.datasets)} dataset(s) across reference genome(s) "
        f"{', '.join(references)} and trait(s) "
        f"{', '.join(_trait_label(t) for t in traits)}.",
        f"- Association signals were grouped only by physical distance "
        f"(`clustering_distance={validated.clustering_distance}` bp for this run) -- "
        "**not** linkage-disequilibrium (LD) blocks and **not** independently "
        "established QTL intervals.",
        f"- Bonferroni family-wise alpha = {validated.alpha}; Benjamini-Hochberg FDR "
        f"level = {validated.fdr_level}, each computed **independently per dataset** "
        "(never pooling datasets, references, or traits into one correction family).",
        "",
        "## Dataset overview",
        "",
        "| Dataset | Reference | Trait | n_tests | Bonferroni | BH | lambda_GC | Signals "
        "| Candidates |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for dataset in validated.datasets:
        lines.append(
            f"| `{dataset.dataset_id}` | {dataset.reference} | {_trait_label(dataset.trait)} "
            f"| {dataset.n_tests:,} | {dataset.bonferroni_discoveries:,} "
            f"| {dataset.bh_discoveries:,} | {dataset.lambda_gc:.6g} "
            f"| {dataset.n_signals:,} | {dataset.n_candidates:,} |"
        )
    lines += [
        "",
        f"**Totals across all {len(validated.datasets)} dataset(s):** "
        f"{total_signals:,} signal(s), {total_candidates:,} candidate(s). This is a simple "
        "sum of independent per-dataset counts, **not** a shared multiple-testing family "
        "and **not** a cross-dataset biological locus count -- Miyagi and Shumari "
        "coordinates are never combined or compared.",
        "",
        "## Candidate overview",
        "",
    ]
    for dataset in validated.datasets:
        lines.append(f"### `{dataset.dataset_id}`")
        lines.append("")
        if dataset.n_candidates == 0:
            lines.append(
                "No candidates passed Bonferroni or Benjamini-Hochberg significance in "
                "this dataset."
            )
            lines.append("")
            continue
        preview = _preview_candidates(dataset, limit=_EXEC_SUMMARY_PREVIEW_ROWS)
        lines.append("| Rank | Tier | Signal | Chr | Pos | p-value | beta | Reasons |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for row in preview:
            lines.append(
                f"| {row['candidate_rank']} | {row['priority_tier']} | `{row['signal_id']}` "
                f"| {row['chr']} | {row['pos']:,} | {row['pval']:.3g} | {row['beta']:.3g} "
                f"| {row['priority_reasons']} |"
            )
        if dataset.n_candidates > len(preview):
            links = _dataset_artifact_links(dataset)
            lines.append(
                f"\n_...and {dataset.n_candidates - len(preview)} more candidate(s); see "
                f"[{links['candidate_ranking']}]({links['candidate_ranking']})._"
            )
        lines.append("")
    lines += [
        "## What is not claimed here",
        "",
        "- A physical-distance signal is **not** an LD block and **not** an "
        "independently established QTL interval (this repository has no "
        "individual-level genotype data, so LD cannot be computed).",
        "- A signal's lead variant is **not** asserted to be the causal variant.",
        "- `priority_tier`/`priority_reasons` describe **downstream validation "
        "priority** derived only from Bonferroni/BH significance flags -- **not** "
        "biological importance, **not** a causal probability, and **not** a "
        "validated or experimentally confirmed breeding marker.",
        "- Miyagi and Shumari are different, non-interchangeable coordinate systems; "
        "no coordinate comparison or cross-reference ranking is performed anywhere "
        "in this report.",
        "",
        "## Recommended next actions",
        "",
        "- Downstream sequence confirmation of top-ranked candidates per dataset.",
        "- Reference-coordinate confirmation before comparing any candidate across "
        "reference genomes.",
        "- Review of neighboring variants around top-ranked candidates.",
        "- Marker-candidate design consideration for confirmed regions.",
        "- Wet-lab validation before treating any candidate as an established marker.",
        "",
        "_(These are suggested next steps, not features implemented by this "
        "repository -- see `analysis_report.md`'s Reproducibility section for "
        "tracked follow-up work.)_",
        "",
        "See `analysis_report.md` for full technical detail and "
        "`reproducibility/run_manifest.json` for the machine-readable audit record.",
    ]
    return "\n".join(lines) + "\n"


def build_analysis_report_markdown(validated: ValidatedAnalysisDir) -> str:
    """Build ``analysis_report.md``: the detailed technical report."""
    lines = [
        "# Analysis Report",
        "",
        "## 1. Analysis scope",
        "",
        "This report is a post-hoc re-analysis and organization of already-published "
        "GWAS summary statistics. **The underlying GWAS was not re-run**, and this "
        "report does not perform linkage-disequilibrium (LD) analysis, LD clumping, "
        "fine-mapping, liftover, or reference-coordinate harmonization.",
        "",
        "## 2. Input datasets and provenance",
        "",
        "| Dataset | Reference | Trait | Source SHA-256 |",
        "|---|---|---|---|",
    ]
    for dataset in validated.datasets:
        lines.append(
            f"| `{dataset.dataset_id}` | {dataset.reference} | {_trait_label(dataset.trait)} "
            f"| `{dataset.source_sha256}` |"
        )
    lines += [
        "",
        "Source checksums above are carried over verbatim from the batch analysis "
        "that produced this report's input (`batch_summary.tsv`); this report does not "
        "re-read or re-hash the original `.assoc.txt` files.",
        "",
        "## 3. Statistical diagnostics",
        "",
        f"- Bonferroni family-wise error rate (FWER) control: alpha = {validated.alpha} "
        "per dataset's own multiple-testing family.",
        f"- Benjamini-Hochberg (BH) false discovery rate (FDR) procedure: "
        f"fdr_level = {validated.fdr_level} per dataset's own family. BH-adjusted "
        "p-values are **not** Storey-style q-values.",
        "- Bonferroni controls the probability of *any* false positive across the "
        "family; BH controls the *expected proportion* of false positives among "
        "discoveries -- these are different guarantees for different purposes, not "
        "interchangeable synonyms.",
        "- The genomic inflation factor (lambda_GC) is reported strictly as a "
        "**diagnostic statistic**. A lambda_GC value alone does not establish or "
        "rule out population stratification, kinship, batch effects, or "
        "over-/under-correction as its cause, and no test statistic or p-value in "
        "this report is adjusted by lambda_GC.",
        "",
        "| Dataset | pvalue column | n_tests | Bonferroni discoveries | BH discoveries "
        "| lambda_GC |",
        "|---|---|---|---|---|---|",
    ]
    for dataset in validated.datasets:
        lines.append(
            f"| `{dataset.dataset_id}` | `{dataset.pvalue_column}` | {dataset.n_tests:,} "
            f"| {dataset.bonferroni_discoveries:,} | {dataset.bh_discoveries:,} "
            f"| {dataset.lambda_gc:.6g} |"
        )
    lines += [
        "",
        "## 4. Association signal organization",
        "",
        f"Association signals in this report are **physical-distance clusters only** "
        f"(`clustering_distance={validated.clustering_distance}` bp for this run): "
        "significant variants on the same chromosome, in the same dataset, are grouped "
        "together only when adjacent enough. This is **not** an LD block and **not** "
        "an independently established QTL interval -- this repository has no "
        "individual-level genotype data and cannot compute LD.",
        "",
        "| Dataset | Signals | Candidates |",
        "|---|---|---|",
    ]
    for dataset in validated.datasets:
        lines.append(
            f"| `{dataset.dataset_id}` | {dataset.n_signals:,} | {dataset.n_candidates:,} |"
        )
    lines += [
        "",
        "## 5. Candidate prioritization",
        "",
        "`priority_tier` is derived only from already-computed Bonferroni/BH "
        "significance flags (tier 1 = Bonferroni-significant, tier 2 = "
        "Benjamini-Hochberg-significant only) -- there is no additional weighting, "
        "scoring formula, or hidden parameter. This is a **downstream validation "
        "priority**, not a biological-importance ranking, not a causal probability, "
        "and not a validated or experimentally confirmed breeding marker.",
        "",
        "## 6. Dataset-specific results",
        "",
    ]
    for dataset in validated.datasets:
        links = _dataset_artifact_links(dataset)
        lines += [
            f"### `{dataset.dataset_id}` ({dataset.reference}, {_trait_label(dataset.trait)})",
            "",
            f"- n_tests: {dataset.n_tests:,}",
            f"- Bonferroni discoveries: {dataset.bonferroni_discoveries:,}",
            f"- BH discoveries: {dataset.bh_discoveries:,}",
            f"- lambda_GC (df={dataset.lambda_gc_df}): {dataset.lambda_gc:.6g}",
            f"- Physical-distance signals: {dataset.n_signals:,}",
            f"- Candidates: {dataset.n_candidates:,}",
            f"- Manhattan plot: [{links['manhattan']}]({links['manhattan']})",
            f"- QQ plot: [{links['qq']}]({links['qq']})",
            "",
        ]
        if dataset.n_candidates == 0:
            lines.append(
                "No candidates passed Bonferroni or Benjamini-Hochberg significance in "
                "this dataset."
            )
            lines.append("")
            continue
        preview = _preview_candidates(dataset, limit=_ANALYSIS_REPORT_PREVIEW_ROWS)
        lines.append("| Rank | Tier | Signal | Chr | Pos | p-value | beta | Reasons |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for row in preview:
            lines.append(
                f"| {row['candidate_rank']} | {row['priority_tier']} | `{row['signal_id']}` "
                f"| {row['chr']} | {row['pos']:,} | {row['pval']:.3g} | {row['beta']:.3g} "
                f"| {row['priority_reasons']} |"
            )
        if dataset.n_candidates > len(preview):
            lines.append(
                f"\n_...and {dataset.n_candidates - len(preview)} more; see "
                f"[{links['candidate_ranking']}]({links['candidate_ranking']})._"
            )
        lines.append("")

    lines += [
        "## 7. Figures",
        "",
        "Manhattan and QQ plots for each dataset are linked in Section 6 above and "
        "copied under `artifacts/<dataset_id>/` in this delivery package.",
        "",
        "## 8. Output artifact inventory",
        "",
        "See `reproducibility/run_manifest.json`'s `artifacts` array for the complete, "
        "checksummed inventory of every file in this delivery package.",
        "",
        "## 9. Scientific limitations",
        "",
        "- Association signals are physical-distance clusters only; this repository "
        "has no individual-level genotype data and cannot compute LD.",
        "- Lead variants are not asserted to be causal.",
        "- Candidates are not validated or experimentally confirmed breeding markers.",
        "- Bonferroni's `m` is the raw per-dataset variant count, not an LD-pruned "
        "effective-test count, so it may be conservative for LD-correlated SNPs.",
        "- The Benjamini-Hochberg guarantee holds under independence or positive "
        "regression dependency on a subset (PRDS); this repository has not "
        "demonstrated that PRDS holds for these correlated SNPs.",
        "- lambda_GC is a diagnostic value only and does not distinguish population "
        "stratification/kinship/batch effects from polygenicity as a cause of "
        "inflation.",
        "- Miyagi and Shumari are different, non-interchangeable coordinate systems; "
        "no cross-reference coordinate comparison or shared candidate ranking is "
        "performed.",
        "",
        "## 10. Reproducibility",
        "",
        "See `reproducibility/input_checksums.tsv` (per-dataset source checksums), "
        "`reproducibility/software_versions.json` (report-generation environment; "
        "the environment that generated the underlying batch/candidate artifacts is "
        "not recorded by those artifacts and is not claimed here), and "
        "`reproducibility/run_manifest.json` (full machine-readable audit record: "
        "parameters, per-dataset counts, scientific scope, and a checksummed "
        "artifact inventory).",
        "",
        "Follow-up work tracked separately, **not** implemented by this report: "
        "reference-genome coordinate contract, flanking-sequence extraction, ARMS "
        "marker candidate design, and GWAS execution on customer-provided "
        "individual-level data.",
    ]
    return "\n".join(lines) + "\n"


def build_run_manifest(
    validated: ValidatedAnalysisDir,
    *,
    provenance: dict[str, object],
    artifacts: list[dict[str, object]],
) -> dict[str, object]:
    """Build the ``run_manifest.json`` document (a plain, JSON-serializable ``dict``)."""
    return {
        "schema_version": RUN_MANIFEST_SCHEMA_VERSION,
        "analysis_type": "public_gwas_summary_statistics_post_hoc_batch_with_candidates",
        "source_batch_summary_schema_version": 2,
        "parameters": {
            "alpha": validated.alpha,
            "fdr_level": validated.fdr_level,
            "visualization_threshold": validated.visualization_threshold,
            "clustering_distance": validated.clustering_distance,
        },
        "scientific_scope": SCIENTIFIC_SCOPE,
        "datasets": [
            {
                "dataset_id": d.dataset_id,
                "reference": d.reference,
                "trait": d.trait,
                "source_sha256": d.source_sha256,
                "pvalue_column": d.pvalue_column,
                "n_tests": d.n_tests,
                "bonferroni_discoveries": d.bonferroni_discoveries,
                "bh_discoveries": d.bh_discoveries,
                "lambda_gc": d.lambda_gc,
                "lambda_gc_df": d.lambda_gc_df,
                "n_signals": d.n_signals,
                "n_candidates": d.n_candidates,
            }
            for d in validated.datasets
        ],
        "provenance": provenance,
        "artifacts": artifacts,
    }


def build_software_versions(
    *, report_generation_environment: dict[str, object]
) -> dict[str, object]:
    """Build the ``software_versions.json`` document.

    ``analysis_generation_environment`` is always the literal string
    :data:`UNAVAILABLE_FROM_SOURCE_ARTIFACTS`: the batch/candidate artifacts this report
    consumes do not themselves record the software versions or Git commit that produced
    them, so this function never substitutes the current (report-generation) environment
    for that unknown, and never guesses.
    """
    return {
        "schema_version": SOFTWARE_VERSIONS_SCHEMA_VERSION,
        "report_generation_environment": report_generation_environment,
        "analysis_generation_environment": UNAVAILABLE_FROM_SOURCE_ARTIFACTS,
        "note": (
            "The GWAS batch/candidate artifacts consumed by this report were not "
            "necessarily generated in this same process invocation, and this "
            "repository's existing batch/candidates outputs do not record their own "
            "generation-time software versions or Git commit. This section therefore "
            "describes only the environment that generated this report, not the "
            "environment that generated the underlying analysis artifacts."
        ),
    }
