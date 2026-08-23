"""Physical-distance association-signal clustering and downstream validation priority.

This module is a pure, in-memory consumer of the same significant-variant population
:func:`adzuki_gwas_analysis.analysis.diagnostics.build_significant_variants_table` already
computes for ``significant_variants.tsv`` (Bonferroni FWER discoveries union
Benjamini-Hochberg FDR discoveries, for one dataset's own multiple-testing family). It never
re-reads a raw ``.assoc.txt`` file and never re-runs Bonferroni/BH itself -- it only groups
and prioritizes rows a caller already validated and corrected.

Scientific scope, deliberately narrow:

* Grouping is **physical-distance clustering only** -- adjacent significant variants on the
  same chromosome, in the same dataset, within a caller-supplied ``clustering_distance`` (bp).
  There is no individual-level genotype data in this repository, so no linkage disequilibrium
  (LD) can be computed. A cluster produced here is not an LD block, not an independently
  defined QTL interval, and adjacency does not establish (or rule out) that two variants
  reflect the same underlying signal.
* A cluster's "lead variant" is the one variant this module's deterministic tie-break rule
  picks first -- it is not asserted to be causal, and no other variant in the cluster is
  asserted to be non-causal or merely a proxy for it.
* ``priority_tier``/``priority_reasons`` describe **downstream validation priority** --
  which candidates to look at first -- built only from already-computed, already-explained
  quantities (the Bonferroni/BH significance flags this dataset's own family already
  established). This is not a biological-importance ranking, not a probability of being a
  true positive, and not a validated or experimentally confirmed breeding marker.
* ``dataset_id``/``reference``/``trait`` are carried on every output row and never combined,
  compared, or re-clustered across datasets: Miyagi and Shumari are different coordinate
  systems, and the 3 traits are different phenotypes. A cluster never spans more than one
  chromosome or more than one dataset.

There is no scientifically justified default for ``clustering_distance`` anywhere in this
repository (the existing ``config/*_regions.toml`` windows are human-picked, post-hoc
visualization windows -- see :mod:`adzuki_gwas_analysis.analysis.region_config` -- not an LD
or QTL-interval estimate), so callers must always supply it explicitly; this module never
assumes a value.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from adzuki_gwas_analysis.analysis.chromosomes import order_chromosomes
from adzuki_gwas_analysis.analysis.diagnostics import SIGNIFICANT_VARIANTS_COLUMNS

#: Bumped only if any of the 3 output tables' column set/meaning changes.
CANDIDATES_SCHEMA_VERSION = 1

#: Self-describing clustering rule, recorded verbatim in every association_peaks.tsv row so
#: the file never needs to be interpreted alongside separate documentation to know what
#: "signal" means for that row.
CLUSTERING_METHOD = (
    "physical_distance_1d: within one dataset_id and one chromosome, significant variants "
    "sorted by pos are merged into the same signal while the gap to the nearest "
    "already-clustered variant is <= clustering_distance bp; not an LD block or an "
    "independently defined QTL interval"
)

ASSOCIATION_PEAKS_COLUMNS: tuple[str, ...] = (
    "schema_version",
    "dataset_id",
    "reference",
    "trait",
    "signal_id",
    "chromosome",
    "start",
    "end",
    "n_significant_variants",
    "lead_pos",
    "lead_allele1",
    "lead_allele0",
    "lead_pval",
    "lead_pval_bonferroni",
    "lead_pval_bh",
    "clustering_method",
    "clustering_distance",
)

CANDIDATE_SNPS_COLUMNS: tuple[str, ...] = (
    "schema_version",
    "dataset_id",
    "reference",
    "trait",
    "signal_id",
    "chr",
    "pos",
    "allele1",
    "allele0",
    "af",
    "beta",
    "pval",
    "pval_bonferroni",
    "pval_bh",
    "bonferroni_significant",
    "bh_significant",
    "is_lead_variant",
)

CANDIDATE_RANKING_COLUMNS: tuple[str, ...] = (
    "schema_version",
    "dataset_id",
    "reference",
    "trait",
    "signal_id",
    "chr",
    "pos",
    "allele1",
    "allele0",
    "candidate_rank",
    "priority_tier",
    "priority_reasons",
)

#: Every column build_significant_variants_table() guarantees is present on its input.
_REQUIRED_INPUT_COLUMNS = frozenset(SIGNIFICANT_VARIANTS_COLUMNS)


def validate_clustering_distance(clustering_distance: int) -> None:
    """Raise ``ValueError`` unless ``clustering_distance`` is a non-negative integer (bp).

    ``0`` is a legal, meaningful value: it merges only variants sharing the exact same
    ``pos`` (e.g. distinct multi-allelic rows at one site), never anything else.
    """
    if isinstance(clustering_distance, bool) or not isinstance(clustering_distance, int):
        raise ValueError(
            f"clustering_distance must be a plain int (base pairs), "
            f"got {clustering_distance!r} ({type(clustering_distance).__name__})"
        )
    if clustering_distance < 0:
        raise ValueError(
            f"clustering_distance must be >= 0 (base pairs), got {clustering_distance}"
        )


@dataclass(frozen=True, slots=True)
class CandidatesResult:
    """The 3 candidate-extraction output tables for one dataset, plus headline counts."""

    association_peaks: pd.DataFrame
    candidate_snps: pd.DataFrame
    candidate_ranking: pd.DataFrame
    n_signals: int
    n_candidates: int


def _empty_result() -> CandidatesResult:
    return CandidatesResult(
        association_peaks=pd.DataFrame(columns=list(ASSOCIATION_PEAKS_COLUMNS)),
        candidate_snps=pd.DataFrame(columns=list(CANDIDATE_SNPS_COLUMNS)),
        candidate_ranking=pd.DataFrame(columns=list(CANDIDATE_RANKING_COLUMNS)),
        n_signals=0,
        n_candidates=0,
    )


def _tie_break_sort(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by (pval asc, |beta| desc, pos asc, source_row_order asc), stably.

    ``pval`` (the manifest's primary likelihood-ratio-test column) is the sole significance
    criterion. On a tie, the variant with the larger-magnitude ``beta`` -- a real, always-
    finite effect-size estimate schema v1 already guarantees for every row (see
    ``README.md``'s "Beta values represent effect size estimates..." note) -- sorts first,
    since it is the more informative of the two remaining candidates at equal significance.
    Further ties fall back to genomic position, then to the row's original position in the
    source ``.assoc.txt`` file (``source_row_order``), so the order is always fully
    deterministic even when every other field matches exactly (e.g. two multi-allelic
    records at the same ``pos``).

    This deliberately differs from
    :func:`adzuki_gwas_analysis.analysis.regions.select_top_variant`, which (via
    ``Series.idxmin``) breaks ties purely by first-occurrence file order: that function picks
    one representative point to plot inside a human-picked, ad-hoc visualization window, with
    no candidate-prioritization claim attached to the choice. A cluster's *lead variant* here
    is presented to a breeder as the headline candidate for that signal, so an
    effect-size-based tie-break is used instead of leaving the choice to incidental file
    order.
    """
    working = df.copy()
    working["_abs_beta"] = working["beta"].abs()
    return working.sort_values(
        by=["pval", "_abs_beta", "pos", "source_row_order"],
        ascending=[True, False, True, True],
        kind="mergesort",
    ).drop(columns=["_abs_beta"])


def _cluster_one_chromosome(chrom_df: pd.DataFrame, *, clustering_distance: int) -> pd.Series:
    """Return a ``cluster_index`` (0-based, local to this chromosome) per row of ``chrom_df``.

    ``chrom_df`` must already be sorted by ``pos`` ascending. A new cluster starts whenever
    the gap between a row's ``pos`` and the running maximum ``pos`` of the currently open
    cluster exceeds ``clustering_distance`` -- i.e. clusters are chains of variants each
    within ``clustering_distance`` of some other already-clustered variant, not merely within
    ``clustering_distance`` of the cluster's first variant.
    """
    positions = chrom_df["pos"].to_numpy()
    cluster_index = np.zeros(len(positions), dtype="int64")
    current_cluster = 0
    cluster_max_pos = positions[0]
    for i in range(1, len(positions)):
        if positions[i] - cluster_max_pos > clustering_distance:
            current_cluster += 1
        cluster_max_pos = max(cluster_max_pos, positions[i])
        cluster_index[i] = current_cluster
    return pd.Series(cluster_index, index=chrom_df.index)


def build_candidates_result(
    significant_df: pd.DataFrame,
    *,
    dataset_id: str,
    reference: str,
    trait: str,
    clustering_distance: int,
) -> CandidatesResult:
    """Cluster ``significant_df`` into signals and build all 3 candidate output tables.

    ``significant_df`` must have exactly the columns
    :data:`adzuki_gwas_analysis.analysis.diagnostics.SIGNIFICANT_VARIANTS_COLUMNS`, in the
    row order :func:`~adzuki_gwas_analysis.analysis.diagnostics.build_significant_variants_table`
    produces (the source file's own row order, restricted to Bonferroni-or-BH discoveries) --
    that row order is used as this function's own final tie-break, via an internal
    ``source_row_order`` column derived from ``significant_df.index`` (0-based, matching
    ``reset_index(drop=True)``'s guarantee on that function's return value).

    Returns header-only, zero-row tables (never omitted files) when ``significant_df`` is
    empty. Every returned row carries ``dataset_id``/``reference``/``trait`` verbatim and a
    cluster never spans more than one chromosome -- clustering and ranking are always scoped
    to this one dataset's own significant-variant population, never combined with any other
    dataset's.
    """
    validate_clustering_distance(clustering_distance)

    missing = _REQUIRED_INPUT_COLUMNS - set(significant_df.columns)
    if missing:
        raise ValueError(
            f"significant_df is missing required column(s) {sorted(missing)}; expected "
            f"exactly the columns build_significant_variants_table() produces"
        )

    if significant_df.empty:
        return _empty_result()

    working = significant_df.reset_index(drop=True).copy()
    working["source_row_order"] = np.arange(len(working), dtype="int64")

    chrom_order = order_chromosomes(working["chr"].astype(str).unique().tolist())

    cluster_frames: list[pd.DataFrame] = []
    signal_counter = 0
    peak_rows: list[dict[str, object]] = []

    for chrom in chrom_order:
        chrom_df = working[working["chr"].astype(str) == chrom].sort_values(
            by="pos", kind="mergesort"
        )
        local_cluster = _cluster_one_chromosome(chrom_df, clustering_distance=clustering_distance)
        for local_id in range(int(local_cluster.max()) + 1):
            signal_counter += 1
            signal_id = f"{dataset_id}_peak_{signal_counter:04d}"
            member_index = local_cluster.index[local_cluster == local_id]
            members = working.loc[member_index]
            lead = _tie_break_sort(members).iloc[0]

            annotated = members.copy()
            annotated["signal_id"] = signal_id
            annotated["is_lead_variant"] = (
                (annotated["pos"] == lead["pos"])
                & (annotated["allele1"] == lead["allele1"])
                & (annotated["allele0"] == lead["allele0"])
                & (annotated["source_row_order"] == lead["source_row_order"])
            )
            cluster_frames.append(annotated)

            peak_rows.append(
                {
                    "schema_version": CANDIDATES_SCHEMA_VERSION,
                    "dataset_id": dataset_id,
                    "reference": reference,
                    "trait": trait,
                    "signal_id": signal_id,
                    "chromosome": chrom,
                    "start": int(members["pos"].min()),
                    "end": int(members["pos"].max()),
                    "n_significant_variants": int(len(members)),
                    "lead_pos": int(lead["pos"]),
                    "lead_allele1": lead["allele1"],
                    "lead_allele0": lead["allele0"],
                    "lead_pval": float(lead["pval"]),
                    "lead_pval_bonferroni": float(lead["pval_bonferroni"]),
                    "lead_pval_bh": float(lead["pval_bh"]),
                    "clustering_method": CLUSTERING_METHOD,
                    "clustering_distance": clustering_distance,
                }
            )

    association_peaks = pd.DataFrame(peak_rows, columns=list(ASSOCIATION_PEAKS_COLUMNS))

    annotated_all = pd.concat(cluster_frames, axis=0)
    # Row order: genome order (the same chromosome/pos scan clustering itself used), not the
    # original source-file order -- grouping adjacent candidates together is more useful to a
    # reader than the incidental order rows happened to appear in the raw .assoc.txt file.
    annotated_all = annotated_all.sort_values(
        by=["signal_id", "pos", "source_row_order"], kind="mergesort"
    )

    candidate_snps = annotated_all.copy()
    candidate_snps["schema_version"] = CANDIDATES_SCHEMA_VERSION
    candidate_snps["dataset_id"] = dataset_id
    candidate_snps["reference"] = reference
    candidate_snps["trait"] = trait
    candidate_snps = candidate_snps.loc[:, list(CANDIDATE_SNPS_COLUMNS)].reset_index(drop=True)

    n_members_by_signal = annotated_all.groupby("signal_id", sort=False).size()
    priority_tier = np.where(annotated_all["bonferroni_significant"], 1, 2)

    reasons_col = []
    for _, row in annotated_all.iterrows():
        reasons = []
        if row["bonferroni_significant"]:
            reasons.append("bonferroni_significant")
        if row["bh_significant"]:
            reasons.append("bh_significant")
        if row["is_lead_variant"]:
            reasons.append("lead_variant_of_signal")
        if n_members_by_signal[row["signal_id"]] > 1:
            reasons.append("member_of_multi_variant_signal")
        else:
            reasons.append("single_variant_signal")
        reasons_col.append(";".join(reasons))

    ranking = annotated_all.copy()
    ranking["schema_version"] = CANDIDATES_SCHEMA_VERSION
    ranking["dataset_id"] = dataset_id
    ranking["reference"] = reference
    ranking["trait"] = trait
    ranking["priority_tier"] = priority_tier
    ranking["priority_reasons"] = reasons_col

    # Rank across the whole dataset's candidate population: priority_tier (Bonferroni-tier
    # candidates first, then BH-only candidates) is the explicit primary sort key -- never
    # left to numeric coincidence -- followed by the same deterministic tie-break used for
    # lead-variant selection. Never combined with any other dataset's candidates.
    ranking["_abs_beta"] = ranking["beta"].abs()
    ranking = ranking.sort_values(
        by=["priority_tier", "pval", "_abs_beta", "pos", "source_row_order"],
        ascending=[True, True, False, True, True],
        kind="mergesort",
    ).drop(columns=["_abs_beta"])
    ranking = ranking.reset_index(drop=True)
    ranking.insert(0, "candidate_rank", np.arange(1, len(ranking) + 1))
    candidate_ranking = ranking.loc[:, list(CANDIDATE_RANKING_COLUMNS)].reset_index(drop=True)

    return CandidatesResult(
        association_peaks=association_peaks,
        candidate_snps=candidate_snps,
        candidate_ranking=candidate_ranking,
        n_signals=len(association_peaks),
        n_candidates=len(candidate_snps),
    )
