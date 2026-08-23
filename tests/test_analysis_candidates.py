"""Unit tests for adzuki_gwas_analysis.analysis.candidates.

Pure in-memory tests -- no file I/O. ``_significant_df`` builds inputs shaped exactly like
:func:`adzuki_gwas_analysis.analysis.diagnostics.build_significant_variants_table`'s return
value, without going through Bonferroni/BH computation itself (that is
tests/test_analysis_statistics.py's and tests/test_analysis_diagnostics.py's job).
"""

from __future__ import annotations

import unittest

import pandas as pd

from adzuki_gwas_analysis.analysis.candidates import (
    ASSOCIATION_PEAKS_COLUMNS,
    CANDIDATE_RANKING_COLUMNS,
    CANDIDATE_SNPS_COLUMNS,
    CANDIDATES_SCHEMA_VERSION,
    build_candidates_result,
    validate_clustering_distance,
)
from adzuki_gwas_analysis.analysis.diagnostics import SIGNIFICANT_VARIANTS_COLUMNS


def _significant_df(rows: list[dict[str, object]]) -> pd.DataFrame:
    """Build a significant_variants-shaped DataFrame from partial per-row overrides.

    Every row starts from a fixed set of defaults (all significant under both corrections)
    and only the fields a test cares about are overridden, keeping each test focused on the
    one thing it is checking.
    """
    defaults = {
        "chr": "Chr01",
        "pos": 1_000_000,
        "allele1": "A",
        "allele0": "G",
        "af": 0.3,
        "beta": 0.1,
        "pval": 1e-8,
        "pval_bonferroni": 1e-3,
        "pval_bh": 1e-3,
        "bonferroni_significant": True,
        "bh_significant": True,
    }
    built = [{**defaults, **row} for row in rows]
    return pd.DataFrame(built, columns=list(SIGNIFICANT_VARIANTS_COLUMNS))


_COMMON_KWARGS = {
    "dataset_id": "miyagi_water_permeability",
    "reference": "Miyagi",
    "trait": "water_permeability",
}


class ValidateClusteringDistanceTests(unittest.TestCase):
    def test_zero_is_allowed(self) -> None:
        validate_clustering_distance(0)

    def test_positive_is_allowed(self) -> None:
        validate_clustering_distance(50_000)

    def test_negative_raises(self) -> None:
        with self.assertRaises(ValueError):
            validate_clustering_distance(-1)

    def test_bool_raises(self) -> None:
        # bool is a subclass of int in Python; True/False must never silently pass as 1/0.
        with self.assertRaises(ValueError):
            validate_clustering_distance(True)  # type: ignore[arg-type]

    def test_float_raises(self) -> None:
        with self.assertRaises(ValueError):
            validate_clustering_distance(50_000.0)  # type: ignore[arg-type]


class InputContractTests(unittest.TestCase):
    def test_missing_required_column_raises(self) -> None:
        df = _significant_df([{}]).drop(columns=["beta"])
        with self.assertRaises(ValueError):
            build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)


class EmptyInputTests(unittest.TestCase):
    def test_zero_candidates_yields_header_only_tables(self) -> None:
        df = _significant_df([]).astype(
            {
                "chr": "object",
                "pos": "int64",
                "bonferroni_significant": "bool",
                "bh_significant": "bool",
            }
        )
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        self.assertEqual(result.n_signals, 0)
        self.assertEqual(result.n_candidates, 0)
        self.assertEqual(len(result.association_peaks), 0)
        self.assertEqual(len(result.candidate_snps), 0)
        self.assertEqual(len(result.candidate_ranking), 0)
        self.assertEqual(list(result.association_peaks.columns), list(ASSOCIATION_PEAKS_COLUMNS))
        self.assertEqual(list(result.candidate_snps.columns), list(CANDIDATE_SNPS_COLUMNS))
        self.assertEqual(list(result.candidate_ranking.columns), list(CANDIDATE_RANKING_COLUMNS))


class SchemaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.df = _significant_df([{"pos": 1_000_000}, {"pos": 1_000_500}])
        self.result = build_candidates_result(self.df, clustering_distance=1000, **_COMMON_KWARGS)

    def test_association_peaks_header(self) -> None:
        self.assertEqual(
            list(self.result.association_peaks.columns), list(ASSOCIATION_PEAKS_COLUMNS)
        )

    def test_candidate_snps_header(self) -> None:
        self.assertEqual(list(self.result.candidate_snps.columns), list(CANDIDATE_SNPS_COLUMNS))

    def test_candidate_ranking_header(self) -> None:
        self.assertEqual(
            list(self.result.candidate_ranking.columns), list(CANDIDATE_RANKING_COLUMNS)
        )

    def test_schema_version_is_stamped_everywhere(self) -> None:
        self.assertTrue((self.result.association_peaks["schema_version"] == 1).all())
        self.assertTrue((self.result.candidate_snps["schema_version"] == 1).all())
        self.assertTrue((self.result.candidate_ranking["schema_version"] == 1).all())
        self.assertEqual(CANDIDATES_SCHEMA_VERSION, 1)

    def test_dataset_reference_trait_carried_on_every_row(self) -> None:
        for table in (
            self.result.association_peaks,
            self.result.candidate_snps,
            self.result.candidate_ranking,
        ):
            self.assertTrue((table["dataset_id"] == "miyagi_water_permeability").all())
            self.assertTrue((table["reference"] == "Miyagi").all())
            self.assertTrue((table["trait"] == "water_permeability").all())


class SingleCandidateTests(unittest.TestCase):
    def test_one_variant_makes_one_signal_and_one_candidate(self) -> None:
        df = _significant_df([{"pos": 5_000_000}])
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        self.assertEqual(result.n_signals, 1)
        self.assertEqual(result.n_candidates, 1)
        self.assertEqual(result.candidate_ranking["candidate_rank"].tolist(), [1])
        self.assertTrue(result.candidate_snps["is_lead_variant"].iloc[0])
        peak = result.association_peaks.iloc[0]
        self.assertEqual(peak["n_significant_variants"], 1)
        self.assertEqual(peak["start"], 5_000_000)
        self.assertEqual(peak["end"], 5_000_000)


class ClusteringBoundaryTests(unittest.TestCase):
    def test_gap_exactly_at_distance_merges(self) -> None:
        df = _significant_df([{"pos": 1_000_000}, {"pos": 1_001_000}])
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        self.assertEqual(result.n_signals, 1)
        self.assertEqual(result.association_peaks.iloc[0]["n_significant_variants"], 2)

    def test_gap_one_over_distance_does_not_merge(self) -> None:
        df = _significant_df([{"pos": 1_000_000}, {"pos": 1_001_001}])
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        self.assertEqual(result.n_signals, 2)

    def test_distance_zero_merges_only_identical_positions(self) -> None:
        df = _significant_df(
            [
                {"pos": 2_000_000, "allele1": "A", "allele0": "G"},
                {"pos": 2_000_000, "allele1": "C", "allele0": "T"},
                {"pos": 2_000_001, "allele1": "A", "allele0": "G"},
            ]
        )
        result = build_candidates_result(df, clustering_distance=0, **_COMMON_KWARGS)
        self.assertEqual(result.n_signals, 2)

    def test_chained_merging_spans_a_run_of_adjacent_gaps(self) -> None:
        # 100 -> 100+500=600 (gap 500, merges) -> 600+500=1100 (gap 500, merges): one signal
        # spanning 1000bp even though the first and last variants are 1000bp apart, further
        # than clustering_distance=500 from each other directly.
        df = _significant_df([{"pos": 100}, {"pos": 600}, {"pos": 1100}])
        result = build_candidates_result(df, clustering_distance=500, **_COMMON_KWARGS)
        self.assertEqual(result.n_signals, 1)
        peak = result.association_peaks.iloc[0]
        self.assertEqual(peak["start"], 100)
        self.assertEqual(peak["end"], 1100)
        self.assertEqual(peak["n_significant_variants"], 3)

    def test_never_clusters_across_chromosomes(self) -> None:
        df = _significant_df([{"chr": "Chr01", "pos": 1000}, {"chr": "Chr02", "pos": 1000}])
        result = build_candidates_result(df, clustering_distance=10_000_000, **_COMMON_KWARGS)
        self.assertEqual(result.n_signals, 2)
        self.assertEqual(set(result.association_peaks["chromosome"]), {"Chr01", "Chr02"})

    def test_signal_ids_use_natural_chromosome_order_not_lexicographic(self) -> None:
        # Chr2 must sort before Chr10 (lexicographic order would put Chr10 first).
        df = _significant_df([{"chr": "Chr10", "pos": 1}, {"chr": "Chr2", "pos": 1}])
        result = build_candidates_result(df, clustering_distance=0, **_COMMON_KWARGS)
        ordered_chroms = result.association_peaks.sort_values("signal_id")["chromosome"].tolist()
        self.assertEqual(ordered_chroms, ["Chr2", "Chr10"])


class LeadVariantTieBreakTests(unittest.TestCase):
    def test_smallest_pval_wins(self) -> None:
        df = _significant_df(
            [
                {"pos": 100, "pval": 1e-4, "beta": 0.1},
                {"pos": 200, "pval": 1e-8, "beta": 0.05},
            ]
        )
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        peak = result.association_peaks.iloc[0]
        self.assertEqual(peak["lead_pos"], 200)

    def test_tie_on_pval_breaks_by_larger_absolute_beta(self) -> None:
        df = _significant_df(
            [
                {"pos": 100, "pval": 1e-8, "beta": -0.05},
                {"pos": 200, "pval": 1e-8, "beta": 0.9},
            ]
        )
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        peak = result.association_peaks.iloc[0]
        self.assertEqual(peak["lead_pos"], 200)

    def test_tie_on_pval_and_beta_breaks_by_smallest_pos(self) -> None:
        df = _significant_df(
            [
                {"pos": 300, "pval": 1e-8, "beta": 0.5},
                {"pos": 100, "pval": 1e-8, "beta": -0.5},
            ]
        )
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        peak = result.association_peaks.iloc[0]
        self.assertEqual(peak["lead_pos"], 100)

    def test_full_tie_breaks_by_original_source_row_order(self) -> None:
        # Two multi-allelic rows at the exact same position, identical pval and |beta|:
        # the row that appeared first in the source file must win.
        df = _significant_df(
            [
                {"pos": 100, "pval": 1e-8, "beta": 0.5, "allele1": "A", "allele0": "G"},
                {"pos": 100, "pval": 1e-8, "beta": -0.5, "allele1": "C", "allele0": "T"},
            ]
        )
        result = build_candidates_result(df, clustering_distance=0, **_COMMON_KWARGS)
        peak = result.association_peaks.iloc[0]
        self.assertEqual((peak["lead_allele1"], peak["lead_allele0"]), ("A", "G"))


class OnlyOneLeadVariantPerSignalTests(unittest.TestCase):
    def test_exactly_one_row_flagged_as_lead_per_signal(self) -> None:
        df = _significant_df([{"pos": 100}, {"pos": 200}, {"pos": 5_000_000}])
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        lead_counts = (
            result.candidate_snps[result.candidate_snps["is_lead_variant"]]
            .groupby("signal_id")
            .size()
        )
        self.assertTrue((lead_counts == 1).all())
        self.assertEqual(len(lead_counts), result.n_signals)


class PriorityTierTests(unittest.TestCase):
    def test_bonferroni_significant_is_tier_one(self) -> None:
        df = _significant_df([{"bonferroni_significant": True, "bh_significant": True}])
        result = build_candidates_result(df, clustering_distance=0, **_COMMON_KWARGS)
        self.assertEqual(result.candidate_ranking["priority_tier"].iloc[0], 1)
        self.assertIn(
            "bonferroni_significant", result.candidate_ranking["priority_reasons"].iloc[0]
        )

    def test_bh_only_is_tier_two(self) -> None:
        df = _significant_df([{"bonferroni_significant": False, "bh_significant": True}])
        result = build_candidates_result(df, clustering_distance=0, **_COMMON_KWARGS)
        self.assertEqual(result.candidate_ranking["priority_tier"].iloc[0], 2)
        reasons = result.candidate_ranking["priority_reasons"].iloc[0]
        self.assertIn("bh_significant", reasons)
        self.assertNotIn("bonferroni_significant", reasons)

    def test_tier_one_ranked_before_tier_two_regardless_of_pval(self) -> None:
        df = _significant_df(
            [
                {
                    "pos": 100,
                    "pval": 1e-2,
                    "bonferroni_significant": True,
                    "bh_significant": True,
                },
                {
                    "pos": 5_000_000,
                    "pval": 1e-12,
                    "bonferroni_significant": False,
                    "bh_significant": True,
                },
            ]
        )
        result = build_candidates_result(df, clustering_distance=0, **_COMMON_KWARGS)
        ranking = result.candidate_ranking.sort_values("candidate_rank")
        self.assertEqual(ranking["pos"].tolist(), [100, 5_000_000])
        self.assertEqual(ranking["priority_tier"].tolist(), [1, 2])

    def test_multi_variant_signal_reason_vs_single_variant_signal_reason(self) -> None:
        df = _significant_df([{"pos": 100}, {"pos": 200}, {"pos": 9_000_000}])
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        by_pos = result.candidate_ranking.set_index("pos")
        self.assertIn("member_of_multi_variant_signal", by_pos.loc[100, "priority_reasons"])
        self.assertIn("member_of_multi_variant_signal", by_pos.loc[200, "priority_reasons"])
        self.assertIn("single_variant_signal", by_pos.loc[9_000_000, "priority_reasons"])


class RankingCompletenessTests(unittest.TestCase):
    def test_rank_is_a_dense_sequence_starting_at_one(self) -> None:
        df = _significant_df([{"pos": p} for p in (100, 200, 300, 9_000_000, 9_000_100)])
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        self.assertEqual(
            sorted(result.candidate_ranking["candidate_rank"].tolist()),
            list(range(1, len(result.candidate_ranking) + 1)),
        )

    def test_ranking_covers_every_candidate_exactly_once(self) -> None:
        df = _significant_df([{"pos": p} for p in (100, 200, 300)])
        result = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        self.assertEqual(len(result.candidate_ranking), result.n_candidates)
        self.assertEqual(
            set(zip(result.candidate_ranking["pos"], strict=False)),
            {(100,), (200,), (300,)},
        )


class DeterminismTests(unittest.TestCase):
    def test_repeated_calls_produce_identical_tables(self) -> None:
        df = _significant_df([{"pos": p} for p in (500, 100, 9_000_000, 300)])
        first = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        second = build_candidates_result(df, clustering_distance=1000, **_COMMON_KWARGS)
        pd.testing.assert_frame_equal(first.association_peaks, second.association_peaks)
        pd.testing.assert_frame_equal(first.candidate_snps, second.candidate_snps)
        pd.testing.assert_frame_equal(first.candidate_ranking, second.candidate_ranking)


class MultipleSignalsTests(unittest.TestCase):
    def test_association_peaks_ordered_by_genome_position(self) -> None:
        df = _significant_df(
            [
                {"chr": "Chr02", "pos": 500},
                {"chr": "Chr01", "pos": 9_000_000},
                {"chr": "Chr01", "pos": 100},
            ]
        )
        result = build_candidates_result(df, clustering_distance=0, **_COMMON_KWARGS)
        self.assertEqual(result.n_signals, 3)
        self.assertEqual(
            list(
                zip(
                    result.association_peaks["chromosome"],
                    result.association_peaks["start"],
                    strict=True,
                )
            ),
            [("Chr01", 100), ("Chr01", 9_000_000), ("Chr02", 500)],
        )


if __name__ == "__main__":
    unittest.main()
