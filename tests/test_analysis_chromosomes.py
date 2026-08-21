"""Unit tests for adzuki_gwas_analysis.analysis.chromosomes.

Pure in-memory DataFrame tests -- no file I/O -- covering chromosome
ordering (including the ``Chr2``/``Chr10`` case a lexicographic sort gets
wrong) and Manhattan x-coordinate/offset computation.
"""

from __future__ import annotations

import unittest

import pandas as pd

from adzuki_gwas_analysis.analysis.chromosomes import (
    compute_manhattan_coordinates,
    natural_chromosome_key,
    order_chromosomes,
)


class NaturalChromosomeKeyTests(unittest.TestCase):
    def test_chr2_sorts_before_chr10(self) -> None:
        # A plain lexicographic sort puts "Chr10" before "Chr2"; natural
        # order must not.
        self.assertEqual(order_chromosomes(["Chr10", "Chr2", "Chr1"]), ["Chr1", "Chr2", "Chr10"])

    def test_zero_padded_current_data_order_is_unaffected(self) -> None:
        # scripts/01_manhattan_plot.py used lexicographic sorted() on the
        # real zero-padded Chr01..Chr11 data; natural order must reproduce
        # that same order for this repository's current data.
        padded = [f"Chr{n:02d}" for n in range(11, 0, -1)]
        self.assertEqual(order_chromosomes(padded), sorted(padded))

    def test_no_trailing_digits_sorts_after_numbered_chromosomes(self) -> None:
        self.assertEqual(order_chromosomes(["ChrUn", "Chr2", "Chr1"]), ["Chr1", "Chr2", "ChrUn"])

    def test_key_is_a_total_order_key(self) -> None:
        self.assertLess(natural_chromosome_key("Chr2"), natural_chromosome_key("Chr10"))
        self.assertLess(natural_chromosome_key("Chr01"), natural_chromosome_key("Chr02"))


class ComputeManhattanCoordinatesTests(unittest.TestCase):
    def test_offsets_and_ticks_match_legacy_accumulation(self) -> None:
        df = pd.DataFrame(
            {
                "chr": ["Chr01", "Chr01", "Chr02", "Chr02"],
                "pos": [100, 300, 50, 200],
                "pval": [1e-3, 1e-4, 1e-2, 1e-5],
            }
        )
        coords = compute_manhattan_coordinates(df)

        self.assertEqual(coords.chrom_order, ["Chr01", "Chr02"])
        self.assertEqual(coords.variant_count, 4)
        self.assertEqual(coords.chrom_variant_counts, {"Chr01": 2, "Chr02": 2})

        # Chr01 has no offset (first chromosome); Chr02's offset is Chr01's
        # max pos (300), matching `offset += sub["pos"].max()` in the
        # original script.
        expected_x = [100, 300, 50 + 300, 200 + 300]
        self.assertEqual(list(coords.x), expected_x)

        self.assertEqual(coords.tick_labels, ["01", "02"])
        self.assertEqual(coords.tick_positions[0], pd.Series([100, 300]).median())
        self.assertEqual(coords.tick_positions[1], pd.Series([350, 500]).median())

    def test_single_chromosome_has_zero_offset(self) -> None:
        df = pd.DataFrame({"chr": ["Chr05", "Chr05"], "pos": [10, 20], "pval": [0.1, 0.2]})
        coords = compute_manhattan_coordinates(df)
        self.assertEqual(list(coords.x), [10, 20])


if __name__ == "__main__":
    unittest.main()
