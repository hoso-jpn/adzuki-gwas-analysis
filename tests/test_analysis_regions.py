"""Unit tests for adzuki_gwas_analysis.analysis.regions.

Pure in-memory DataFrame tests -- no file I/O.
"""

from __future__ import annotations

import unittest

import pandas as pd

from adzuki_gwas_analysis.analysis.regions import select_top_variant, subset_region
from adzuki_gwas_analysis.errors import EmptyRegionError


def _df(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows)


class SubsetRegionTests(unittest.TestCase):
    def test_filters_by_chrom_and_inclusive_bounds(self) -> None:
        df = _df(
            [
                {"chr": "Chr01", "pos": 100},
                {"chr": "Chr01", "pos": 200},
                {"chr": "Chr01", "pos": 300},
                {"chr": "Chr02", "pos": 200},
            ]
        )
        sub = subset_region(df, chrom="Chr01", start=100, end=200)
        self.assertEqual(list(sub["pos"]), [100, 200])

    def test_empty_result_for_no_matching_rows(self) -> None:
        df = _df([{"chr": "Chr01", "pos": 100}])
        sub = subset_region(df, chrom="Chr02", start=1, end=1000)
        self.assertTrue(sub.empty)


class SelectTopVariantTests(unittest.TestCase):
    def _rows(self) -> list[dict[str, object]]:
        return [
            {
                "chr": "Chr01",
                "pos": 100,
                "allele1": "A",
                "allele0": "G",
                "af": 0.3,
                "beta": 0.05,
                "pval": 1e-4,
            },
            {
                "chr": "Chr01",
                "pos": 200,
                "allele1": "C",
                "allele0": "T",
                "af": 0.4,
                "beta": -0.02,
                "pval": 1e-4,
            },
            {
                "chr": "Chr01",
                "pos": 300,
                "allele1": "G",
                "allele0": "A",
                "af": 0.2,
                "beta": 0.03,
                "pval": 5e-2,
            },
        ]

    def test_picks_smallest_pval(self) -> None:
        df = _df(self._rows())
        top = select_top_variant(
            df, dataset_id="d", region_id="r", chrom="Chr01", start=1, end=1000
        )
        self.assertEqual(top.pos, 100)
        self.assertEqual(top.pval, 1e-4)

    def test_tie_breaks_to_first_occurrence_in_row_order_not_smallest_pos(self) -> None:
        # Rows at pos=100 and pos=200 tie on pval (1e-4). idxmin's
        # first-occurrence semantics must return pos=100 (it comes first in
        # the DataFrame), not pos=200 (the smaller position).
        df = _df(self._rows())
        top = select_top_variant(
            df, dataset_id="d", region_id="r", chrom="Chr01", start=1, end=1000
        )
        self.assertEqual(top.pos, 100)

    def test_tie_break_respects_row_order_even_when_reversed(self) -> None:
        rows = list(reversed(self._rows()))
        df = _df(rows)
        top = select_top_variant(
            df, dataset_id="d", region_id="r", chrom="Chr01", start=1, end=1000
        )
        # Now the pos=200 tied row comes first in the (reversed) frame.
        self.assertEqual(top.pos, 200)

    def test_empty_region_raises_empty_region_error(self) -> None:
        df = _df(self._rows()).iloc[0:0]
        with self.assertRaises(EmptyRegionError) as ctx:
            select_top_variant(
                df, dataset_id="d", region_id="my_region", chrom="Chr01", start=1, end=1000
            )
        self.assertEqual(ctx.exception.region_id, "my_region")


if __name__ == "__main__":
    unittest.main()
