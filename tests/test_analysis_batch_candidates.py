"""Unit tests for candidate-extraction integration in adzuki_gwas_analysis.analysis.batch.

A separate file from tests/test_analysis_batch.py / tests/test_analysis_batch_failures.py
(left unmodified by this Issue, aside from the pre-existing subcommand-count pin in
tests/test_analysis_cli_batch.py) so those modules' diffs stay empty. Every test uses tiny
synthetic fixtures under a TemporaryDirectory -- none reads data/raw/.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from adzuki_gwas_analysis.analysis import batch
from tests.analysis_support import write_full_manifest, write_six_dataset_files

_DATASET_IDS_IN_MANIFEST_ORDER: tuple[str, ...] = (
    "miyagi_water_permeability",
    "miyagi_red_seedcoat",
    "miyagi_mottled_black_seedcoat",
    "shumari_water_permeability",
    "shumari_red_seedcoat",
    "shumari_mottled_black_seedcoat",
)

_ROW_TEMPLATE = (
    "Chr01\t.\t{pos1}\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t{p1}\t1e-3",
    "Chr01\t.\t{pos2}\t0\tC\tT\t0.45\t-0.02\t0.02\t71.0\t31.0\t29.0\t5e-2\t{p2}\t3e-2",
    "Chr02\t.\t{pos3}\t1\tT\tA\t0.10\t0.11\t0.03\t72.0\t32.0\t30.0\t9e-3\t{p3}\t7e-3",
)


def _six_dataset_rows() -> dict[str, list[str]]:
    """3 rows per dataset, each with at least one Bonferroni-significant p-value.

    Distinct ``pos`` values per dataset (offset by index) so a candidate file's content can
    be checked against the dataset directory it was written under.
    """
    rows_by_id: dict[str, list[str]] = {}
    for i, dataset_id in enumerate(_DATASET_IDS_IN_MANIFEST_ORDER):
        rows_by_id[dataset_id] = [
            row.format(
                pos1=1_000_000 + i,
                pos2=1_000_500 + i,
                pos3=9_000_000 + i,
                p1="1e-9",
                p2="2e-9",
                p3="9e-1",
            )
            for row in _ROW_TEMPLATE
        ]
    return rows_by_id


class BatchCandidatesTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.data_dir = self.tmp_path / "data"
        self.data_dir.mkdir()
        self.manifest_path = self.tmp_path / "manifest.toml"
        self.output_dir = self.tmp_path / "out"

    def _write_all_six_valid(self) -> None:
        dataset_paths = write_six_dataset_files(self.data_dir, _six_dataset_rows())
        write_full_manifest(self.manifest_path, dataset_paths)


class BackwardCompatibilityTests(BatchCandidatesTestCase):
    def test_omitting_clustering_distance_produces_exactly_twenty_five_files(self) -> None:
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        all_files = [p for p in self.output_dir.rglob("*") if p.is_file()]
        self.assertEqual(len(all_files), 25)

    def test_omitting_clustering_distance_keeps_schema_version_one(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        summary = pd.read_csv(outcome.summary_path, sep="\t")
        self.assertTrue((summary["schema_version"] == 1).all())
        for candidate_column in (
            "n_signals",
            "n_candidates",
            "clustering_distance",
            "association_peaks_path",
            "candidate_snps_path",
            "candidate_ranking_path",
        ):
            self.assertNotIn(candidate_column, summary.columns)

    def test_omitting_clustering_distance_leaves_dataset_result_fields_none(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        for result in outcome.datasets:
            self.assertIsNone(result.n_signals)
            self.assertIsNone(result.n_candidates)
            self.assertIsNone(result.association_peaks_path)


class CandidateExtractionEnabledTests(BatchCandidatesTestCase):
    def test_supplying_clustering_distance_produces_forty_three_files(self) -> None:
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        all_files = [p for p in self.output_dir.rglob("*") if p.is_file()]
        self.assertEqual(len(all_files), 43)

    def test_each_dataset_directory_has_exactly_seven_files(self) -> None:
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        for dataset_id in _DATASET_IDS_IN_MANIFEST_ORDER:
            files = {p.name for p in (self.output_dir / dataset_id).iterdir() if p.is_file()}
            self.assertEqual(
                files,
                {
                    f"{dataset_id}_manhattan.png",
                    f"{dataset_id}_qq.png",
                    "statistical_diagnostics.tsv",
                    "significant_variants.tsv",
                    "association_peaks.tsv",
                    "candidate_snps.tsv",
                    "candidate_ranking.tsv",
                },
            )

    def test_batch_summary_schema_version_is_two_and_has_candidate_columns(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        summary = pd.read_csv(outcome.summary_path, sep="\t")
        self.assertTrue((summary["schema_version"] == 2).all())
        self.assertTrue((summary["clustering_distance"] == 1000).all())
        for candidate_column in (
            "n_signals",
            "n_candidates",
            "association_peaks_path",
            "candidate_snps_path",
            "candidate_ranking_path",
        ):
            self.assertIn(candidate_column, summary.columns)

    def test_candidates_are_never_pooled_across_datasets(self) -> None:
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        for dataset_id in _DATASET_IDS_IN_MANIFEST_ORDER:
            table = pd.read_csv(self.output_dir / dataset_id / "candidate_snps.tsv", sep="\t")
            self.assertTrue((table["dataset_id"] == dataset_id).all())

    def test_each_dataset_gets_its_own_two_candidates_one_signal(self) -> None:
        # Each dataset's 2 tiny p-values (pos1, pos2, 500bp apart) cluster into 1 signal
        # under clustering_distance=1000; the 3rd row (p=0.9) is never significant.
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        for result in outcome.datasets:
            self.assertEqual(result.n_signals, 1)
            self.assertEqual(result.n_candidates, 2)

    def test_root_batch_summary_still_present_alongside_dataset_directories(self) -> None:
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        self.assertTrue((self.output_dir / "batch_summary.tsv").is_file())


class InvalidClusteringDistanceTests(BatchCandidatesTestCase):
    def test_negative_clustering_distance_fails_before_anything_is_created(self) -> None:
        self._write_all_six_valid()
        with self.assertRaises(ValueError):
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
                clustering_distance=-1,
            )
        self.assertFalse(self.output_dir.exists())


if __name__ == "__main__":
    unittest.main()
