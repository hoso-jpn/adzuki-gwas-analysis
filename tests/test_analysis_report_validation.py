"""Unit tests for adzuki_gwas_analysis.analysis.report_validation.validate_analysis_dir.

Every test uses a real, small candidate-enabled ``batch`` output built by
:func:`tests.analysis_support.build_candidate_enabled_batch_fixture` (which itself calls the
real ``run_batch``) -- never a hand-authored ``batch_summary.tsv`` that could silently drift
from the real schema. Tampering scenarios load the real fixture, mutate one file with
pandas, and confirm the specific violation is caught.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from adzuki_gwas_analysis.analysis import batch
from adzuki_gwas_analysis.analysis.report_validation import validate_analysis_dir
from adzuki_gwas_analysis.errors import (
    ReportRequiresCandidateEnabledBatchError,
    ReportSourceInconsistentError,
    ReportSourcePathUnsafeError,
)
from tests.analysis_support import (
    build_candidate_enabled_batch_fixture,
    write_full_manifest,
    write_six_dataset_files,
)


class ReportValidationTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.analysis_dir = build_candidate_enabled_batch_fixture(self.tmp_path)

    def _read_summary(self) -> pd.DataFrame:
        return pd.read_csv(self.analysis_dir / "batch_summary.tsv", sep="\t")

    def _write_summary(self, table: pd.DataFrame) -> None:
        table.to_csv(self.analysis_dir / "batch_summary.tsv", sep="\t", index=False)


class HappyPathTests(ReportValidationTestCase):
    def test_validates_all_six_datasets(self) -> None:
        validated = validate_analysis_dir(self.analysis_dir)
        self.assertEqual(len(validated.datasets), 6)

    def test_zero_candidate_dataset_validates_cleanly(self) -> None:
        validated = validate_analysis_dir(self.analysis_dir)
        zero = next(
            d for d in validated.datasets if d.dataset_id == "miyagi_mottled_black_seedcoat"
        )
        self.assertEqual(zero.n_candidates, 0)
        self.assertEqual(zero.n_signals, 0)
        self.assertEqual(len(zero.candidate_snps_df), 0)
        self.assertEqual(len(zero.candidate_ranking_df), 0)
        self.assertEqual(len(zero.association_peaks_df), 0)

    def test_global_parameters_are_uniform(self) -> None:
        validated = validate_analysis_dir(self.analysis_dir)
        self.assertEqual(validated.alpha, 0.05)
        self.assertEqual(validated.fdr_level, 0.05)
        self.assertEqual(validated.clustering_distance, 1000)


class RequiresCandidateEnabledBatchTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)

    def test_schema_v1_batch_output_is_rejected(self) -> None:
        data_dir = self.tmp_path / "data"
        data_dir.mkdir()
        manifest_path = self.tmp_path / "manifest.toml"
        analysis_dir = self.tmp_path / "no_candidates"
        rows = ["Chr01\t.\t1000000\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-3\t1e-3"]
        dataset_paths = write_six_dataset_files(
            data_dir,
            {
                did: rows
                for did in (
                    "miyagi_water_permeability",
                    "miyagi_red_seedcoat",
                    "miyagi_mottled_black_seedcoat",
                    "shumari_water_permeability",
                    "shumari_red_seedcoat",
                    "shumari_mottled_black_seedcoat",
                )
            },
        )
        write_full_manifest(manifest_path, dataset_paths)
        batch.run_batch(manifest_path=manifest_path, data_dir=data_dir, output_dir=analysis_dir)

        with self.assertRaises(ReportRequiresCandidateEnabledBatchError):
            validate_analysis_dir(analysis_dir)


class MissingArtifactTests(ReportValidationTestCase):
    def test_missing_batch_summary_raises_path_unsafe(self) -> None:
        (self.analysis_dir / "batch_summary.tsv").unlink()
        with self.assertRaises(ReportSourcePathUnsafeError):
            validate_analysis_dir(self.analysis_dir)

    def test_missing_referenced_artifact_file_raises_path_unsafe(self) -> None:
        target = self.analysis_dir / "miyagi_water_permeability" / "candidate_ranking.tsv"
        target.unlink()
        with self.assertRaises(ReportSourcePathUnsafeError):
            validate_analysis_dir(self.analysis_dir)


class PathSafetyTests(ReportValidationTestCase):
    def test_absolute_artifact_path_is_rejected(self) -> None:
        summary = self._read_summary()
        summary.loc[0, "candidate_ranking_path"] = "/etc/passwd"
        self._write_summary(summary)
        with self.assertRaises(ReportSourcePathUnsafeError):
            validate_analysis_dir(self.analysis_dir)

    def test_parent_traversal_artifact_path_is_rejected(self) -> None:
        summary = self._read_summary()
        summary.loc[0, "candidate_ranking_path"] = "../outside.tsv"
        self._write_summary(summary)
        with self.assertRaises(ReportSourcePathUnsafeError):
            validate_analysis_dir(self.analysis_dir)

    def test_symlink_artifact_path_is_rejected(self) -> None:
        real_target = self.analysis_dir / "miyagi_water_permeability" / "candidate_ranking.tsv"
        symlink_path = self.analysis_dir / "miyagi_water_permeability" / "sneaky_link.tsv"
        symlink_path.symlink_to(real_target)
        summary = self._read_summary()
        summary.loc[0, "candidate_ranking_path"] = "miyagi_water_permeability/sneaky_link.tsv"
        self._write_summary(summary)
        with self.assertRaises(ReportSourcePathUnsafeError):
            validate_analysis_dir(self.analysis_dir)

    def test_broken_symlink_artifact_path_is_rejected(self) -> None:
        broken_link = self.analysis_dir / "miyagi_water_permeability" / "broken_link.tsv"
        broken_link.symlink_to(self.analysis_dir / "does_not_exist.tsv")
        summary = self._read_summary()
        summary.loc[0, "candidate_ranking_path"] = "miyagi_water_permeability/broken_link.tsv"
        self._write_summary(summary)
        with self.assertRaises(ReportSourcePathUnsafeError):
            validate_analysis_dir(self.analysis_dir)

    def test_symlinked_directory_escape_is_rejected(self) -> None:
        outside_dir = self.tmp_path / "outside"
        outside_dir.mkdir()
        (outside_dir / "candidate_ranking.tsv").write_text("not real data", encoding="utf-8")
        escape_link = self.analysis_dir / "escape_dir"
        escape_link.symlink_to(outside_dir)
        summary = self._read_summary()
        summary.loc[0, "candidate_ranking_path"] = "escape_dir/candidate_ranking.tsv"
        self._write_summary(summary)
        with self.assertRaises(ReportSourcePathUnsafeError):
            validate_analysis_dir(self.analysis_dir)


class ConsistencyTests(ReportValidationTestCase):
    def test_nonuniform_alpha_across_rows_is_rejected(self) -> None:
        summary = self._read_summary()
        summary.loc[0, "alpha"] = 0.01
        self._write_summary(summary)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)

    def test_n_tests_mismatch_between_summary_and_diagnostics_is_rejected(self) -> None:
        summary = self._read_summary()
        summary.loc[0, "n_tests"] = 999999
        self._write_summary(summary)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)

    def test_dataset_id_mismatch_in_diagnostics_is_rejected(self) -> None:
        diagnostics_path = (
            self.analysis_dir / "miyagi_water_permeability" / "statistical_diagnostics.tsv"
        )
        table = pd.read_csv(diagnostics_path, sep="\t")
        table.loc[0, "dataset_id"] = "wrong_dataset_id"
        table.to_csv(diagnostics_path, sep="\t", index=False)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)

    def test_n_signals_mismatch_with_association_peaks_row_count_is_rejected(self) -> None:
        peaks_path = self.analysis_dir / "miyagi_red_seedcoat" / "association_peaks.tsv"
        table = pd.read_csv(peaks_path, sep="\t")
        table = pd.concat([table, table.iloc[[0]]], ignore_index=True)
        table.to_csv(peaks_path, sep="\t", index=False)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)

    def test_wrong_dataset_id_in_candidate_snps_is_rejected(self) -> None:
        snps_path = self.analysis_dir / "miyagi_water_permeability" / "candidate_snps.tsv"
        table = pd.read_csv(snps_path, sep="\t")
        table.loc[0, "dataset_id"] = "miyagi_red_seedcoat"
        table.to_csv(snps_path, sep="\t", index=False)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)

    def test_reference_mismatch_in_candidate_ranking_is_rejected(self) -> None:
        ranking_path = self.analysis_dir / "miyagi_water_permeability" / "candidate_ranking.tsv"
        table = pd.read_csv(ranking_path, sep="\t")
        table.loc[0, "reference"] = "Shumari"
        table.to_csv(ranking_path, sep="\t", index=False)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)

    def test_trait_mismatch_in_association_peaks_is_rejected(self) -> None:
        peaks_path = self.analysis_dir / "miyagi_water_permeability" / "association_peaks.tsv"
        table = pd.read_csv(peaks_path, sep="\t")
        table.loc[0, "trait"] = "red_seedcoat"
        table.to_csv(peaks_path, sep="\t", index=False)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)

    def test_non_dense_candidate_rank_sequence_is_rejected(self) -> None:
        ranking_path = self.analysis_dir / "shumari_red_seedcoat" / "candidate_ranking.tsv"
        table = pd.read_csv(ranking_path, sep="\t")
        table.loc[0, "candidate_rank"] = table["candidate_rank"].max() + 5
        table.to_csv(ranking_path, sep="\t", index=False)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)

    def test_unknown_signal_id_reference_is_rejected(self) -> None:
        ranking_path = self.analysis_dir / "shumari_red_seedcoat" / "candidate_ranking.tsv"
        table = pd.read_csv(ranking_path, sep="\t")
        table.loc[0, "signal_id"] = "shumari_red_seedcoat_peak_9999"
        table.to_csv(ranking_path, sep="\t", index=False)
        with self.assertRaises(ReportSourceInconsistentError):
            validate_analysis_dir(self.analysis_dir)


if __name__ == "__main__":
    unittest.main()
