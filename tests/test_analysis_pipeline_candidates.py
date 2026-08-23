"""Unit tests for adzuki_gwas_analysis.analysis.pipeline.run_candidates.

A separate file from tests/test_analysis_pipeline_diagnostics.py (left unmodified by this
Issue) so that module's diff stays empty. Every test uses tiny synthetic fixtures under a
TemporaryDirectory -- none reads data/raw/ (the real, un-tracked Dryad file).
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from adzuki_gwas_analysis.analysis import pipeline
from adzuki_gwas_analysis.errors import DatasetValidationFailedError
from tests.analysis_support import write_dataset_file, write_manifest

# pval column (14th field): 1e-4, 4e-2, 2e-4, 8e-3. With alpha=0.05/m=4, the Bonferroni
# threshold is 0.0125: row 2 (pval=4e-2) fails it but every row passes BH at fdr_level=0.05
# (all 4 discovered) -- giving one Bonferroni-and-BH candidate row and one BH-only row,
# useful for exercising both priority tiers from real orchestration output.
_VALID_ROWS = [
    "Chr01\t.\t1000000\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-4\t1e-5",
    "Chr01\t.\t2000000\t0\tC\tT\t0.45\t-0.02\t0.02\t71.0\t31.0\t29.0\t5e-2\t4e-2\t3e-2",
    "Chr01\t.\t6000000\t0\tG\tA\t0.20\t0.03\t0.01\t70.5\t30.5\t28.5\t2e-3\t2e-4\t2e-5",
    "Chr02\t.\t1500000\t1\tT\tA\t0.10\t0.11\t0.03\t72.0\t32.0\t30.0\t9e-3\t8e-3\t7e-3",
]
_INVALID_ROWS = _VALID_ROWS + [
    "Chr02\t.\t9000000\t0\tG\tC\t0.5\t0.01\t0.01\t70.0\t30.0\t28.0\t1.5\t1.5\t1.5",
]


class CandidatesPipelineTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.data_dir = self.tmp_path / "data"
        self.data_dir.mkdir()
        self.manifest_path = self.tmp_path / "manifest.toml"
        self.output_dir = self.tmp_path / "out"

    def _write_valid_dataset(self) -> None:
        dataset_path = write_dataset_file(self.data_dir, _VALID_ROWS)
        write_manifest(self.manifest_path, dataset_path)

    def _write_invalid_dataset(self) -> None:
        dataset_path = write_dataset_file(self.data_dir, _INVALID_ROWS)
        write_manifest(self.manifest_path, dataset_path)


class RunCandidatesTests(CandidatesPipelineTestCase):
    def test_produces_exactly_five_files(self) -> None:
        self._write_valid_dataset()
        pipeline.run_candidates(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        produced = {p.name for p in self.output_dir.iterdir()}
        self.assertEqual(
            produced,
            {
                "statistical_diagnostics.tsv",
                "significant_variants.tsv",
                "association_peaks.tsv",
                "candidate_snps.tsv",
                "candidate_ranking.tsv",
            },
        )

    def test_bonferroni_and_bh_only_candidates_both_present(self) -> None:
        self._write_valid_dataset()
        outcome = pipeline.run_candidates(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        ranking = outcome.candidates_result.candidate_ranking
        self.assertEqual(set(ranking["priority_tier"]), {1, 2})
        self.assertEqual(outcome.candidates_result.n_candidates, 4)

    def test_dataset_context_is_correct_on_every_output_row(self) -> None:
        self._write_valid_dataset()
        pipeline.run_candidates(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
            clustering_distance=1000,
        )
        for filename in ("association_peaks.tsv", "candidate_snps.tsv", "candidate_ranking.tsv"):
            table = pd.read_csv(self.output_dir / filename, sep="\t")
            self.assertTrue((table["dataset_id"] == "miyagi_water_permeability").all())
            self.assertTrue((table["reference"] == "Miyagi").all())
            self.assertTrue((table["trait"] == "water_permeability").all())

    def test_never_clusters_across_the_two_chromosomes_present(self) -> None:
        self._write_valid_dataset()
        outcome = pipeline.run_candidates(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
            clustering_distance=10_000_000,
        )
        peaks = outcome.candidates_result.association_peaks
        self.assertEqual(set(peaks["chromosome"]), {"Chr01", "Chr02"})
        # Even with a huge clustering_distance, Chr01's 3 far-apart significant variants
        # here (1e6, 2e6, 6e6) are all within 10e6bp of each other and merge into one
        # Chr01 signal, but Chr02's variant never joins it.
        self.assertEqual(len(peaks), 2)

    def test_produces_no_output_when_validation_fails(self) -> None:
        self._write_invalid_dataset()
        with self.assertRaises(DatasetValidationFailedError):
            pipeline.run_candidates(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_dir=self.output_dir,
                clustering_distance=1000,
            )
        self.assertFalse(self.output_dir.exists())

    def test_rejects_negative_clustering_distance(self) -> None:
        self._write_valid_dataset()
        with self.assertRaises(ValueError):
            pipeline.run_candidates(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                dataset_id="miyagi_water_permeability",
                output_dir=self.output_dir,
                clustering_distance=-1,
            )

    def test_alpha_and_fdr_level_are_forwarded_not_hardcoded(self) -> None:
        self._write_valid_dataset()
        # A very strict alpha (well below the 4-test Bonferroni threshold for every p-value
        # here) must drop bonferroni_significant to zero, cutting priority_tier 1 to nothing.
        outcome = pipeline.run_candidates(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=self.output_dir,
            clustering_distance=1000,
            alpha=1e-10,
            fdr_level=1e-10,
        )
        self.assertEqual(outcome.candidates_result.n_candidates, 0)

    def test_output_parent_directory_created(self) -> None:
        self._write_valid_dataset()
        nested = self.output_dir / "nested" / "deeper"
        pipeline.run_candidates(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            dataset_id="miyagi_water_permeability",
            output_dir=nested,
            clustering_distance=1000,
        )
        self.assertTrue((nested / "candidate_ranking.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
