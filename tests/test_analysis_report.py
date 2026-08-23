"""Unit tests for adzuki_gwas_analysis.analysis.report.run_report (happy-path behavior).

Every test uses a real, small candidate-enabled ``batch`` output built by
:func:`tests.analysis_support.build_candidate_enabled_batch_fixture`. Failure/atomicity
scenarios live in tests/test_analysis_report_failures.py.
"""

from __future__ import annotations

import getpass
import json
import os
import socket
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.analysis.report import run_report
from adzuki_gwas_analysis.loader import compute_sha256
from tests.analysis_support import build_candidate_enabled_batch_fixture

_ALLOWLISTED_PER_DATASET_SUFFIXES = (
    "_manhattan.png",
    "_qq.png",
    "statistical_diagnostics.tsv",
    "significant_variants.tsv",
    "association_peaks.tsv",
    "candidate_snps.tsv",
    "candidate_ranking.tsv",
)


class ReportTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.analysis_dir = build_candidate_enabled_batch_fixture(self.tmp_path)
        self.output_dir = self.tmp_path / "delivery"

    def _all_files(self, root: Path) -> list[Path]:
        return [p for p in root.rglob("*") if p.is_file()]


class OutputStructureTests(ReportTestCase):
    def test_top_level_files_and_directories(self) -> None:
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        top_level = {p.name for p in self.output_dir.iterdir()}
        self.assertEqual(
            top_level,
            {"executive_summary.md", "analysis_report.md", "artifacts", "reproducibility"},
        )

    def test_reproducibility_has_exactly_three_files(self) -> None:
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        files = {p.name for p in (self.output_dir / "reproducibility").iterdir()}
        self.assertEqual(
            files, {"input_checksums.tsv", "software_versions.json", "run_manifest.json"}
        )

    def test_artifacts_has_batch_summary_and_six_dataset_directories(self) -> None:
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        artifacts_dir = self.output_dir / "artifacts"
        self.assertTrue((artifacts_dir / "batch_summary.tsv").is_file())
        dataset_dirs = {p.name for p in artifacts_dir.iterdir() if p.is_dir()}
        self.assertEqual(len(dataset_dirs), 6)

    def test_no_raw_assoc_txt_file_is_ever_copied(self) -> None:
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        for path in self._all_files(self.output_dir):
            self.assertFalse(path.name.endswith(".assoc.txt"), path)

    def test_only_allowlisted_filenames_appear_under_a_dataset_directory(self) -> None:
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        artifacts_dir = self.output_dir / "artifacts"
        for dataset_dir in (p for p in artifacts_dir.iterdir() if p.is_dir()):
            for file_path in dataset_dir.iterdir():
                self.assertTrue(
                    any(
                        file_path.name.endswith(suffix)
                        for suffix in _ALLOWLISTED_PER_DATASET_SUFFIXES
                    ),
                    file_path,
                )


class ChecksumIntegrityTests(ReportTestCase):
    def test_every_manifest_artifact_checksum_matches_the_actual_delivered_file(self) -> None:
        outcome = run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        manifest = json.loads(outcome.run_manifest_path.read_text(encoding="utf-8"))
        for entry in manifest["artifacts"]:
            delivered = self.output_dir / entry["relative_path"]
            self.assertTrue(delivered.is_file(), entry["relative_path"])
            self.assertEqual(compute_sha256(delivered), entry["sha256"], entry["relative_path"])

    def test_run_manifest_itself_is_not_in_its_own_artifact_inventory(self) -> None:
        outcome = run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        manifest = json.loads(outcome.run_manifest_path.read_text(encoding="utf-8"))
        relative_paths = {entry["relative_path"] for entry in manifest["artifacts"]}
        self.assertNotIn("reproducibility/run_manifest.json", relative_paths)

    def test_copied_artifact_content_matches_source_analysis_dir_content(self) -> None:
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        src = self.analysis_dir / "batch_summary.tsv"
        dst = self.output_dir / "artifacts" / "batch_summary.tsv"
        self.assertEqual(src.read_bytes(), dst.read_bytes())


class ReadOnlySourceTests(ReportTestCase):
    def test_analysis_dir_is_byte_for_byte_unchanged_after_report_generation(self) -> None:
        before = {path: compute_sha256(path) for path in self._all_files(self.analysis_dir)}
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        after = {path: compute_sha256(path) for path in self._all_files(self.analysis_dir)}
        self.assertEqual(before, after)

    def test_analysis_dir_file_count_is_unchanged(self) -> None:
        before_count = len(self._all_files(self.analysis_dir))
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        after_count = len(self._all_files(self.analysis_dir))
        self.assertEqual(before_count, after_count)


class ConfidentialityTests(ReportTestCase):
    def test_no_hostname_username_or_home_path_anywhere_in_the_delivery_package(self) -> None:
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        hostname = socket.gethostname()
        username = getpass.getuser()
        home = os.path.expanduser("~")
        text_extensions = {".md", ".json", ".tsv"}
        for path in self._all_files(self.output_dir):
            if path.suffix not in text_extensions:
                continue
            content = path.read_text(encoding="utf-8", errors="replace")
            self.assertNotIn(hostname, content, path)
            if username and len(username) > 2:
                self.assertNotIn(username, content, path)
            self.assertNotIn(home, content, path)
            self.assertNotIn(str(self.tmp_path), content, path)

    def test_software_versions_has_no_environ_dump(self) -> None:
        run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        doc = json.loads(
            (self.output_dir / "reproducibility" / "software_versions.json").read_text(
                encoding="utf-8"
            )
        )
        env = doc["report_generation_environment"]
        self.assertEqual(set(env.keys()), {"python_version", "platform", "packages", "git_commit"})

    def test_run_manifest_has_no_absolute_paths_in_artifact_inventory(self) -> None:
        outcome = run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        manifest = json.loads(outcome.run_manifest_path.read_text(encoding="utf-8"))
        for entry in manifest["artifacts"]:
            self.assertFalse(Path(entry["relative_path"]).is_absolute(), entry["relative_path"])


class DeterminismTests(ReportTestCase):
    def test_two_independent_runs_produce_identical_markdown_and_manifest_content(self) -> None:
        output_dir_a = self.tmp_path / "delivery_a"
        output_dir_b = self.tmp_path / "delivery_b"
        run_report(analysis_dir=self.analysis_dir, output_dir=output_dir_a)
        run_report(analysis_dir=self.analysis_dir, output_dir=output_dir_b)
        self.assertEqual(
            (output_dir_a / "executive_summary.md").read_text(),
            (output_dir_b / "executive_summary.md").read_text(),
        )
        self.assertEqual(
            (output_dir_a / "analysis_report.md").read_text(),
            (output_dir_b / "analysis_report.md").read_text(),
        )
        manifest_a = json.loads(
            (output_dir_a / "reproducibility" / "run_manifest.json").read_text()
        )
        manifest_b = json.loads(
            (output_dir_b / "reproducibility" / "run_manifest.json").read_text()
        )
        self.assertEqual(manifest_a, manifest_b)


class OutputDirSafetyTests(ReportTestCase):
    def test_existing_nonempty_output_dir_is_rejected(self) -> None:
        self.output_dir.mkdir(parents=True)
        (self.output_dir / "stray.txt").write_text("x", encoding="utf-8")
        with self.assertRaises(Exception):  # noqa: B017
            run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)
        self.assertEqual([p.name for p in self.output_dir.iterdir()], ["stray.txt"])

    def test_symlink_output_dir_is_rejected(self) -> None:
        real_dir = self.tmp_path / "real_target"
        real_dir.mkdir()
        self.output_dir.symlink_to(real_dir)
        with self.assertRaises(Exception):  # noqa: B017
            run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)


if __name__ == "__main__":
    unittest.main()
