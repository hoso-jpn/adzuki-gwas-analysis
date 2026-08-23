"""Failure/staging/transaction tests for adzuki_gwas_analysis.analysis.report.run_report.

Every scenario here must leave behind no partial delivery package: ``--output-dir`` must
not exist (or, if it pre-existed empty, must remain empty) and the staging directory
created under ``output_dir.parent`` must be removed.

Failure injection targets one specific call by inspecting its arguments (which file it is
about to touch), never a global function unconditionally -- an earlier version of a
similar test in this repository (PR #17's batch failure-path test) patched ``os.replace``
unconditionally and the injected failure fired on the very first unrelated write instead
of the step it claimed to test. See ``PublishFailureTests`` below for the exact,
destination-checked pattern that avoids that mistake, reused here for every injected
failure.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from adzuki_gwas_analysis.analysis import report
from tests.analysis_support import build_candidate_enabled_batch_fixture

#: The full delivery package has exactly this many files once every step succeeds --
#: 1 batch_summary.tsv + 6 datasets x 7 files + 2 markdown reports + 3 reproducibility
#: files. Verified directly against a real run in tests/test_analysis_report.py; used here
#: only to confirm a failure fires at the intended step, not to duplicate that coverage.
_FULL_PACKAGE_FILE_COUNT = 48


class ReportFailureTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.analysis_dir = build_candidate_enabled_batch_fixture(self.tmp_path)
        self.output_dir = self.tmp_path / "delivery"
        self._recorded_staging_dirs: list[Path] = []

    def _run_report_recording_staging_dir(self) -> None:
        original_mkdtemp = report.tempfile.mkdtemp

        def _recording_mkdtemp(*args: object, **kwargs: object) -> str:
            path = original_mkdtemp(*args, **kwargs)
            self._recorded_staging_dirs.append(Path(path))
            return path

        with mock.patch(
            "adzuki_gwas_analysis.analysis.report.tempfile.mkdtemp",
            side_effect=_recording_mkdtemp,
        ):
            report.run_report(analysis_dir=self.analysis_dir, output_dir=self.output_dir)

    def _assert_no_partial_output(self) -> None:
        self.assertFalse(self.output_dir.exists(), "output_dir must not exist after a failure")
        self.assertEqual(len(self._recorded_staging_dirs), 1)
        self.assertFalse(
            self._recorded_staging_dirs[0].exists(),
            "the staging directory must be removed after a failure",
        )


class MarkdownWriteFailureTests(ReportFailureTestCase):
    def test_executive_summary_write_failure_leaves_no_partial_output(self) -> None:
        original = report._atomic_write_text

        def _fail_only_for_executive_summary(content: str, output_path: Path) -> None:
            if output_path.name == "executive_summary.md":
                raise OSError("simulated executive_summary.md write failure")
            original(content, output_path)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.report._atomic_write_text",
                side_effect=_fail_only_for_executive_summary,
            ),
            self.assertRaises(OSError),
        ):
            self._run_report_recording_staging_dir()
        self._assert_no_partial_output()


class ArtifactCopyFailureTests(ReportFailureTestCase):
    def test_one_dataset_manhattan_png_copy_failure_leaves_no_partial_output(self) -> None:
        original = report.shutil.copyfile
        target_name = "shumari_red_seedcoat_manhattan.png"

        def _fail_only_for_target(src: Path, dst: Path) -> None:
            if Path(dst).name == target_name:
                raise OSError("simulated artifact copy failure")
            original(src, dst)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.report.shutil.copyfile",
                side_effect=_fail_only_for_target,
            ),
            self.assertRaises(OSError),
        ):
            self._run_report_recording_staging_dir()
        self._assert_no_partial_output()


class ChecksumFailureTests(ReportFailureTestCase):
    def test_checksum_failure_on_analysis_report_md_leaves_no_partial_output(self) -> None:
        original = report.compute_sha256
        call_count = 0

        def _fail_on_analysis_report_checksum(path: Path, **kwargs: object) -> str:
            nonlocal call_count
            call_count += 1
            if Path(path).name == "analysis_report.md":
                raise OSError("simulated checksum failure")
            return original(path, **kwargs)  # type: ignore[arg-type]

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.report.compute_sha256",
                side_effect=_fail_on_analysis_report_checksum,
            ),
            self.assertRaises(OSError),
        ):
            self._run_report_recording_staging_dir()
        self._assert_no_partial_output()
        self.assertGreater(call_count, 0)


class ManifestWriteFailureTests(ReportFailureTestCase):
    def test_run_manifest_json_write_failure_leaves_no_partial_output(self) -> None:
        # run_manifest.json is written last, after every other file has already succeeded
        # -- the staging tree holds 47 files (everything except run_manifest.json itself)
        # at the moment this failure fires.
        original = report._atomic_write_text
        seen_staging_dir: list[Path] = []

        def _fail_only_for_run_manifest(content: str, output_path: Path) -> None:
            if output_path.name == "run_manifest.json":
                staging_dir = output_path.parent.parent
                seen_staging_dir.append(staging_dir)
                generated_files = [p for p in staging_dir.rglob("*") if p.is_file()]
                self.assertEqual(len(generated_files), _FULL_PACKAGE_FILE_COUNT - 1)
                raise OSError("simulated run_manifest.json write failure")
            original(content, output_path)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.report._atomic_write_text",
                side_effect=_fail_only_for_run_manifest,
            ),
            self.assertRaises(OSError),
        ):
            self._run_report_recording_staging_dir()
        self.assertEqual(len(seen_staging_dir), 1)
        self._assert_no_partial_output()


class PublishFailureTests(ReportFailureTestCase):
    def test_final_publish_failure_leaves_no_partial_output(self) -> None:
        # os.replace(staging_dir, output_dir) is the very last step of run_report, after
        # every artifact copy, both Markdown reports, and all 3 reproducibility files have
        # already succeeded. os.replace is one process-wide function -- _atomic_write_text
        # calls it too, once per file, for its own tmp-then-rename writes. An unconditional
        # mock would therefore fire on the very first Markdown/JSON write instead of the
        # final publish (exactly the false-positive PR #17 fixed for batch.py's own
        # equivalent test). This version only fails the rename whose destination is
        # output_dir itself, letting every other os.replace call run for real, and asserts
        # the staging tree already holds the complete, correct file count at that point.
        original_replace = report.os.replace
        publish_attempts = 0

        def _fail_only_on_final_publish(src: Path, dst: Path) -> None:
            nonlocal publish_attempts
            if Path(dst) == self.output_dir:
                publish_attempts += 1
                generated_files = [p for p in Path(src).rglob("*") if p.is_file()]
                self.assertEqual(len(generated_files), _FULL_PACKAGE_FILE_COUNT)
                raise OSError("simulated final publish failure")
            original_replace(src, dst)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.report.os.replace",
                side_effect=_fail_only_on_final_publish,
            ),
            self.assertRaises(OSError),
        ):
            self._run_report_recording_staging_dir()

        self.assertEqual(publish_attempts, 1)
        self._assert_no_partial_output()


class InvalidAnalysisDirTests(ReportFailureTestCase):
    def test_missing_analysis_dir_never_creates_a_staging_directory(self) -> None:
        with (
            mock.patch("adzuki_gwas_analysis.analysis.report.tempfile.mkdtemp") as mocked,
            self.assertRaises(Exception),  # noqa: B017
        ):
            report.run_report(
                analysis_dir=self.tmp_path / "does_not_exist", output_dir=self.output_dir
            )
        mocked.assert_not_called()
        self.assertFalse(self.output_dir.exists())

    def test_schema_v1_analysis_dir_never_creates_a_staging_directory(self) -> None:
        from adzuki_gwas_analysis.analysis import batch
        from tests.analysis_support import write_full_manifest, write_six_dataset_files

        data_dir = self.tmp_path / "no_candidates_data"
        data_dir.mkdir()
        manifest_path = self.tmp_path / "no_candidates_manifest.toml"
        no_candidates_analysis_dir = self.tmp_path / "no_candidates_analysis_dir"
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
        batch.run_batch(
            manifest_path=manifest_path, data_dir=data_dir, output_dir=no_candidates_analysis_dir
        )

        with (
            mock.patch("adzuki_gwas_analysis.analysis.report.tempfile.mkdtemp") as mocked,
            self.assertRaises(Exception),  # noqa: B017
        ):
            report.run_report(analysis_dir=no_candidates_analysis_dir, output_dir=self.output_dir)
        mocked.assert_not_called()
        self.assertFalse(self.output_dir.exists())


if __name__ == "__main__":
    unittest.main()
