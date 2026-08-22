"""Failure/staging/transaction tests for adzuki_gwas_analysis.analysis.batch.run_batch.

Every scenario here must leave behind no partial batch output: --output-dir must not
exist (or, if it pre-existed empty, must remain empty) and the staging directory created
under output_dir.parent must be removed.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

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
    rows_by_id: dict[str, list[str]] = {}
    for i, dataset_id in enumerate(_DATASET_IDS_IN_MANIFEST_ORDER):
        rows_by_id[dataset_id] = [
            row.format(
                pos1=1_000_000 + i,
                pos2=2_000_000 + i,
                pos3=1_500_000 + i,
                p1="2e-2",
                p2="4e-2",
                p3="8e-3",
            )
            for row in _ROW_TEMPLATE
        ]
    return rows_by_id


class BatchFailureTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.data_dir = self.tmp_path / "data"
        self.data_dir.mkdir()
        self.manifest_path = self.tmp_path / "manifest.toml"
        self.output_dir = self.tmp_path / "out"
        self._recorded_staging_dirs: list[Path] = []

    def _write_all_six_valid(self) -> None:
        dataset_paths = write_six_dataset_files(self.data_dir, _six_dataset_rows())
        write_full_manifest(self.manifest_path, dataset_paths)

    def _run_batch_recording_staging_dir(self, **kwargs: object) -> None:
        original_mkdtemp = batch.tempfile.mkdtemp

        def _recording_mkdtemp(*args: object, **inner_kwargs: object) -> str:
            path = original_mkdtemp(*args, **inner_kwargs)
            self._recorded_staging_dirs.append(Path(path))
            return path

        with mock.patch(
            "adzuki_gwas_analysis.analysis.batch.tempfile.mkdtemp",
            side_effect=_recording_mkdtemp,
        ):
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
                **kwargs,  # type: ignore[arg-type]
            )

    def _assert_no_partial_output(self) -> None:
        self.assertFalse(self.output_dir.exists(), "output_dir must not exist after a failure")
        self.assertEqual(len(self._recorded_staging_dirs), 1)
        self.assertFalse(
            self._recorded_staging_dirs[0].exists(),
            "the staging directory must be removed after a failure",
        )


class ValidationFailureTests(BatchFailureTestCase):
    def test_second_dataset_validation_failure_leaves_no_partial_output(self) -> None:
        self._write_all_six_valid()
        # Corrupt the 2nd dataset's file after the manifest already recorded its real
        # checksum, so validate_dataset() fails specifically for miyagi_red_seedcoat.
        second_dataset_path = self.data_dir / "mapped_to_Miyagi_red_seedcoat.maf_0.05.assoc.txt"
        with second_dataset_path.open("a", encoding="utf-8") as fh:
            fh.write(
                "Chr01\t.\t9999999\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t1e-3\t1e-3\n"
            )

        with self.assertRaises(Exception):  # noqa: B017
            self._run_batch_recording_staging_dir()
        self._assert_no_partial_output()


class LoadFailureTests(BatchFailureTestCase):
    def test_fourth_dataset_load_failure_leaves_no_partial_output(self) -> None:
        self._write_all_six_valid()
        call_count = 0
        original_loader = batch.load_analysis_frame

        def _failing_on_fourth(path: Path):  # type: ignore[no-untyped-def]
            nonlocal call_count
            call_count += 1
            if call_count == 4:
                raise OSError("simulated load failure for the 4th dataset")
            return original_loader(path)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.batch.load_analysis_frame",
                side_effect=_failing_on_fourth,
            ),
            self.assertRaises(OSError),
        ):
            self._run_batch_recording_staging_dir()
        self._assert_no_partial_output()
        self.assertEqual(call_count, 4)


class PlottingFailureTests(BatchFailureTestCase):
    def test_manhattan_plot_failure_midway_leaves_no_partial_output(self) -> None:
        self._write_all_six_valid()
        call_count = 0
        original = batch.plot_manhattan

        def _failing_on_third(*args, **kwargs):  # type: ignore[no-untyped-def]
            nonlocal call_count
            call_count += 1
            if call_count == 3:
                raise RuntimeError("simulated Manhattan plot failure")
            return original(*args, **kwargs)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.batch.plot_manhattan", side_effect=_failing_on_third
            ),
            self.assertRaises(RuntimeError),
        ):
            self._run_batch_recording_staging_dir()
        self._assert_no_partial_output()

    def test_qq_plot_failure_midway_leaves_no_partial_output(self) -> None:
        self._write_all_six_valid()
        call_count = 0
        original = batch.plot_qq

        def _failing_on_third(*args, **kwargs):  # type: ignore[no-untyped-def]
            nonlocal call_count
            call_count += 1
            if call_count == 3:
                raise RuntimeError("simulated QQ plot failure")
            return original(*args, **kwargs)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.batch.plot_qq", side_effect=_failing_on_third
            ),
            self.assertRaises(RuntimeError),
        ):
            self._run_batch_recording_staging_dir()
        self._assert_no_partial_output()


class TsvWriteFailureTests(BatchFailureTestCase):
    def test_per_dataset_diagnostics_tsv_write_failure_leaves_no_partial_output(self) -> None:
        self._write_all_six_valid()
        call_count = 0
        original = batch.atomic_write_tsv

        def _failing_on_third_call(table, output_path):  # type: ignore[no-untyped-def]
            nonlocal call_count
            call_count += 1
            # Each dataset writes 2 TSVs (diagnostics, significant_variants); the 3rd call
            # overall is the 2nd dataset's diagnostics TSV.
            if call_count == 3:
                raise OSError("simulated diagnostics TSV write failure")
            return original(table, output_path)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.batch.atomic_write_tsv",
                side_effect=_failing_on_third_call,
            ),
            self.assertRaises(OSError),
        ):
            self._run_batch_recording_staging_dir()
        self._assert_no_partial_output()

    def test_batch_summary_write_failure_leaves_no_partial_output(self) -> None:
        self._write_all_six_valid()
        call_count = 0
        original = batch.atomic_write_tsv

        def _failing_on_last_call(table, output_path):  # type: ignore[no-untyped-def]
            nonlocal call_count
            call_count += 1
            # 6 datasets x 2 TSVs each = 12 calls, then the 13th is batch_summary.tsv.
            if call_count == 13:
                raise OSError("simulated batch_summary.tsv write failure")
            return original(table, output_path)

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.batch.atomic_write_tsv",
                side_effect=_failing_on_last_call,
            ),
            self.assertRaises(OSError),
        ):
            self._run_batch_recording_staging_dir()
        self._assert_no_partial_output()
        self.assertEqual(call_count, 13)


class PublishFailureTests(BatchFailureTestCase):
    def test_final_publish_failure_leaves_no_partial_output(self) -> None:
        # os.replace(staging_dir, output_dir) is the very last step of run_batch, after
        # every dataset and batch_summary.tsv already succeeded. If it fails (a race, a
        # permission error, a cross-device rename), the staging directory -- which by this
        # point holds a fully-built, 25-file tree -- must still be cleaned up, not left
        # behind as an orphaned near-complete batch.
        self._write_all_six_valid()
        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.batch.os.replace",
                side_effect=OSError("simulated publish failure"),
            ),
            self.assertRaises(OSError),
        ):
            self._run_batch_recording_staging_dir()
        self._assert_no_partial_output()


class InvalidParameterTests(BatchFailureTestCase):
    def test_invalid_alpha_fails_before_anything_is_created(self) -> None:
        self._write_all_six_valid()
        with self.assertRaises(ValueError):
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
                alpha=0.0,
            )
        self.assertFalse(self.output_dir.exists())

    def test_invalid_fdr_level_fails_before_anything_is_created(self) -> None:
        self._write_all_six_valid()
        with self.assertRaises(ValueError):
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
                fdr_level=1.5,
            )
        self.assertFalse(self.output_dir.exists())

    def test_invalid_threshold_fails_before_anything_is_created(self) -> None:
        self._write_all_six_valid()
        with self.assertRaises(ValueError):
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
                threshold=-1.0,
            )
        self.assertFalse(self.output_dir.exists())

    def test_invalid_alpha_never_creates_a_staging_directory(self) -> None:
        self._write_all_six_valid()
        with (
            mock.patch("adzuki_gwas_analysis.analysis.batch.tempfile.mkdtemp") as mocked_mkdtemp,
            self.assertRaises(ValueError),
        ):
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
                alpha=2.0,
            )
        mocked_mkdtemp.assert_not_called()


if __name__ == "__main__":
    unittest.main()
