"""Unit tests for adzuki_gwas_analysis.analysis.batch.run_batch (happy-path behavior).

Every test uses tiny synthetic fixtures for all 6 canonical datasets under a
TemporaryDirectory -- none reads data/raw/ (the real, un-tracked Dryad files). Failure/
staging/transaction scenarios live in tests/test_analysis_batch_failures.py.
"""

from __future__ import annotations

import dataclasses
import gc
import tempfile
import unittest
import weakref
from pathlib import Path
from unittest import mock

import pandas as pd

from adzuki_gwas_analysis.analysis import batch
from adzuki_gwas_analysis.manifest import load_manifest
from tests.analysis_support import write_full_manifest, write_six_dataset_files

_DATASET_IDS_IN_MANIFEST_ORDER: tuple[str, ...] = (
    "miyagi_water_permeability",
    "miyagi_red_seedcoat",
    "miyagi_mottled_black_seedcoat",
    "shumari_water_permeability",
    "shumari_red_seedcoat",
    "shumari_mottled_black_seedcoat",
)

_EXPECTED_MANHATTAN_TITLES = {
    "miyagi_water_permeability": "Water Permeability GWAS (Miyagi reference)",
    "miyagi_red_seedcoat": "Red Seed Coat Color GWAS (Miyagi reference)",
    "miyagi_mottled_black_seedcoat": "Mottled Black Seed Coat Color GWAS (Miyagi reference)",
    "shumari_water_permeability": "Water Permeability GWAS (Shumari reference)",
    "shumari_red_seedcoat": "Red Seed Coat Color GWAS (Shumari reference)",
    "shumari_mottled_black_seedcoat": "Mottled Black Seed Coat Color GWAS (Shumari reference)",
}

_EXPECTED_QQ_TITLES = {
    dataset_id: f"QQ Plot: {title}" for dataset_id, title in _EXPECTED_MANHATTAN_TITLES.items()
}

_ROW_TEMPLATE = (
    "Chr01\t.\t{pos1}\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t{p1}\t1e-3",
    "Chr01\t.\t{pos2}\t0\tC\tT\t0.45\t-0.02\t0.02\t71.0\t31.0\t29.0\t5e-2\t{p2}\t3e-2",
    "Chr02\t.\t{pos3}\t1\tT\tA\t0.10\t0.11\t0.03\t72.0\t32.0\t30.0\t9e-3\t{p3}\t7e-3",
)


def _six_dataset_rows(*, tiny_p_for_first_dataset: bool = True) -> dict[str, list[str]]:
    """3 valid, plottable, multi-chromosome rows per canonical dataset_id.

    Distinct ``pos`` values per dataset (offset by index) make it easy to assert that a
    given output file's content actually corresponds to the dataset directory it was
    written under. ``miyagi_water_permeability`` optionally gets one very small p-value
    so it has a nonzero, easily-distinguished Bonferroni/BH discovery count.
    """
    rows_by_id: dict[str, list[str]] = {}
    for i, dataset_id in enumerate(_DATASET_IDS_IN_MANIFEST_ORDER):
        p1 = "1e-9" if (tiny_p_for_first_dataset and i == 0) else "2e-2"
        rows_by_id[dataset_id] = [
            row.format(
                pos1=1_000_000 + i,
                pos2=2_000_000 + i,
                pos3=1_500_000 + i,
                p1=p1,
                p2="4e-2",
                p3="8e-3",
            )
            for row in _ROW_TEMPLATE
        ]
    return rows_by_id


class BatchTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        self.data_dir = self.tmp_path / "data"
        self.data_dir.mkdir()
        self.manifest_path = self.tmp_path / "manifest.toml"
        self.output_dir = self.tmp_path / "out"

    def _write_all_six_valid(self, **kwargs: bool) -> None:
        dataset_paths = write_six_dataset_files(self.data_dir, _six_dataset_rows(**kwargs))
        write_full_manifest(self.manifest_path, dataset_paths)


class ManifestOrderTests(BatchTestCase):
    def test_manifest_declares_datasets_in_expected_order(self) -> None:
        # Sanity check on the fixture helper itself: manifest.toml's own real dataset order
        # is Miyagi x 3 traits then Shumari x 3 traits (see manifest.py's
        # SCHEMA_V1_DATASETS), and the test fixture must reproduce that, not alphabetical
        # order (e.g. "miyagi_mottled..." would sort before "miyagi_red..." alphabetically,
        # which is NOT manifest.toml's actual order).
        self._write_all_six_valid()
        manifest = load_manifest(self.manifest_path)
        self.assertEqual(
            tuple(e.dataset_id for e in manifest.datasets), _DATASET_IDS_IN_MANIFEST_ORDER
        )

    def test_batch_processes_datasets_in_manifest_order(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        self.assertEqual(
            tuple(r.dataset_id for r in outcome.datasets), _DATASET_IDS_IN_MANIFEST_ORDER
        )

    def test_batch_summary_row_order_matches_manifest_order_not_alphabetical(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        summary = pd.read_csv(outcome.summary_path, sep="\t")
        self.assertEqual(list(summary["dataset_id"]), list(_DATASET_IDS_IN_MANIFEST_ORDER))
        self.assertNotEqual(list(summary["dataset_id"]), sorted(_DATASET_IDS_IN_MANIFEST_ORDER))


class ValidationAndLoadCountTests(BatchTestCase):
    def test_validation_runs_exactly_once_per_dataset(self) -> None:
        self._write_all_six_valid()
        with mock.patch(
            "adzuki_gwas_analysis.analysis.batch.ensure_validated_with_result",
            side_effect=batch.ensure_validated_with_result,
        ) as mocked:
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
            )
        self.assertEqual(mocked.call_count, 6)

    def test_dataframe_load_runs_exactly_once_per_dataset(self) -> None:
        self._write_all_six_valid()
        with mock.patch(
            "adzuki_gwas_analysis.analysis.batch.load_analysis_frame",
            side_effect=batch.load_analysis_frame,
        ) as mocked:
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
            )
        self.assertEqual(mocked.call_count, 6)

    def test_primary_pvalue_column_is_read_from_manifest(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        for result in outcome.datasets:
            self.assertEqual(result.pvalue_column, "pval")

    def test_n_tests_equals_validated_row_count_for_every_dataset(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        for result in outcome.datasets:
            self.assertEqual(result.n_tests, 3)


class MemoryDisciplineTests(BatchTestCase):
    def test_no_dataframe_from_a_previous_dataset_remains_referenced(self) -> None:
        self._write_all_six_valid()
        seen_refs: list[weakref.ReferenceType[pd.DataFrame]] = []
        original_loader = batch.load_analysis_frame

        def _tracking_loader(path: Path) -> pd.DataFrame:
            df = original_loader(path)
            seen_refs.append(weakref.ref(df))
            return df

        with mock.patch(
            "adzuki_gwas_analysis.analysis.batch.load_analysis_frame",
            side_effect=_tracking_loader,
        ):
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
            )

        self.assertEqual(len(seen_refs), 6)
        gc.collect()
        still_alive = [ref for ref in seen_refs if ref() is not None]
        self.assertEqual(
            still_alive,
            [],
            "a DataFrame from a previous dataset is still referenced after run_batch "
            "returned -- batch results must not keep any dataset's DataFrame alive",
        )

    def test_batch_dataset_result_holds_only_scalars_and_paths(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        for result in outcome.datasets:
            for field in dataclasses.fields(result):
                value = getattr(result, field.name)
                self.assertIsInstance(
                    value,
                    (str, int, float, type(None)),
                    f"BatchDatasetResult.{field.name} is a {type(value)!r}, not a scalar",
                )

    def test_batch_dataset_result_field_names_are_all_scalar_typed(self) -> None:
        # Structural cross-check against the module's own exported field-name tuple, so a
        # future field addition can't silently skip the scalar-only check above.
        self.assertEqual(
            batch.BATCH_DATASET_RESULT_FIELD_NAMES,
            tuple(f.name for f in dataclasses.fields(batch.BatchDatasetResult)),
        )


class PlotTitleTests(BatchTestCase):
    def test_manhattan_titles_for_all_six_datasets(self) -> None:
        self._write_all_six_valid()
        manifest = load_manifest(self.manifest_path)
        for entry in manifest.datasets:
            self.assertEqual(
                batch.format_manhattan_title(entry), _EXPECTED_MANHATTAN_TITLES[entry.dataset_id]
            )

    def test_qq_titles_for_all_six_datasets(self) -> None:
        self._write_all_six_valid()
        manifest = load_manifest(self.manifest_path)
        for entry in manifest.datasets:
            self.assertEqual(batch.format_qq_title(entry), _EXPECTED_QQ_TITLES[entry.dataset_id])

    def test_titles_are_passed_to_the_plotting_functions(self) -> None:
        # Assert on the title argument actually given to plot_manhattan/plot_qq (not by
        # OCR-reading a PNG, and not by re-deriving the title a second time), for one
        # representative dataset.
        self._write_all_six_valid()
        seen_manhattan_titles = []
        seen_qq_titles = []
        original_plot_manhattan = batch.plot_manhattan
        original_plot_qq = batch.plot_qq

        def _spy_manhattan(*args: object, **kwargs: object) -> object:
            seen_manhattan_titles.append(kwargs["title"])
            return original_plot_manhattan(*args, **kwargs)  # type: ignore[arg-type]

        def _spy_qq(*args: object, **kwargs: object) -> object:
            seen_qq_titles.append(kwargs["title"])
            return original_plot_qq(*args, **kwargs)  # type: ignore[arg-type]

        with (
            mock.patch(
                "adzuki_gwas_analysis.analysis.batch.plot_manhattan", side_effect=_spy_manhattan
            ),
            mock.patch("adzuki_gwas_analysis.analysis.batch.plot_qq", side_effect=_spy_qq),
        ):
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
            )
        self.assertEqual(
            seen_manhattan_titles,
            [_EXPECTED_MANHATTAN_TITLES[d] for d in _DATASET_IDS_IN_MANIFEST_ORDER],
        )
        self.assertEqual(
            seen_qq_titles, [_EXPECTED_QQ_TITLES[d] for d in _DATASET_IDS_IN_MANIFEST_ORDER]
        )

    def test_trait_label_mapping_is_exact(self) -> None:
        self.assertEqual(batch.format_trait_label("water_permeability"), "Water Permeability")
        self.assertEqual(batch.format_trait_label("red_seedcoat"), "Red Seed Coat Color")
        self.assertEqual(
            batch.format_trait_label("mottled_black_seedcoat"), "Mottled Black Seed Coat Color"
        )

    def test_unknown_trait_raises_rather_than_guessing(self) -> None:
        with self.assertRaises(ValueError):
            batch.format_trait_label("not_a_real_trait")


class OutputStructureTests(BatchTestCase):
    def test_exactly_twenty_five_files_are_produced(self) -> None:
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        all_files = [p for p in self.output_dir.rglob("*") if p.is_file()]
        self.assertEqual(len(all_files), 25)

    def test_each_dataset_directory_has_exactly_four_files(self) -> None:
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        for dataset_id in _DATASET_IDS_IN_MANIFEST_ORDER:
            dataset_dir = self.output_dir / dataset_id
            files = {p.name for p in dataset_dir.iterdir() if p.is_file()}
            self.assertEqual(
                files,
                {
                    f"{dataset_id}_manhattan.png",
                    f"{dataset_id}_qq.png",
                    "statistical_diagnostics.tsv",
                    "significant_variants.tsv",
                },
            )

    def test_root_batch_summary_exists(self) -> None:
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        self.assertTrue((self.output_dir / "batch_summary.tsv").is_file())

    def test_dataset_outputs_do_not_mix(self) -> None:
        # Each dataset's own statistical_diagnostics.tsv reports its own dataset_id, not a
        # neighboring dataset's.
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        for dataset_id in _DATASET_IDS_IN_MANIFEST_ORDER:
            table = pd.read_csv(
                self.output_dir / dataset_id / "statistical_diagnostics.tsv", sep="\t"
            )
            self.assertEqual(table["dataset_id"].iloc[0], dataset_id)

    def test_artifact_paths_in_summary_are_relative_posix_and_not_absolute(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        summary = pd.read_csv(outcome.summary_path, sep="\t")
        for column in (
            "manhattan_path",
            "qq_path",
            "diagnostics_path",
            "significant_variants_path",
        ):
            for value in summary[column]:
                self.assertFalse(value.startswith("/"))
                self.assertNotIn(str(self.output_dir), value)
                self.assertNotIn("\\", value)
                self.assertTrue((self.output_dir / value).is_file())

    def test_summary_column_order(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        summary = pd.read_csv(outcome.summary_path, sep="\t")
        self.assertEqual(
            list(summary.columns),
            [
                "schema_version",
                "dataset_id",
                "reference",
                "trait",
                "source_sha256",
                "pvalue_column",
                "pvalue_semantics",
                "family_scope",
                "n_tests",
                "visualization_threshold",
                "alpha",
                "bonferroni_threshold",
                "bonferroni_discoveries",
                "fdr_level",
                "bh_raw_p_cutoff",
                "bh_discoveries",
                "lambda_gc_df",
                "expected_chi2_median",
                "lambda_gc",
                "manhattan_path",
                "qq_path",
                "diagnostics_path",
                "significant_variants_path",
            ],
        )

    def test_summary_has_exactly_six_rows(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        summary = pd.read_csv(outcome.summary_path, sep="\t")
        self.assertEqual(len(summary), 6)

    def test_family_scope_is_independent_per_dataset_row(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        summary = pd.read_csv(outcome.summary_path, sep="\t")
        scopes = list(summary["family_scope"])
        self.assertEqual(len(scopes), len(set(scopes)), "family_scope must differ per dataset")
        for dataset_id, scope in zip(summary["dataset_id"], scopes, strict=True):
            self.assertIn(dataset_id, scope)

    def test_total_n_tests_is_sum_and_equals_18_for_this_fixture(self) -> None:
        self._write_all_six_valid()
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        self.assertEqual(outcome.total_n_tests, sum(r.n_tests for r in outcome.datasets))
        self.assertEqual(outcome.total_n_tests, 18)  # 6 datasets x 3 rows each, this fixture only

    def test_significant_variants_header_only_when_zero_discoveries(self) -> None:
        self._write_all_six_valid(tiny_p_for_first_dataset=False)
        batch.run_batch(
            manifest_path=self.manifest_path,
            data_dir=self.data_dir,
            output_dir=self.output_dir,
            alpha=1e-300,
            fdr_level=1e-300,
        )
        for dataset_id in _DATASET_IDS_IN_MANIFEST_ORDER:
            table = pd.read_csv(self.output_dir / dataset_id / "significant_variants.tsv", sep="\t")
            self.assertEqual(len(table), 0)
            self.assertEqual(
                list(table.columns),
                [
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
                ],
            )


class OutputDirSafetyTests(BatchTestCase):
    def test_rejects_existing_nonempty_output_dir(self) -> None:
        self._write_all_six_valid()
        self.output_dir.mkdir(parents=True)
        (self.output_dir / "stray.txt").write_text("pre-existing", encoding="utf-8")
        with self.assertRaises(Exception):  # noqa: B017 - BatchOutputDirectoryUnsafeError
            batch.run_batch(
                manifest_path=self.manifest_path,
                data_dir=self.data_dir,
                output_dir=self.output_dir,
            )
        # The pre-existing content must survive untouched.
        self.assertTrue((self.output_dir / "stray.txt").is_file())

    def test_accepts_existing_empty_output_dir(self) -> None:
        self._write_all_six_valid()
        self.output_dir.mkdir(parents=True)
        outcome = batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        self.assertEqual(len(outcome.datasets), 6)

    def test_rejects_symlink_output_dir(self) -> None:
        self._write_all_six_valid()
        real_dir = self.tmp_path / "real"
        real_dir.mkdir()
        symlink_dir = self.tmp_path / "link"
        symlink_dir.symlink_to(real_dir, target_is_directory=True)
        with self.assertRaises(Exception):  # noqa: B017 - BatchOutputDirectoryUnsafeError
            batch.run_batch(
                manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=symlink_dir
            )


class HeadlessRenderingTests(BatchTestCase):
    def test_succeeds_under_agg_backend(self) -> None:
        # This suite already runs with MPLBACKEND=Agg in CI; assert explicitly that batch
        # completes and both PNGs per dataset are non-empty files.
        self._write_all_six_valid()
        batch.run_batch(
            manifest_path=self.manifest_path, data_dir=self.data_dir, output_dir=self.output_dir
        )
        for dataset_id in _DATASET_IDS_IN_MANIFEST_ORDER:
            manhattan = self.output_dir / dataset_id / f"{dataset_id}_manhattan.png"
            qq = self.output_dir / dataset_id / f"{dataset_id}_qq.png"
            self.assertGreater(manhattan.stat().st_size, 0)
            self.assertGreater(qq.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
