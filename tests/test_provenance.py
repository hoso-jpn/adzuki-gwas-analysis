import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from adzuki_gwas_analysis.errors import BatchOutputDirectoryUnsafeError
from adzuki_gwas_analysis.provenance import (
    CONTRACT_FILE,
    PROVENANCE_FILE,
    finish_provenance,
    generation_environment,
    input_checksums,
    output_transaction,
    run_audited_analysis,
    validate_provenance,
)


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "input.tsv"
        self.source.write_text("synthetic\n")
        self.inputs = {"summary": self.source}
        self.initial = input_checksums(self.inputs)

    def emit(self, destination):
        with output_transaction(destination) as stage:
            (stage / "result.tsv").write_text("value\n1\n")
            finish_provenance(
                stage,
                operation="synthetic",
                inputs=self.inputs,
                initial=self.initial,
                parameters={"alpha": 0.05},
                environment=generation_environment(),
                families=[],
            )

    def test_identity_is_deterministic_and_paths_are_not_recorded(self):
        first, second = self.root / "one", self.root / "two"
        self.emit(first)
        self.emit(second)
        self.assertEqual(
            (first / PROVENANCE_FILE).read_bytes(), (second / PROVENANCE_FILE).read_bytes()
        )
        self.assertNotIn(str(self.root), (first / PROVENANCE_FILE).read_text())
        self.assertIsNotNone(validate_provenance(first, required=True))

    def test_tampering_and_missing_metadata(self):
        target = self.root / "out"
        self.emit(target)
        (target / "result.tsv").write_text("value\n2\n")
        with self.assertRaisesRegex(ValueError, "checksum"):
            validate_provenance(target)
        (target / PROVENANCE_FILE).unlink()
        with self.assertRaisesRegex(ValueError, "incomplete"):
            validate_provenance(target)
        (target / CONTRACT_FILE).unlink()
        self.assertIsNone(validate_provenance(target))
        with self.assertRaises(ValueError):
            validate_provenance(target, required=True)

    def test_changed_input_aborts_and_cleans_stage(self):
        self.source.write_text("changed\n")
        with self.assertRaisesRegex(ValueError, "input changed"):
            self.emit(self.root / "out")
        self.assertFalse((self.root / "out").exists())
        self.assertFalse(list(self.root.glob(".analysis-stage-*")))

    def test_failure_does_not_replace_existing_results(self):
        target = self.root / "out"
        self.emit(target)
        original = (target / PROVENANCE_FILE).read_bytes()
        with self.assertRaises(BatchOutputDirectoryUnsafeError):
            self.emit(target)
        self.assertEqual((target / PROVENANCE_FILE).read_bytes(), original)
        with (
            self.assertRaisesRegex(RuntimeError, "injected"),
            output_transaction(self.root / "failed") as stage,
        ):
            (stage / "partial").write_text("partial")
            raise RuntimeError("injected")
        self.assertFalse((self.root / "failed").exists())


class AuditedAnalysisTests(unittest.TestCase):
    def test_batch_to_report_preserves_original_environment(self):
        from adzuki_gwas_analysis.analysis.report import run_report
        from tests.analysis_support import build_candidate_enabled_batch_fixture

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            legacy = build_candidate_enabled_batch_fixture(root)
            audited, delivery = root / "audited", root / "delivery"
            run_audited_analysis(
                command="batch",
                manifest_path=root / "fixture_manifest.toml",
                data_dir=root / "fixture_data",
                output_dir=audited,
                dataset_id=None,
                alpha=0.05,
                fdr_level=0.05,
                clustering_distance=1000,
            )
            for path in legacy.rglob("*"):
                if path.is_file():
                    self.assertEqual(
                        path.read_bytes(), (audited / path.relative_to(legacy)).read_bytes()
                    )
            generation = validate_provenance(audited, required=True)
            self.assertEqual(len(generation["families"]), 6)
            with patch(
                "adzuki_gwas_analysis.analysis.report._report_generation_environment",
                return_value={"python_version": "different-report-interpreter"},
            ):
                run_report(analysis_dir=audited, output_dir=delivery, require_provenance=True)
            versions = json.loads((delivery / "reproducibility/software_versions.json").read_text())
            self.assertEqual(versions["analysis_generation_environment"], generation["environment"])
            self.assertEqual(
                versions["report_generation_environment"]["python_version"],
                "different-report-interpreter",
            )
            validate_provenance(delivery / "artifacts", required=True)
            with self.assertRaises(ValueError):
                run_report(
                    analysis_dir=legacy, output_dir=root / "rejected", require_provenance=True
                )

    def test_candidate_cli_records_provenance(self):
        from adzuki_gwas_analysis.analysis.cli import main
        from tests.analysis_support import write_manifest

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # The helper binds the canonical member name; materialize this synthetic file only.
            import shutil

            from tests.analysis_support import MIYAGI_WATER_PERMEABILITY_FILENAME

            shutil.copyfile(
                "tests/fixtures/valid.assoc.txt", root / MIYAGI_WATER_PERMEABILITY_FILENAME
            )
            manifest = write_manifest(
                root / "manifest.toml", root / MIYAGI_WATER_PERMEABILITY_FILENAME
            )
            result = main(
                [
                    "candidates",
                    "--manifest",
                    str(manifest),
                    "--data-dir",
                    str(root),
                    "--output-dir",
                    str(root / "out"),
                    "--clustering-distance",
                    "1000",
                    "--record-provenance",
                ]
            )
            self.assertEqual(result, 0)
            record = validate_provenance(root / "out", required=True)
            self.assertEqual(len(record["families"]), 1)
