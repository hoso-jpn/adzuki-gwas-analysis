"""Unit tests for adzuki_gwas_analysis.cli, exercised in-process (no subprocess)."""

from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.cli import main
from adzuki_gwas_analysis.loader import compute_sha256, count_data_rows

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def _manifest_toml_for(filename: str) -> str:
    path = FIXTURES_DIR / filename
    sha256 = compute_sha256(path)
    row_count = count_data_rows(path)
    return f"""
schema_version = 1

[dryad]
doi = "10.5061/dryad.8w9ghx3xv"
dataset_id = 149675
version_id = 356599
version_number = 6
publication_doi = "10.1126/science.ads2871"

[archive]
filename = "adzuki_GWAS_data.zip"
size_bytes = 257167472
sha256 = "92f0bc6a235c3ee863e96736a831f8652569269a8f61c8699cc9f562743b0ed5"

[pvalue_columns]
p_wald = "Wald test p-value"
pval = "Likelihood ratio test (LRT) p-value"
p_score = "Score test p-value"
primary = "pval"

[[datasets]]
dataset_id = "cli_test_dataset"
reference = "Miyagi"
trait = "water_permeability"
member_filename = "{filename}"
member_sha256 = "{sha256}"
row_count = {row_count}
"""


class CliMainTests(unittest.TestCase):
    def test_success_exits_zero_and_reports_all_passed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            manifest_path = Path(tmp_dir) / "manifest.toml"
            manifest_path.write_text(_manifest_toml_for("valid.assoc.txt"), encoding="utf-8")

            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "--manifest",
                        str(manifest_path),
                        "--data-dir",
                        str(FIXTURES_DIR),
                    ]
                )
        self.assertEqual(exit_code, 0)
        self.assertIn("All 1 dataset(s) passed validation.", stdout.getvalue())
        self.assertIn("[OK]", stdout.getvalue())

    def test_failure_exits_nonzero(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            manifest_path = Path(tmp_dir) / "manifest.toml"
            manifest_path.write_text(
                _manifest_toml_for("bad_numeric_pos.assoc.txt"), encoding="utf-8"
            )

            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                exit_code = main(
                    [
                        "--manifest",
                        str(manifest_path),
                        "--data-dir",
                        str(FIXTURES_DIR),
                    ]
                )
        self.assertEqual(exit_code, 1)
        self.assertIn("[FAIL]", stdout.getvalue())
        self.assertIn("failed validation", stderr.getvalue())

    def test_missing_manifest_exits_nonzero_without_a_false_success_line(self) -> None:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            exit_code = main(
                [
                    "--manifest",
                    "/nonexistent/manifest.toml",
                    "--data-dir",
                    str(FIXTURES_DIR),
                ]
            )
        self.assertEqual(exit_code, 1)
        self.assertNotIn("passed validation", stdout.getvalue())
        self.assertIn("ERROR", stderr.getvalue())

    def test_output_flag_writes_machine_readable_json_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            manifest_path = Path(tmp_dir) / "manifest.toml"
            manifest_path.write_text(_manifest_toml_for("valid.assoc.txt"), encoding="utf-8")
            output_path = Path(tmp_dir) / "result.json"

            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "--manifest",
                        str(manifest_path),
                        "--data-dir",
                        str(FIXTURES_DIR),
                        "--output",
                        str(output_path),
                    ]
                )

            self.assertEqual(exit_code, 0)
            payload = json.loads(output_path.read_text(encoding="utf-8"))
        self.assertEqual(len(payload), 1)
        self.assertEqual(payload[0]["dataset_id"], "cli_test_dataset")
        self.assertTrue(payload[0]["success"])
        self.assertEqual(payload[0]["row_count"], 3)


if __name__ == "__main__":
    unittest.main()
