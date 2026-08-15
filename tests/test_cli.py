"""Unit tests for adzuki_gwas_analysis.cli, exercised in-process (no subprocess).

Schema v1 requires a manifest to declare exactly the 6 canonical
(reference, trait) datasets (PR #2 review, P1-1), so every CLI test that
expects a *successful* run must supply all 6 -- a single-dataset manifest is
no longer a valid fixture for that case. ``_prepare_data_dir`` copies the
small ``valid.assoc.txt`` fixture's content under each of the 6 canonical
member filenames (optionally substituting a different fixture for one
dataset, to test a single failing dataset among 5 passing ones), and
``_build_manifest_toml`` builds a matching manifest from the real
checksums/row counts of what was just written -- nothing here is a
hand-copied value that could drift from the fixture files.
"""

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

# (dataset_id, reference, trait, member_filename) for all 6 schema v1 datasets.
CANONICAL_ENTRIES: tuple[tuple[str, str, str, str], ...] = (
    (
        "miyagi_water_permeability",
        "Miyagi",
        "water_permeability",
        "mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt",
    ),
    (
        "miyagi_red_seedcoat",
        "Miyagi",
        "red_seedcoat",
        "mapped_to_Miyagi_red_seedcoat.maf_0.05.assoc.txt",
    ),
    (
        "miyagi_mottled_black_seedcoat",
        "Miyagi",
        "mottled_black_seedcoat",
        "mapped_to_Miyagi_mottled_black_seedcoat.maf_0.05.assoc.txt",
    ),
    (
        "shumari_water_permeability",
        "Shumari",
        "water_permeability",
        "mapped_to_Shumari_water_permeability.maf_0.05.assoc.txt",
    ),
    (
        "shumari_red_seedcoat",
        "Shumari",
        "red_seedcoat",
        "mapped_to_Shumari_red_seedcoat.maf_0.05.assoc.txt",
    ),
    (
        "shumari_mottled_black_seedcoat",
        "Shumari",
        "mottled_black_seedcoat",
        "mapped_to_Shumari_mottled_black_seedcoat.maf_0.05.assoc.txt",
    ),
)


def _prepare_data_dir(
    data_dir: Path,
    entries: tuple[tuple[str, str, str, str], ...],
    *,
    fixture_overrides: dict[str, str] | None = None,
) -> dict[str, tuple[str, int]]:
    """Copy fixture content into ``data_dir`` under each canonical member filename.

    Returns ``{dataset_id: (sha256, row_count)}`` computed from what was
    actually written, for use building a matching manifest.
    """
    fixture_overrides = fixture_overrides or {}
    info: dict[str, tuple[str, int]] = {}
    for dataset_id, _reference, _trait, member_filename in entries:
        source = FIXTURES_DIR / fixture_overrides.get(dataset_id, "valid.assoc.txt")
        dest = data_dir / member_filename
        dest.write_bytes(source.read_bytes())
        info[dataset_id] = (compute_sha256(dest), count_data_rows(dest))
    return info


def _build_manifest_toml(
    entries: tuple[tuple[str, str, str, str], ...], info: dict[str, tuple[str, int]]
) -> str:
    blocks = []
    for dataset_id, reference, trait, member_filename in entries:
        sha256, row_count = info[dataset_id]
        blocks.append(
            f"""
[[datasets]]
dataset_id = "{dataset_id}"
reference = "{reference}"
trait = "{trait}"
member_filename = "{member_filename}"
member_sha256 = "{sha256}"
row_count = {row_count}
"""
        )
    return """
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
""" + "".join(blocks)


class CliMainTests(unittest.TestCase):
    def test_success_exits_zero_and_reports_all_six_passed(self) -> None:
        with tempfile.TemporaryDirectory() as data_dir_str:
            data_dir = Path(data_dir_str)
            info = _prepare_data_dir(data_dir, CANONICAL_ENTRIES)
            manifest_path = data_dir / "manifest.toml"
            manifest_path.write_text(
                _build_manifest_toml(CANONICAL_ENTRIES, info), encoding="utf-8"
            )

            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                exit_code = main(["--manifest", str(manifest_path), "--data-dir", str(data_dir)])
        self.assertEqual(exit_code, 0)
        self.assertIn("All 6 dataset(s) passed validation.", stdout.getvalue())
        self.assertEqual(stdout.getvalue().count("[OK]"), 6)

    def test_one_failing_dataset_among_six_exits_nonzero(self) -> None:
        with tempfile.TemporaryDirectory() as data_dir_str:
            data_dir = Path(data_dir_str)
            info = _prepare_data_dir(
                data_dir,
                CANONICAL_ENTRIES,
                fixture_overrides={"miyagi_water_permeability": "bad_numeric_pos.assoc.txt"},
            )
            manifest_path = data_dir / "manifest.toml"
            manifest_path.write_text(
                _build_manifest_toml(CANONICAL_ENTRIES, info), encoding="utf-8"
            )

            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                exit_code = main(["--manifest", str(manifest_path), "--data-dir", str(data_dir)])
        self.assertEqual(exit_code, 1)
        self.assertIn("[FAIL]", stdout.getvalue())
        self.assertEqual(stdout.getvalue().count("[OK]"), 5)
        self.assertIn("failed validation", stderr.getvalue())

    def test_missing_manifest_exits_nonzero_without_a_false_success_line(self) -> None:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            exit_code = main(
                ["--manifest", "/nonexistent/manifest.toml", "--data-dir", str(FIXTURES_DIR)]
            )
        self.assertEqual(exit_code, 1)
        self.assertNotIn("passed validation", stdout.getvalue())
        self.assertIn("ERROR", stderr.getvalue())

    def test_invalid_toml_manifest_exits_nonzero_with_diagnosable_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            manifest_path = Path(tmp_dir) / "manifest.toml"
            manifest_path.write_text("this is not [[[ valid toml", encoding="utf-8")

            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                exit_code = main(
                    ["--manifest", str(manifest_path), "--data-dir", str(FIXTURES_DIR)]
                )
        self.assertEqual(exit_code, 1)
        self.assertNotIn("passed validation", stdout.getvalue())
        self.assertIn("ERROR", stderr.getvalue())
        self.assertNotIn("Traceback", stderr.getvalue())

    def test_incomplete_manifest_exits_nonzero_and_never_reports_success(self) -> None:
        # Schema v1 requires all 6 datasets (P1-1); a 5-dataset manifest
        # must be rejected at manifest-load time, before any file is
        # touched, so no [OK]/[FAIL] lines or "passed validation" text
        # should appear at all.
        with tempfile.TemporaryDirectory() as data_dir_str:
            data_dir = Path(data_dir_str)
            incomplete_entries = CANONICAL_ENTRIES[:5]
            info = _prepare_data_dir(data_dir, incomplete_entries)
            manifest_path = data_dir / "manifest.toml"
            manifest_path.write_text(
                _build_manifest_toml(incomplete_entries, info), encoding="utf-8"
            )

            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                exit_code = main(["--manifest", str(manifest_path), "--data-dir", str(data_dir)])
        self.assertEqual(exit_code, 1)
        self.assertNotIn("passed validation", stdout.getvalue())
        self.assertNotIn("[OK]", stdout.getvalue())
        self.assertIn("ERROR", stderr.getvalue())

    def test_output_flag_writes_machine_readable_json_summary(self) -> None:
        with tempfile.TemporaryDirectory() as data_dir_str:
            data_dir = Path(data_dir_str)
            info = _prepare_data_dir(data_dir, CANONICAL_ENTRIES)
            manifest_path = data_dir / "manifest.toml"
            manifest_path.write_text(
                _build_manifest_toml(CANONICAL_ENTRIES, info), encoding="utf-8"
            )
            output_path = data_dir / "result.json"

            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "--manifest",
                        str(manifest_path),
                        "--data-dir",
                        str(data_dir),
                        "--output",
                        str(output_path),
                    ]
                )

            self.assertEqual(exit_code, 0)
            payload = json.loads(output_path.read_text(encoding="utf-8"))
        self.assertEqual(len(payload), 6)
        dataset_ids = {row["dataset_id"] for row in payload}
        self.assertEqual(dataset_ids, {entry[0] for entry in CANONICAL_ENTRIES})
        self.assertTrue(all(row["success"] for row in payload))


if __name__ == "__main__":
    unittest.main()
