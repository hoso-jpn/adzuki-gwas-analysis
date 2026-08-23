"""Unit tests for adzuki_gwas_analysis.analysis.report_content.

Pure content-builder tests -- no file I/O beyond the read-only fixture setup shared with
tests/test_analysis_report_validation.py. Uses a real, small candidate-enabled ``batch``
output (built via :func:`tests.analysis_support.build_candidate_enabled_batch_fixture`) run
through :func:`~adzuki_gwas_analysis.analysis.report_validation.validate_analysis_dir` so
every check exercises the real, already-validated data shape these builders actually
receive in production.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.analysis.report_content import (
    UNAVAILABLE_FROM_SOURCE_ARTIFACTS,
    build_analysis_report_markdown,
    build_executive_summary_markdown,
    build_run_manifest,
    build_software_versions,
)
from adzuki_gwas_analysis.analysis.report_validation import validate_analysis_dir
from tests.analysis_support import build_candidate_enabled_batch_fixture

#: Phrases that must never appear as an unqualified, affirmative claim -- every line
#: containing one of these must also contain a negation marker on the same line.
_RISK_PHRASES = (
    "LD block",
    "independent QTL",
    "causal variant",
    "causal SNP",
    "validated breeding marker",
    "breeding marker",
    "validated marker",
)
_NEGATION_MARKERS = ("not ", "Not ", "never ", "Never ", "no ", "No ", "**not**")


def _assert_no_unqualified_risk_claims(testcase: unittest.TestCase, text: str) -> None:
    for line in text.splitlines():
        for phrase in _RISK_PHRASES:
            if phrase in line:
                testcase.assertTrue(
                    any(marker in line for marker in _NEGATION_MARKERS),
                    f"line contains {phrase!r} with no negation marker: {line!r}",
                )


class ReportContentTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)
        analysis_dir = build_candidate_enabled_batch_fixture(self.tmp_path)
        self.validated = validate_analysis_dir(analysis_dir)


class ExecutiveSummaryTests(ReportContentTestCase):
    def test_no_unqualified_risk_claims(self) -> None:
        text = build_executive_summary_markdown(self.validated)
        _assert_no_unqualified_risk_claims(self, text)

    def test_mentions_all_six_dataset_ids(self) -> None:
        text = build_executive_summary_markdown(self.validated)
        for dataset in self.validated.datasets:
            self.assertIn(dataset.dataset_id, text)

    def test_zero_candidate_dataset_says_no_candidates(self) -> None:
        text = build_executive_summary_markdown(self.validated)
        self.assertIn("No candidates passed Bonferroni or Benjamini-Hochberg significance", text)

    def test_totals_are_sums_not_hardcoded(self) -> None:
        text = build_executive_summary_markdown(self.validated)
        total_signals = sum(d.n_signals for d in self.validated.datasets)
        total_candidates = sum(d.n_candidates for d in self.validated.datasets)
        self.assertIn(f"{total_signals:,} signal(s)", text)
        self.assertIn(f"{total_candidates:,} candidate(s)", text)

    def test_totals_are_explicitly_not_a_shared_family(self) -> None:
        text = build_executive_summary_markdown(self.validated)
        self.assertIn("not** a shared multiple-testing family", text)

    def test_preview_truncation_pointer_appears_when_a_dataset_exceeds_the_limit(self) -> None:
        text = build_executive_summary_markdown(self.validated)
        if any(d.n_candidates > 3 for d in self.validated.datasets):
            self.assertIn("more candidate(s)", text)
            self.assertIn("candidate_ranking.tsv", text)

    def test_no_absolute_paths(self) -> None:
        text = build_executive_summary_markdown(self.validated)
        self.assertNotIn(str(self.tmp_path), text)
        self.assertNotIn("/tmp", text)

    def test_deterministic_across_calls(self) -> None:
        first = build_executive_summary_markdown(self.validated)
        second = build_executive_summary_markdown(self.validated)
        self.assertEqual(first, second)


class AnalysisReportTests(ReportContentTestCase):
    def test_no_unqualified_risk_claims(self) -> None:
        text = build_analysis_report_markdown(self.validated)
        _assert_no_unqualified_risk_claims(self, text)

    def test_bonferroni_and_bh_purposes_are_distinguished(self) -> None:
        text = build_analysis_report_markdown(self.validated)
        self.assertIn("family-wise error rate", text.lower())
        self.assertIn("false discovery rate", text.lower())

    def test_lambda_gc_is_called_a_diagnostic_only(self) -> None:
        text = build_analysis_report_markdown(self.validated)
        self.assertIn("diagnostic statistic", text)

    def test_zero_candidate_dataset_section_present(self) -> None:
        text = build_analysis_report_markdown(self.validated)
        self.assertIn("### `miyagi_mottled_black_seedcoat`", text)
        self.assertIn("No candidates passed Bonferroni or Benjamini-Hochberg significance", text)

    def test_no_absolute_paths(self) -> None:
        text = build_analysis_report_markdown(self.validated)
        self.assertNotIn(str(self.tmp_path), text)

    def test_deterministic_across_calls(self) -> None:
        first = build_analysis_report_markdown(self.validated)
        second = build_analysis_report_markdown(self.validated)
        self.assertEqual(first, second)


class RunManifestTests(ReportContentTestCase):
    def test_schema_and_parameters_reflect_validated_input(self) -> None:
        manifest = build_run_manifest(self.validated, provenance={}, artifacts=[])
        self.assertEqual(manifest["source_batch_summary_schema_version"], 2)
        self.assertEqual(
            manifest["parameters"],
            {
                "alpha": self.validated.alpha,
                "fdr_level": self.validated.fdr_level,
                "visualization_threshold": self.validated.visualization_threshold,
                "clustering_distance": self.validated.clustering_distance,
            },
        )

    def test_datasets_array_has_one_entry_per_dataset(self) -> None:
        manifest = build_run_manifest(self.validated, provenance={}, artifacts=[])
        self.assertEqual(len(manifest["datasets"]), len(self.validated.datasets))

    def test_scientific_scope_never_asserts_ld_qtl_causal_or_marker_claims(self) -> None:
        manifest = build_run_manifest(self.validated, provenance={}, artifacts=[])
        scope = manifest["scientific_scope"]
        self.assertTrue(scope["post_hoc_analysis_of_published_summary_statistics"])
        self.assertFalse(scope["gwas_rerun"])
        self.assertFalse(scope["association_signal_is_ld_block"])
        self.assertFalse(scope["association_signal_is_independent_qtl"])
        self.assertFalse(scope["lead_variant_is_causal_variant_claim"])
        self.assertFalse(scope["candidate_is_causal_variant_claim"])
        self.assertFalse(scope["candidate_is_validated_breeding_marker_claim"])
        self.assertFalse(scope["priority_tier_is_biological_importance_ranking"])
        self.assertFalse(scope["cross_reference_coordinate_integration"])
        self.assertFalse(scope["cross_dataset_shared_candidate_ranking"])

    def test_artifacts_passed_through_verbatim(self) -> None:
        artifacts = [{"relative_path": "x", "sha256": "y", "role": "z"}]
        manifest = build_run_manifest(self.validated, provenance={}, artifacts=artifacts)
        self.assertEqual(manifest["artifacts"], artifacts)


class SoftwareVersionsTests(unittest.TestCase):
    def test_analysis_generation_environment_is_the_unavailable_sentinel(self) -> None:
        doc = build_software_versions(report_generation_environment={"python_version": "3.11"})
        self.assertEqual(doc["analysis_generation_environment"], UNAVAILABLE_FROM_SOURCE_ARTIFACTS)

    def test_report_generation_environment_passed_through_verbatim(self) -> None:
        env = {"python_version": "3.11.15", "platform": "Linux", "packages": {}}
        doc = build_software_versions(report_generation_environment=env)
        self.assertEqual(doc["report_generation_environment"], env)

    def test_never_claims_analysis_and_report_environments_are_the_same(self) -> None:
        doc = build_software_versions(report_generation_environment={"python_version": "3.11"})
        note = str(doc["note"])
        self.assertIn("not necessarily generated in this same process", note)


if __name__ == "__main__":
    unittest.main()
