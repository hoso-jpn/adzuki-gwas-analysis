import json
import math
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.customer_summary import (
    adjust_and_cluster,
    load_metadata,
    neg_log10_pvalue,
    normalize_row,
    run_customer_summary,
)
from adzuki_gwas_analysis.provenance import validate_provenance
from adzuki_gwas_analysis.reference import load_reference_bundle
from adzuki_gwas_analysis.tables import read_tsv
from tests.reference_support import write_bundle
from tests.summary_support import SUMMARY_HEADER, write_summary_metadata


class SummaryFixture:
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle_path = write_bundle(self.root / "ref")
        self.bundle = load_reference_bundle(self.bundle_path)
        self.metadata = write_summary_metadata(self.root / "metadata.toml")
        self.meta = load_metadata(self.metadata, self.bundle)
        self.row = dict(
            CHROM="chrA",
            BP="1",
            REF="A",
            ALT="C",
            EA="A",
            OA="C",
            B="2",
            SE="0.5",
            P="1e-20",
            N="25",
            EAF="0.8",
        )

    def normalize(self):
        return normalize_row(self.row, self.meta, self.bundle, row_number=2, source_sha256="a" * 64)


class SummaryCoreTests(SummaryFixture, unittest.TestCase):
    def test_declared_allele_swap_and_frequency(self):
        result = self.normalize()
        self.assertEqual(result["beta"], -2)
        self.assertAlmostEqual(result["effect_allele_frequency"], 0.2)
        self.assertEqual(result["effect_allele"], "C")
        self.assertEqual(result["source_effect_allele"], "A")

    def test_unknown_strand_and_palindrome_remain_unresolved(self):
        self.meta["effect_strand"] = "unknown"
        self.meta["trait_coding"] = "unknown"
        self.row.update(ALT="T", EA="A", OA="T")
        result = self.normalize()
        self.assertEqual(result["beta"], 2)
        self.assertEqual(result["effect_orientation"], "unresolved_strand")
        self.assertEqual(result["favorable_allele"], "not_inferred")

    def test_invalid_values_are_not_coerced(self):
        for key, value in (
            ("SE", "0"),
            ("B", "NaN"),
            ("N", "2.5"),
            ("P", "0"),
            ("EAF", "1.1"),
            ("REF", "T"),
            ("OA", "G"),
        ):
            original = self.row[key]
            self.row[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.normalize()
            self.row[key] = original

    def test_extreme_and_zero_pvalue_contract(self):
        self.assertEqual(neg_log10_pvalue("1e-400")[0], 400)
        self.assertEqual(neg_log10_pvalue("0", "450")[0], 450)
        for p, log in (("0", None), ("-1", None), ("Inf", None), ("0.1", "2")):
            with self.subTest(p=p), self.assertRaises(ValueError):
                neg_log10_pvalue(p, log)

    def test_log_space_bh_and_family(self):
        rows = [
            {
                "candidate_id": str(i),
                "neg_log10_pvalue": -math.log10(p),
                "chr": "chrA",
                "pos": i + 1,
            }
            for i, p in enumerate((0.01, 0.04, 0.03, 0.2))
        ]
        selected = adjust_and_cluster(
            rows, alpha=0.05, fdr_level=0.05, clustering_distance=10, contig_order=["chrA"]
        )
        for row, expected in zip(
            rows, (0.04, 0.05333333333333334, 0.05333333333333334, 0.2), strict=True
        ):
            self.assertAlmostEqual(10 ** -row["neg_log10_p_bh"], expected)
        self.assertEqual(len(selected), 1)

    def test_unknown_test_cannot_declare_chi_square_df(self):
        write_summary_metadata(
            self.metadata, test="unknown", chi_square_df=1, chi_square_basis="unknown"
        )
        with self.assertRaises(ValueError):
            load_metadata(self.metadata, self.bundle)


class CustomerSummaryIntegrationTests(SummaryFixture, unittest.TestCase):
    def test_custom_columns_to_report_and_context_ready_candidates(self):
        summary = self.root / "summary.tsv"
        summary.write_text(
            SUMMARY_HEADER + "chrA\t1\tA\tC\tA\tC\t2\t0.5\t1e-400\t25\t0.8\n"
            "chrA\t2\tC\tG\tG\tC\t1\t0.2\t0.7\t25\t0.3\n"
        )
        output = self.root / "out"
        run_customer_summary(
            summary=summary,
            metadata_path=self.metadata,
            bundle_path=self.bundle_path,
            output_dir=output,
            clustering_distance=50,
        )
        rows = list(read_tsv(output / "normalized_summary.tsv"))
        self.assertEqual(len(rows), 2)
        self.assertEqual(float(rows[0]["neg_log10_pvalue"]), 400)
        diagnostics = json.loads((output / "statistical_diagnostics.json").read_text())
        self.assertEqual(diagnostics["n_tests"], 2)
        self.assertEqual(diagnostics["lambda_gc"]["status"], "not_assessed")
        validate_provenance(output, required=True)
        self.assertTrue((output / "manhattan.png").is_file())
        self.assertTrue((output / "analysis_report.md").is_file())

    def test_exclusion_accounting_and_duplicates(self):
        write_summary_metadata(self.metadata, invalid_row_policy="exclude")
        summary = self.root / "summary.tsv"
        good = "chrA\t1\tA\tC\tA\tC\t2\t0.5\t0.001\t25\t0.8\n"
        summary.write_text(SUMMARY_HEADER + good + "chrA\t2\tC\tG\tG\tC\t1\t0\t0.7\t25\t0.3\n")
        output = self.root / "out"
        run_customer_summary(
            summary=summary,
            metadata_path=self.metadata,
            bundle_path=self.bundle_path,
            output_dir=output,
            clustering_distance=50,
        )
        record = validate_provenance(output, required=True)
        self.assertEqual(
            (record["parameters"]["n_input"], record["parameters"]["n_excluded"]), (2, 1)
        )
        summary.write_text(SUMMARY_HEADER + good + good)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            run_customer_summary(
                summary=summary,
                metadata_path=self.metadata,
                bundle_path=self.bundle_path,
                output_dir=self.root / "duplicate",
                clustering_distance=50,
            )
