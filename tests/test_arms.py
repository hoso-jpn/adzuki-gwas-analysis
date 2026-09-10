import random
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.arms import PrimerSettings, design_candidate, run_arms
from adzuki_gwas_analysis.context import run_context
from adzuki_gwas_analysis.primer_metrics import complementary_run, tm_nn
from adzuki_gwas_analysis.provenance import validate_provenance
from adzuki_gwas_analysis.reference import reverse_complement
from adzuki_gwas_analysis.tables import read_tsv
from tests.reference_support import write_bundle


class PrimerMetricsTests(unittest.TestCase):
    def test_published_biopython_nearest_neighbor_example(self):
        # Official Bio.SeqUtils.MeltingTemp documentation: DNA_NN3, Na=50, dnac1=dnac2=25.
        self.assertAlmostEqual(tm_nn("CGTTCCAAAGATGTGGGCATGAGCTTAC"), 60.32, places=2)
        self.assertGreater(tm_nn("CGTTCCAAAGATGTGGGCATGAGCTTAC", sodium_mM=100), 60.32)

    def test_complementarity_and_invalid_conditions(self):
        self.assertEqual(complementary_run("AAAA", "TTTT"), 4)
        self.assertEqual(complementary_run("AAAA", "AAAA"), 0)
        with self.assertRaises(ValueError):
            tm_nn("ACNT")
        with self.assertRaises(ValueError):
            tm_nn("ACGT", sodium_mM=0)
        with self.assertRaises(ValueError):
            PrimerSettings(min_length=30, max_length=20).validate()


class ArmsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        generator = random.Random(1)
        self.sequence = "".join(generator.choice("ACGT") for _ in range(500))
        self.bundle = write_bundle(self.root / "ref", sequence=self.sequence)
        self.ref = self.sequence[200]
        self.alt = next(base for base in "ACGT" if base != self.ref)
        self.candidates = self.root / "candidates.tsv"
        self.candidates.write_text(
            "candidate_id\tchr\tpos\tref\talt\ttrait\n"
            f"c1\tchrA\t201\t{self.ref}\t{self.alt}\ttrait\n"
        )
        self.context = self.root / "context"
        run_context(
            candidate_table=self.candidates,
            bundle_path=self.bundle,
            dataset_id="synthetic_trait",
            output_dir=self.context,
            flank_bp=180,
            neighbor_window_bp=180,
        )
        self.config = self.root / "primer.toml"
        self.config.write_text(
            "schema_version = 1\n[primer]\nmin_length = 18\nmax_length = 20\n"
            "min_product_bp = 75\nmax_product_bp = 90\nmin_tm_C = 35\n"
            "max_tm_C = 80\nmax_tm_difference_C = 12\n"
        )

    def test_design_e2e_multiple_alternatives_and_three_prime_alleles(self):
        output = self.root / "design"
        run_arms(context_dir=self.context, config=self.config, output_dir=output)
        rows = list(read_tsv(output / "primer_candidates.tsv"))
        self.assertEqual(len(rows), 3)
        self.assertEqual([row["rank"] for row in rows], ["1", "2", "3"])
        for row in rows:
            self.assertEqual(row["specific_ref_5to3"][-1], self.ref)
            self.assertEqual(row["specific_alt_5to3"][-1], self.alt)
            self.assertEqual(row["specificity"], "unverified")
            self.assertEqual(row["cohort_binding_check"], "unavailable")
            self.assertEqual(
                int(row["common_end"]) - int(row["specific_start"]) + 1, int(row["product_bp"])
            )
        validate_provenance(output, required=True)

    def test_context_tamper_is_rejected(self):
        (self.context / "flanking_sequences.fasta").write_text(">wrong\nAAAA\n")
        with self.assertRaisesRegex(ValueError, "checksum"):
            run_arms(context_dir=self.context, config=self.config, output_dir=self.root / "bad")

    def test_reverse_strand_coordinates_and_partial_cohort_coverage(self):
        row = next(read_tsv(self.context / "candidate_context.tsv"))
        row["strand"] = "-"
        settings = PrimerSettings(
            min_length=18,
            max_length=20,
            min_product_bp=75,
            max_product_bp=90,
            min_tm_C=35,
            max_tm_C=80,
            max_tm_difference_C=12,
        )
        designs, _ = design_candidate(
            row,
            reverse_complement(self.sequence[20:381]),
            [],
            settings,
            cohort_available=True,
            neighbor_window=1,
        )
        self.assertTrue(designs)
        for design in designs:
            self.assertEqual(design["specific_ref_5to3"][-1], reverse_complement(self.ref))
            self.assertEqual(
                design["specific_end"] - design["common_start"] + 1, design["product_bp"]
            )
            self.assertEqual(design["cohort_binding_check"], "partial_window_unverified")

    def test_ineligible_and_short_flanks_have_reasons(self):
        row = next(read_tsv(self.context / "candidate_context.tsv"))
        sequence = self.sequence[20:381]
        row["alt"] = "CG"
        designs, reasons = design_candidate(
            row, sequence, [], PrimerSettings(), cohort_available=False, neighbor_window=0
        )
        self.assertFalse(designs)
        self.assertIn("not_biallelic_snp", reasons)
        row["alt"] = self.alt
        row["snp_offset_0based"] = "0"
        designs, reasons = design_candidate(
            row, self.ref + "ACGT", [], PrimerSettings(), cohort_available=False, neighbor_window=0
        )
        self.assertIn("insufficient_flanks", reasons)

    def test_cohort_overlap_excludes_designs(self):
        row = next(read_tsv(self.context / "candidate_context.tsv"))
        settings = PrimerSettings(
            min_length=18,
            max_length=20,
            min_product_bp=75,
            max_product_bp=90,
            min_tm_C=35,
            max_tm_C=80,
            max_tm_difference_C=12,
        )
        neighbor = {"source": "cohort_vcf", "pos": "200", "ref": self.sequence[199]}
        designs, reasons = design_candidate(
            row,
            self.sequence[20:381],
            [neighbor],
            settings,
            cohort_available=True,
            neighbor_window=180,
        )
        self.assertFalse(designs)
        self.assertIn("cohort_variant_in_binding_site", reasons)
