import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.context import main, run_context
from adzuki_gwas_analysis.provenance import validate_provenance
from adzuki_gwas_analysis.tables import read_tsv, write_tsv
from tests.reference_support import write_bundle


class ContextTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle = write_bundle(self.root / "ref")
        self.candidates = self.root / "candidates.tsv"
        self.output = self.root / "out"
        self.rows = [
            {
                "candidate_id": "edge id",
                "dataset_id": "synthetic_trait",
                "reference": "Synthetic",
                "assembly_id": "synthetic-v1",
                "trait": "trait",
                "chr": "chrA",
                "pos": 1,
                "ref": "A",
                "alt": "C",
                "strand": "+",
            },
            {
                "candidate_id": "reverse",
                "dataset_id": "synthetic_trait",
                "reference": "Synthetic",
                "assembly_id": "synthetic-v1",
                "trait": "trait",
                "chr": "chrA",
                "pos": 9,
                "ref": "A",
                "alt": "G",
                "strand": "-",
            },
        ]

    def run_context(self, **options):
        write_tsv(
            self.candidates,
            self.rows[0] if self.rows else ("candidate_id", "chr", "pos", "ref", "alt", "trait"),
            self.rows,
        )
        run_context(
            candidate_table=self.candidates,
            bundle_path=self.bundle,
            dataset_id="synthetic_trait",
            output_dir=self.output,
            flank_bp=5,
            neighbor_window_bp=4,
            **options,
        )

    def test_edges_reverse_and_no_neighbors(self):
        self.run_context()
        rows = list(read_tsv(self.output / "candidate_context.tsv"))
        self.assertEqual(rows[0]["interval_start"], "1")
        self.assertEqual(rows[0]["flanks_complete"], "False")
        self.assertEqual(rows[1]["snp_offset_0based"], "5")
        self.assertIn(">edge%20id ", (self.output / "flanking_sequences.fasta").read_text())
        self.assertEqual(list(read_tsv(self.output / "neighboring_variants.tsv")), [])
        validate_provenance(self.output, required=True)

    def test_source_separation_and_indel_overlap(self):
        summary = self.root / "summary.tsv"
        summary.write_text("chr\tpos\tallele1\tallele0\nchrA\t3\tA\tG\n")
        vcf = self.root / "cohort.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n##reference=synthetic-v1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chrA\t4\t.\tTA\tT\t.\tPASS\t.\n"
        )
        self.run_context(summary_table=summary, cohort_vcf=vcf)
        rows = list(read_tsv(self.output / "neighboring_variants.tsv"))
        self.assertEqual({row["source"] for row in rows}, {"summary", "cohort_vcf"})
        self.assertTrue(all(len(row["source_sha256"]) == 64 for row in rows))
        self.assertEqual(sum(row["candidate_id"] == "reverse" for row in rows), 1)

    def test_empty_candidates(self):
        self.rows = []
        self.run_context()
        self.assertEqual(list(read_tsv(self.output / "candidate_context.tsv")), [])
        self.assertEqual((self.output / "flanking_sequences.fasta").read_text(), "")

    def test_mixed_reference_or_wrong_ref_fails_without_output(self):
        for key, value in (("reference", "Shumari"), ("assembly_id", "other"), ("ref", "T")):
            original = self.rows[0][key]
            self.rows[0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.run_context()
            self.assertFalse(self.output.exists())
            self.rows[0][key] = original

    def test_vcf_requires_matching_assembly(self):
        vcf = self.root / "bad.vcf"
        vcf.write_text("##reference=other\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        with self.assertRaises(ValueError):
            self.run_context(cohort_vcf=vcf)
        self.assertFalse(self.output.exists())

    def test_cli_and_legacy_alleles(self):
        self.candidates.write_text(
            "dataset_id\treference\ttrait\tchr\tpos\tallele1\tallele0\n"
            "synthetic_trait\tSynthetic\ttrait\tchrA\t20\tT\tC\n"
        )
        result = main(
            [
                "--candidate-table",
                str(self.candidates),
                "--bundle-path",
                str(self.bundle),
                "--dataset-id",
                "synthetic_trait",
                "--output-dir",
                str(self.output),
                "--flank-bp",
                "10",
                "--neighbor-window-bp",
                "10",
            ]
        )
        self.assertEqual(result, 0)
        row = next(read_tsv(self.output / "candidate_context.tsv"))
        self.assertEqual((row["ref"], row["alt"]), ("T", "C"))
        self.assertIn("effect_orientation_not_inferred", row["allele_source"])
