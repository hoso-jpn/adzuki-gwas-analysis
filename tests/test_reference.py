import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.loader import compute_sha256
from adzuki_gwas_analysis.reference import ReferenceError, load_reference_bundle
from tests.reference_support import write_bundle


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = write_bundle(self.root)

    def test_indexed_sequence_edges_and_strands(self):
        bundle = load_reference_bundle(self.path)
        bundle.check_dataset("synthetic_trait", "Synthetic", "synthetic-v1")
        self.assertEqual(bundle.sequence("chrA", 39, 44), "GTACGT")
        self.assertEqual(bundle.sequence("chrA", 1, 3, "-"), "CGT")
        self.assertEqual(bundle.sequence("chrA", 320, 320), "T")
        bundle.check_snp("chrA", 1, "A", "C")

    def test_reject_identity_and_coordinates(self):
        bundle = load_reference_bundle(self.path)
        for dataset, reference, assembly in (
            ("other", "Synthetic", "synthetic-v1"),
            ("synthetic_trait", "Miyagi", "synthetic-v1"),
            ("synthetic_trait", "Synthetic", "other"),
        ):
            with (
                self.subTest(dataset=dataset, reference=reference, assembly=assembly),
                self.assertRaises(ReferenceError),
            ):
                bundle.check_dataset(dataset, reference, assembly)
        for chrom, pos, ref, alt in (
            ("chrZ", 1, "A", "C"),
            ("chrA", 0, "A", "C"),
            ("chrA", 321, "A", "C"),
            ("chrA", 1, "T", "C"),
            ("chrA", 1, "A", "A"),
            ("chrA", 1, "A", "C,G"),
        ):
            with (
                self.subTest(chrom=chrom, pos=pos, ref=ref, alt=alt),
                self.assertRaises(ReferenceError),
            ):
                bundle.check_snp(chrom, pos, ref, alt)

    def test_missing_and_changed_assets(self):
        fasta = self.root / "reference.fasta"
        fasta.write_text(">chrA\nAAAA\n")
        with self.assertRaisesRegex(ReferenceError, "checksum"):
            load_reference_bundle(self.path)
        fasta.unlink()
        with self.assertRaises(ReferenceError):
            load_reference_bundle(self.path)

    def test_fai_offsets_validated_even_with_correct_checksum(self):
        fai = self.root / "reference.fasta.fai"
        previous = compute_sha256(fai)
        fai.write_text("chrA\t320\t7\t40\t41\n")
        self.path.write_text(self.path.read_text().replace(previous, compute_sha256(fai)))
        with self.assertRaisesRegex(ReferenceError, "FAI"):
            load_reference_bundle(self.path)

    def test_annotation_assembly_and_intervals(self):
        original = self.path.read_text()
        before, after = original.rsplit('assembly_id = "synthetic-v1"', 1)
        self.path.write_text(before + 'assembly_id = "wrong"' + after)
        with self.assertRaisesRegex(ReferenceError, "assembly"):
            load_reference_bundle(self.path)
        annotation = self.root / "genes.gff3"
        previous = compute_sha256(annotation)
        annotation.write_text("chrUnknown\tx\tgene\t1\t4\t.\t+\t.\tID=x\n")
        self.path.write_text(original.replace(previous, compute_sha256(annotation)))
        with self.assertRaisesRegex(ReferenceError, "interval"):
            load_reference_bundle(self.path)

    def test_unsafe_asset_path(self):
        self.path.write_text(
            self.path.read_text().replace('path = "reference.fasta"', 'path = "../reference.fasta"')
        )
        with self.assertRaisesRegex(ReferenceError, "relative"):
            load_reference_bundle(self.path)


if __name__ == "__main__":
    unittest.main()
