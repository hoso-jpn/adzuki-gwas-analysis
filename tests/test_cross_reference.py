"""Synthetic reversed/gapped/duplicated alignments with complete candidate accounting."""

import json
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.chain_mapping import load_chain_metadata, map_chain
from adzuki_gwas_analysis.cross_reference import CANDIDATE_FIELDS, run_cross_reference
from adzuki_gwas_analysis.loader import compute_sha256
from adzuki_gwas_analysis.provenance import validate_provenance
from adzuki_gwas_analysis.reference import load_reference_bundle
from adzuki_gwas_analysis.tables import read_tsv, write_tsv
from tests.reference_support import write_bundle


class CrossReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source_sequence = "ACGT" * 10
        target_sequence = list("A" * 40)
        for pos, base in {1: "C", 7: "C", 9: "C", 15: "G", 30: "T", 35: "G", 40: "T"}.items():
            target_sequence[pos - 1] = base
        self.target_sequence = "".join(target_sequence)
        self.source_bundle = write_bundle(
            self.root / "source",
            sequence=self.source_sequence,
            assembly="source-v1",
            dataset="source-trait",
            reference="Source",
        )
        self.target_bundle = write_bundle(
            self.root / "target",
            sequence=self.target_sequence,
            assembly="target-v1",
            dataset="target-trait",
            reference="Target",
        )
        self.source_table, self.target_table = self.root / "source.tsv", self.root / "target.tsv"
        self.source_rows = [self.candidate("source", pos) for pos in (1, 5, 6, 8, 15, 18, 21, 40)]
        self.target_rows = [self.candidate("target", pos) for pos in (1, 7, 30, 40, 33)]
        self.target_rows[0]["beta"] = "-2"
        self.target_rows[2]["alt"] = self.target_rows[2]["effect_allele"] = "G"
        self.write_tables()
        self.chain = self.root / "alignment.chain"
        self.chain.write_text(
            "chain 100 chrA 40 + 0 12 chrA 40 + 0 11 1\n4 2 1\n6\n\n"
            "chain 80 chrA 40 + 20 30 chrA 40 - 10 20 2\n10\n\n"
            "chain 60 chrA 40 + 14 15 chrA 40 + 14 15 3\n1\n\n"
            "chain 50 chrA 40 + 34 35 chrA 40 + 14 15 4\n1\n\n"
            "chain 40 chrA 40 + 39 40 chrA 40 + 39 40 5\n1\n"
        )
        self.metadata = self.root / "chain.toml"
        self.write_metadata()

    def candidate(self, side, pos):
        sequence = self.source_sequence if side == "source" else self.target_sequence
        ref = sequence[pos - 1]
        alt = "C" if ref == "A" else "A"
        return dict(
            candidate_id=f"{side}-{pos}",
            dataset_id=f"{side}-trait",
            cohort_id="same-cohort",
            analysis_id=side + "-analysis",
            reference=side.title(),
            assembly_id=side + "-v1",
            chr="chrA",
            pos=str(pos),
            ref=ref,
            alt=alt,
            trait="trait",
            trait_unit="g",
            trait_coding="mass",
            test="wald",
            effect_scale="g",
            effect_allele=alt,
            other_allele=ref,
            effect_orientation="aligned_to_ALT",
            beta="2",
            se="0.2",
            neg_log10_pvalue="8" if side == "source" else "6",
        )

    def write_tables(self):
        write_tsv(self.source_table, CANDIDATE_FIELDS, self.source_rows)
        write_tsv(self.target_table, CANDIDATE_FIELDS, self.target_rows)

    def write_metadata(self, **overrides):
        values = dict(
            schema_version=1,
            source_assembly="source-v1",
            target_assembly="target-v1",
            source_uri="synthetic fixture",
            version="fixture-v1",
            license="CC0",
            validation_reference="constructed mini-assembly mapping",
            reviewer_id="synthetic-reviewer",
            chain_path="alignment.chain",
            chain_sha256=compute_sha256(self.chain),
            reviewed=True,
            direction="query_to_target",
        )
        values.update(overrides)
        self.metadata.write_text(
            "\n".join(f"{key} = {json.dumps(value)}" for key, value in values.items())
        )

    def run_comparison(self, output, **overrides):
        values = dict(
            source_table=self.source_table,
            source_bundle=self.source_bundle,
            target_table=self.target_table,
            target_bundle=self.target_bundle,
            trait="trait",
            output_dir=output,
            chain_metadata=self.metadata,
        )
        values.update(overrides)
        run_cross_reference(**values)

    def test_reversal_indels_ends_multimap_unmapped_and_accounting(self):
        output = self.root / "comparison"
        self.run_comparison(output)
        self.assertIsNotNone(validate_provenance(output, required=True))
        accounting = list(read_tsv(output / "candidate_accounting.tsv"))
        self.assertEqual(len(accounting), len(self.source_rows) + len(self.target_rows))
        statuses = {row["candidate_id"]: row["status"] for row in accounting}
        self.assertEqual(statuses["source-5"], "unmapped")
        self.assertEqual(statuses["source-18"], "unmapped")
        self.assertEqual(statuses["source-15"], "multimap")
        self.assertEqual(statuses["source-8"], "ambiguous")
        self.assertEqual(statuses["target-33"], "not_compared")
        comparisons = {
            row["source_candidate_id"]: row
            for row in read_tsv(output / "cross_reference_comparison.tsv")
        }
        self.assertEqual(set(comparisons), {"source-1", "source-6", "source-21", "source-40"})
        self.assertEqual(comparisons["source-6"]["target_pos"], "7")
        reverse = comparisons["source-21"]
        self.assertEqual((reverse["target_pos"], reverse["strand"]), ("30", "-"))
        self.assertEqual(comparisons["source-40"]["target_pos"], "40")
        swap = comparisons["source-1"]
        self.assertEqual(swap["allele_alignment"], "REF_ALT_swapped")
        self.assertEqual(swap["effect_direction_comparison"], "same_sign_after_alignment")
        self.assertEqual(swap["cohort_relationship"], "same_cohort")
        self.assertEqual(
            (swap["source_neg_log10_pvalue"], swap["target_neg_log10_pvalue"]), ("8", "6")
        )
        annotations = list(read_tsv(output / "candidate_annotations.tsv"))
        self.assertEqual(
            next(r for r in annotations if r["candidate_id"] == "source-21")["distance_bp"], "1"
        )

    def test_no_mapping_asset_runs_reference_local_annotation(self):
        output = self.root / "local"
        self.run_comparison(output, chain_metadata=None)
        summary = json.loads((output / "comparison_summary.json").read_text())
        self.assertEqual(summary["cross_reference_status"], "not_assessed")
        self.assertFalse(list(read_tsv(output / "cross_reference_comparison.tsv")))
        self.assertEqual(summary["source_status_counts"]["not_assessed"], 8)
        self.assertTrue(list(read_tsv(output / "candidate_annotations.tsv")))

    def test_assembly_chain_integrity_and_terminal_blocks(self):
        source, target = (
            load_reference_bundle(self.source_bundle),
            load_reference_bundle(self.target_bundle),
        )
        self.write_metadata(source_assembly="wrong")
        with self.assertRaisesRegex(ValueError, "assembly"):
            load_chain_metadata(self.metadata, source, target)
        self.write_metadata()
        self.chain.write_text(self.chain.read_text().replace("chrA 40 +", "chrA 41 +", 1))
        with self.assertRaisesRegex(ValueError, "checksum"):
            load_chain_metadata(self.metadata, source, target)
        self.write_metadata()
        with self.assertRaisesRegex(ValueError, "length mismatch"):
            map_chain(self.chain, source, target, self.source_rows)
        self.chain.write_text("chain 5 chrA 40 + 0 8 chrA 40 + 0 8 1\n4\n")
        with self.assertRaisesRegex(ValueError, "endpoints"):
            map_chain(self.chain, source, target, self.source_rows)

    def test_different_trait_and_reference_allele_fail_before_output(self):
        output = self.root / "bad"
        self.target_rows[0]["trait"] = "different"
        self.write_tables()
        with self.assertRaisesRegex(ValueError, "trait"):
            self.run_comparison(output)
        self.assertFalse(output.exists())
        self.target_rows[0]["trait"] = "trait"
        self.source_rows[0]["ref"] = "G"
        self.write_tables()
        with self.assertRaisesRegex(ValueError, "reference allele"):
            self.run_comparison(output)

    def test_annotation_bundle_assembly_mismatch_fails(self):
        text = self.target_bundle.read_text()
        before, annotation = text.split("[annotation]")
        self.target_bundle.write_text(
            before + "[annotation]" + annotation.replace("target-v1", "source-v1")
        )
        with self.assertRaisesRegex(ValueError, "assembly"):
            self.run_comparison(self.root / "bad-annotation")

    def test_unknown_effect_orientation_is_not_compared_as_a_sign(self):
        self.source_rows[0]["effect_orientation"] = "unresolved_strand"
        self.write_tables()
        output = self.root / "unknown-effect"
        self.run_comparison(output)
        comparison = next(read_tsv(output / "cross_reference_comparison.tsv"))
        self.assertEqual(
            comparison["effect_direction_comparison"], "not_assessed_unresolved_strand"
        )
