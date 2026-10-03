"""Small synthetic references and a published Primer3 thermodynamic calculation."""

import importlib.util
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from adzuki_gwas_analysis.arms import PRIMER_FIELDS, PrimerSettings, design_candidate
from adzuki_gwas_analysis.assay_review import run_assay_review
from adzuki_gwas_analysis.primer_specificity import (
    PASS,
    assess_design,
    assess_thermodynamics,
    find_hits,
    load_settings,
    run_specificity,
)
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    validate_provenance,
)
from adzuki_gwas_analysis.reference import load_reference_bundle, reverse_complement
from adzuki_gwas_analysis.tables import read_tsv, write_tsv
from tests.reference_support import write_bundle


class SpecificityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.settings = json.loads(
            (Path(__file__).parents[1] / "config/specificity.example.json").read_text()
        )
        self.settings["thermodynamics"]["enabled"] = False
        self.settings_path = self.root / "settings.json"
        self.settings_path.write_text(json.dumps(self.settings))
        rng = random.Random(83)
        self.sequence = "".join(rng.choice("ACGT") for _ in range(500))
        self.bundle_path = write_bundle(self.root / "reference", sequence=self.sequence)
        ref = self.sequence[200]
        alt = next(base for base in "ACGT" if base != ref)
        self.design = {key: "" for key in PRIMER_FIELDS}
        self.design.update(
            candidate_id="synthetic-candidate",
            dataset_id="synthetic_trait",
            reference="Synthetic",
            design_id="synthetic-design",
            primer_version="synthetic-version",
            assembly_id="synthetic-v1",
            chr="chrA",
            pos="201",
            ref=ref,
            alt=alt,
            strand="+",
            specific_start="182",
            specific_end="201",
            common_start="251",
            common_end="270",
            product_bp="89",
            cohort_binding_check="unavailable",
            specific_ref_5to3=self.sequence[181:201],
            specific_alt_5to3=self.sequence[181:200] + alt,
            common_5to3=reverse_complement(self.sequence[250:270]),
        )

    def test_both_reactions_have_one_intended_product_on_declared_haplotypes(self):
        _, products, states = assess_design(
            self.design, load_reference_bundle(self.bundle_path), self.settings
        )
        self.assertEqual([s["specificity"] for s in states], [PASS, PASS])
        self.assertEqual([p["product_bp"] for p in products], [89, 89])
        self.assertTrue(all(s["assay_validation"] == "not_assessed" for s in states))
        self.assertTrue(all(s["data_scope"] == "synthetic" for s in states))
        self.assertTrue(all(s["reference_completeness"] == "not_assessed" for s in states))
        self.assertTrue(
            all(s["reference_search_scope"] == "supplied_reference_contigs_only" for s in states)
        )
        inconsistent = {**self.design, "strand": "-"}
        with self.assertRaisesRegex(ValueError, "3-prime coordinate"):
            assess_design(inconsistent, load_reference_bundle(self.bundle_path), self.settings)

    def test_generated_arms_candidates_match_on_both_strands(self):
        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                candidate = {
                    **self.design,
                    "strand": strand,
                    "interval_start": "1",
                    "interval_end": "500",
                    "snp_offset_0based": "200" if strand == "+" else "299",
                }
                sequence = self.sequence if strand == "+" else reverse_complement(self.sequence)
                designs, _ = design_candidate(
                    candidate,
                    sequence,
                    [],
                    PrimerSettings(
                        min_length=20,
                        max_length=20,
                        min_product_bp=75,
                        max_product_bp=90,
                        min_tm_C=0,
                        max_tm_C=100,
                        max_tm_difference_C=100,
                        min_gc_fraction=0,
                        max_gc_fraction=1,
                        max_homopolymer=20,
                        max_self_match=20,
                        max_pair_match=20,
                    ),
                    cohort_available=False,
                    neighbor_window=0,
                )
                self.assertTrue(designs)
                for design in designs:
                    _, _, states = assess_design(
                        {key: str(value) for key, value in design.items()},
                        load_reference_bundle(self.bundle_path),
                        self.settings,
                    )
                    self.assertEqual([s["specificity"] for s in states], [PASS, PASS])

    def test_inconsistent_inner_endpoint_and_allele_are_rejected(self):
        for changes in (
            {"specific_end": "202"},
            {"common_start": "252"},
            {"specific_alt_5to3": self.design["specific_ref_5to3"]},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                assess_design(
                    {**self.design, **changes},
                    load_reference_bundle(self.bundle_path),
                    self.settings,
                )

    def test_paralog_and_reverse_copy_are_off_target_products(self):
        for i, copy in enumerate((self.sequence, reverse_complement(self.sequence))):
            with self.subTest(i=i):
                path = write_bundle(
                    self.root / f"repeat{i}", sequence=self.sequence + "G" * 20 + copy
                )
                _, _, states = assess_design(
                    self.design, load_reference_bundle(path), self.settings
                )
                self.assertEqual(states[0]["specificity"], "failed")
                self.assertGreater(states[0]["n_off_target_products"], 0)

    def test_ambiguous_reference_is_not_a_pass(self):
        path = write_bundle(self.root / "ambiguous", sequence=self.sequence + "N")
        _, _, states = assess_design(self.design, load_reference_bundle(path), self.settings)
        self.assertTrue(all(s["specificity"] == "ambiguous" for s in states))

    def test_same_primer_can_generate_an_off_target_amplicon(self):
        common = self.design["common_5to3"]
        sequence = self.sequence + "G" * 600 + common + "T" * 60 + reverse_complement(common)
        path = write_bundle(self.root / "single-primer-product", sequence=sequence)
        _, products, states = assess_design(self.design, load_reference_bundle(path), self.settings)
        self.assertTrue(any(p["primer1_role"] == p["primer2_role"] == "common" for p in products))
        self.assertTrue(all(s["specificity"] == "failed" for s in states))

    def test_terminal_mismatch_and_reverse_coordinates(self):
        primer = "ACGTTCAGGATC"
        sequence = "GG" + primer + "CCC" + reverse_complement(primer) + "AA"
        hits = find_hits(sequence, primer, chrom="c", settings=self.settings)
        self.assertIn((3, "+"), [(h["start"], h["strand"]) for h in hits])
        self.assertIn((18, "-"), [(h["start"], h["strand"]) for h in hits])
        mismatch = primer[:-1] + ("A" if primer[-1] != "A" else "C")
        self.assertFalse(find_hits(mismatch, primer, chrom="c", settings=self.settings))
        interior = ("T" if primer[0] != "T" else "A") + primer[1:]
        self.assertEqual(
            find_hits(interior, primer, chrom="c", settings=self.settings)[0]["mismatches"], 1
        )

    def test_search_limit_prevents_partial_claims(self):
        settings = {**self.settings, "max_hits_per_primer": 2}
        with self.assertRaisesRegex(ValueError, "hit limit"):
            find_hits("A" * 60, "A" * 20, chrom="repeat", settings=settings)

    def test_candidate_window_budget_counts_nonqualifying_hits(self):
        # An exact terminal seed can be frequent even when no full primer qualifies.
        settings = {**self.settings, "max_candidate_windows": 3, "max_mismatches": 0}
        with self.assertRaisesRegex(ValueError, "candidate window budget"):
            find_hits("A" * 60, "C" * 17 + "AAA", chrom="repeat", settings=settings)

    def test_disabled_and_missing_thermodynamics_are_never_passes(self):
        with patch(
            "adzuki_gwas_analysis.primer_specificity.importlib.import_module",
            side_effect=ImportError,
        ):
            self.assertEqual(
                assess_thermodynamics(self.design, self.settings)["status"], "not_checked"
            )
            self.settings["thermodynamics"]["enabled"] = True
            with self.assertRaisesRegex(ValueError, "locked thermodynamics extra"):
                assess_thermodynamics(self.design, self.settings)

    @unittest.skipUnless(
        importlib.util.find_spec("primer3"), "optional thermodynamics extra not installed"
    )
    def test_published_primer3_hairpin_example_and_threshold(self):
        # Primer3-py 2.3.1 official docs: https://libnano.github.io/primer3-py/
        # Published hairpin example reports Tm=66.75 C with the library default conditions.
        design = {**self.design, "specific_ref_5to3": "CCCCCATCCGATCAGGGGG"}
        self.settings["thermodynamics"]["enabled"] = True
        result = assess_thermodynamics(design, self.settings)
        self.assertAlmostEqual(result["metrics"]["ref_hairpin"]["tm"], 66.75, delta=0.01)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["license"], "GPL-2.0-or-later")
        self.assertEqual(result["assay_validation"], "not_assessed")

    def test_bundle_provenance_and_failed_configuration(self):
        context_dir = self.root / "context"
        context_dir.mkdir()
        finish_provenance(
            context_dir,
            operation="candidate_context",
            inputs={},
            initial={},
            parameters={"reference": load_reference_bundle(self.bundle_path).metadata()},
            environment=generation_environment(),
            families=[],
        )
        context_run = validate_provenance(context_dir, required=True)
        design_dir = self.root / "design"
        design_dir.mkdir()
        write_tsv(design_dir / "primer_candidates.tsv", PRIMER_FIELDS, [self.design])
        write_tsv(
            design_dir / "arms_marker_candidates.tsv",
            ("candidate_id", "design_count", "rejection_reasons"),
            [
                {
                    "candidate_id": self.design["candidate_id"],
                    "design_count": 1,
                    "rejection_reasons": "",
                }
            ],
        )
        finish_provenance(
            design_dir,
            operation="arms_design",
            inputs={},
            initial={},
            parameters={"context_run_id": context_run["run_id"]},
            environment=generation_environment(),
            families=[],
        )
        output = self.root / "review"
        run_specificity(
            design_dir=design_dir,
            context_dir=context_dir,
            bundle_path=self.bundle_path,
            settings_path=self.settings_path,
            output_dir=output,
        )
        self.assertIsNotNone(validate_provenance(output, required=True))
        self.assertEqual(len(list(read_tsv(output / "specificity_review.tsv"))), 2)
        handoff = self.root / "handoff"
        run_assay_review(design_dir=design_dir, output_dir=handoff, specificity_dir=output)
        self.assertTrue((handoff / "primer_confirmation/specificity_review.tsv").is_file())
        self.assertEqual(
            list(read_tsv(handoff / "marker_states.tsv"))[0]["state"], "computational_candidate"
        )
        self.settings["max_reference_bases"] = 10
        self.settings_path.write_text(json.dumps(self.settings))
        with self.assertRaisesRegex(ValueError, "operating limit"):
            run_specificity(
                design_dir=design_dir,
                context_dir=context_dir,
                bundle_path=self.bundle_path,
                settings_path=self.settings_path,
                output_dir=self.root / "partial",
            )
        self.assertFalse((self.root / "partial").exists())
        self.settings["unknown"] = True
        self.settings_path.write_text(json.dumps(self.settings))
        with self.assertRaises(ValueError):
            load_settings(self.settings_path)
