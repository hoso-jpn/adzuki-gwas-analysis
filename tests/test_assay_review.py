"""Synthetic assay handoff/results, denominator checks and guarded state transitions."""

import json
import random
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.arms import run_arms
from adzuki_gwas_analysis.assay_contract import IDENTITY_FIELDS, RESULT_FIELDS
from adzuki_gwas_analysis.assay_metrics import assay_metrics, review_state
from adzuki_gwas_analysis.assay_review import HISTORY_FIELDS, run_assay_review
from adzuki_gwas_analysis.context import run_context
from adzuki_gwas_analysis.provenance import validate_provenance
from adzuki_gwas_analysis.tables import read_tsv, write_tsv
from tests.reference_support import write_bundle


class AssayReviewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        rng = random.Random(1)
        sequence = "".join(rng.choice("ACGT") for _ in range(500))
        bundle = write_bundle(self.root / "reference", sequence=sequence)
        ref, alt = sequence[200], next(base for base in "ACGT" if base != sequence[200])
        candidates = self.root / "candidates.tsv"
        candidates.write_text(
            f"candidate_id\tchr\tpos\tref\talt\ttrait\nc1\tchrA\t201\t{ref}\t{alt}\ttrait\n"
        )
        context = self.root / "context"
        run_context(
            candidate_table=candidates,
            bundle_path=bundle,
            dataset_id="synthetic_trait",
            output_dir=context,
            flank_bp=180,
            neighbor_window_bp=180,
        )
        config = self.root / "primer.toml"
        config.write_text(
            "schema_version = 1\n[primer]\nmin_length = 18\nmax_length = 20\n"
            "min_product_bp = 75\nmax_product_bp = 90\nmin_tm_C = 35\n"
            "max_tm_C = 80\nmax_tm_difference_C = 12\n"
        )
        self.design_dir = self.root / "design"
        run_arms(context_dir=context, config=config, output_dir=self.design_dir)
        self.designs = list(read_tsv(self.design_dir / "primer_candidates.tsv"))
        self.design = self.designs[0]
        self.results = self.root / "results.tsv"
        self.metadata = self.root / "metadata.toml"
        self.meta = dict(
            schema_version=1,
            validation_dataset_id="synthetic-assay-1",
            conditions_id="conditions-1",
            measurement_conditions="synthetic PCR fixture; no physical experiment",
            target_population="synthetic population",
            reference_method="synthetic known genotype",
            analyst_review_reference="synthetic-review-1",
            reviewer_id="fixture-reviewer",
            population_evidence_reference="synthetic-population-evidence",
            minimum_samples=3,
            minimum_comparable_samples=3,
            minimum_replicated_samples=1,
            minimum_call_rate=1.0,
            minimum_concordance=1.0,
            data_scope="synthetic",
            requested_state="assay_validated",
            analyst_approved=True,
            population_scope_confirmed=True,
        )
        self.rows = [
            self.row("001", "0", truth="0"),
            self.row("NA", "1", truth="1"),
            self.row("s3", "2", truth="2"),
            self.row("001", "0", truth="0", replicate="2"),
            *(
                self.row("positive-" + call, call, control="positive", expected=call)
                for call in ("0", "1", "2")
            ),
            self.row("negative", "NA", control="negative"),
        ]
        self.write_inputs()

    def row(self, sample, call, *, truth="NA", control="sample", expected="NA", replicate="1"):
        return {
            **{key: self.design[key] for key in IDENTITY_FIELDS},
            "allele_encoding": "ALT_0_1_2_forward",
            "conditions_id": "conditions-1",
            "plate_id": "plate-1",
            "batch_id": "batch-1",
            "sample_id": sample,
            "replicate_id": replicate,
            "control_type": control,
            "call": call,
            "no_call_reason": "",
            "expected_genotype": expected,
            "reference_genotype": truth,
        }

    def write_inputs(self, **overrides):
        values = {**self.meta, **overrides}
        self.metadata.write_text(
            "\n".join(f"{key} = {json.dumps(value)}" for key, value in values.items())
        )
        write_tsv(self.results, RESULT_FIELDS, self.rows)

    def run_review(self, output, **overrides):
        values = dict(
            design_dir=self.design_dir,
            output_dir=output,
            results=self.results,
            metadata_path=self.metadata,
        )
        values.update(overrides)
        run_assay_review(**values)
        return {row["design_id"]: row for row in read_tsv(output / "marker_states.tsv")}

    def test_handoff_without_measurements_cannot_be_validated(self):
        output = self.root / "handoff"
        states = self.run_review(output, results=None, metadata_path=None)
        self.assertTrue(all(row["state"] == "computational_candidate" for row in states.values()))
        self.assertFalse(list(read_tsv(output / "assay_return_template.tsv")))
        self.assertFalse(list(read_tsv(output / "population_validated_panel.tsv")))
        self.assertEqual(len(list(read_tsv(output / "handoff.tsv"))), 3)
        self.assertIsNotNone(validate_provenance(output, required=True))

    def test_synthetic_evidence_and_population_transition_are_labeled(self):
        first = self.root / "assay"
        states = self.run_review(first)
        state = states[self.design["design_id"]]
        self.assertEqual(
            (state["state"], state["simulated_state"]),
            ("computational_candidate", "assay_validated"),
        )
        self.assertIn("synthetic_only", state["reasons"])
        self.write_inputs(
            requested_state="population_validated", validation_dataset_id="synthetic-population-2"
        )
        second = self.root / "population"
        states = self.run_review(second, previous_review=first)
        self.assertEqual(
            states[self.design["design_id"]]["simulated_state"], "population_validated"
        )
        self.assertFalse(list(read_tsv(second / "population_validated_panel.tsv")))
        self.assertIsNotNone(validate_provenance(second, required=True))

    def test_no_calls_concordance_and_opaque_sample_ids(self):
        self.rows[2].update(call="NA", no_call_reason="weak_band")
        self.write_inputs()
        output = self.root / "no-call"
        states = self.run_review(output)
        values = next(read_tsv(output / "assay_metrics.tsv"))
        self.assertEqual(
            (values["n_sample_attempts"], values["n_called"], values["n_no_calls"]), ("4", "3", "1")
        )
        self.assertEqual(float(values["call_rate"]), 0.75)
        self.assertEqual((values["n_comparable_calls"], float(values["concordance"])), ("3", 1.0))
        self.assertEqual(
            states[self.design["design_id"]]["simulated_state"], "computational_candidate"
        )
        returned = list(read_tsv(output / "assay_results.tsv"))
        self.assertEqual([row["sample_id"] for row in returned[:2]], ["001", "NA"])

    def test_duplicate_old_primer_and_reversed_alleles_are_rejected(self):
        original = [row.copy() for row in self.rows]
        cases = [original + [original[0].copy()]]
        for field, value in (
            ("primer_version", "old-version"),
            ("allele_encoding", "REF_0_1_2"),
            ("ref", self.design["alt"]),
        ):
            modified = [row.copy() for row in original]
            modified[0][field] = value
            cases.append(modified)
        for number, rows in enumerate(cases):
            with self.subTest(number=number):
                self.rows = rows
                self.write_inputs()
                output = self.root / f"invalid-{number}"
                with self.assertRaises(ValueError):
                    self.run_review(output)
                self.assertFalse(output.exists())

    def test_failed_controls_and_replicate_disagreement_block_acceptance(self):
        self.rows[-1]["call"] = "0"
        self.rows[3]["call"] = "1"
        self.write_inputs()
        output = self.root / "failed-controls"
        state = self.run_review(output)[self.design["design_id"]]
        metrics = next(read_tsv(output / "assay_metrics.tsv"))
        self.assertEqual(
            (metrics["n_failed_controls"], metrics["n_replicate_disagreements"]), ("1", "1")
        )
        self.assertEqual(state["simulated_state"], "computational_candidate")
        self.assertIn("n_failed_controls", state["reasons"])
        self.assertTrue(list(read_tsv(output / "alternative_designs.tsv")))

    def test_missing_controls_and_absent_reference_genotypes_are_not_passes(self):
        self.rows = self.rows[:4]
        for row in self.rows:
            row["reference_genotype"] = "NA"
        values = assay_metrics(self.rows)
        self.assertIsNone(values["concordance"])
        self.assertEqual(values["n_incomplete_control_plates"], 1)
        state = review_state(values, self.meta, None)
        self.assertEqual(state["simulated_state"], "computational_candidate")

    def test_review_and_prior_evidence_are_required_for_promotion(self):
        values = assay_metrics(self.rows)
        # Exercise operational rules in memory; every persisted test result remains synthetic.
        actual = {**self.meta, "data_scope": "public"}
        self.assertEqual(review_state(values, actual, None)["state"], "assay_validated")
        self.assertEqual(
            review_state(values, {**actual, "analyst_approved": False}, None)["state"],
            "computational_candidate",
        )
        requested = {**actual, "requested_state": "population_validated"}
        self.assertIn(
            "prior_assay_validation_required", review_state(values, requested, None)["reasons"]
        )
        prior = {
            "state": "assay_validated",
            "data_scope": "public",
            "validation_dataset_id": "previous-dataset",
        }
        self.assertEqual(review_state(values, requested, prior)["state"], "population_validated")
        self.assertIn(
            "distinct_population",
            review_state(
                values,
                requested,
                {**prior, "validation_dataset_id": actual["validation_dataset_id"]},
            )["reasons"],
        )

    def test_redesign_history_binds_prior_and_replacement_versions(self):
        first = self.root / "first"
        self.run_review(first)
        history = self.root / "history.tsv"
        write_tsv(
            history,
            HISTORY_FIELDS,
            [
                {
                    "candidate_id": self.design["candidate_id"],
                    "prior_design_id": self.design["design_id"],
                    "replacement_design_id": self.designs[1]["design_id"],
                    "reason": "synthetic redesign exercise",
                    "evidence_reference": "synthetic assay",
                    "analyst_id": "fixture-reviewer",
                }
            ],
        )
        second = self.root / "replacement"
        self.run_review(second, previous_review=first, redesign_history=history)
        self.assertEqual(len(list(read_tsv(second / "redesign_history.tsv"))), 1)
        (first / "marker_states.tsv").write_text("tampered\n")
        with self.assertRaisesRegex(ValueError, "checksum"):
            self.run_review(self.root / "bad-prior", previous_review=first)
