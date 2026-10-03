"""Independent LOCO algebra and fail-closed statistical evidence contracts."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from adzuki_gwas_analysis.gwas_validation import (
    SCENARIOS,
    compare_results,
    simulate,
    wilson_interval,
)
from adzuki_gwas_analysis.kinship import kinship_for_chromosome
from adzuki_gwas_analysis.loader import compute_sha256
from adzuki_gwas_analysis.provenance import validate_provenance
from adzuki_gwas_analysis.tables import write_tsv


class LocoTests(unittest.TestCase):
    def test_excludes_whole_test_chromosome(self):
        g = np.array([[0, 1, 0, 2], [1, 2, 1, 1], [2, 0, 2, 0]], dtype=float)
        p = g.mean(axis=0) / 2
        k, count = kinship_for_chromosome(g, p, ["a", "a", "b", "c"], "a")
        centered = g[:, 2:] - g[:, 2:].mean(axis=0)
        np.testing.assert_allclose(k, centered @ centered.T / sum(2 * p[2:] * (1 - p[2:])))
        changed = g.copy()
        changed[:, :2] = 2 - changed[:, :2]
        k2, _ = kinship_for_chromosome(changed, changed.mean(axis=0) / 2, ["a", "a", "b", "c"], "a")
        np.testing.assert_allclose(k, k2)
        self.assertEqual(count, 2)

    def test_no_global_fallback(self):
        g = np.array([[0, 1], [1, 2], [2, 0]], dtype=float)
        with self.assertRaisesRegex(ValueError, "other chromosomes"):
            kinship_for_chromosome(g, g.mean(axis=0) / 2, ["a", "a"], "a")


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_simulation_runs_all_prespecified_scenarios_without_certification(self):
        plan = dict(
            schema_version=1,
            seed=483,
            replicates=3,
            n_samples=16,
            n_markers=15,
            alpha=0.05,
            max_null_upper95=0.2,
            min_power_lower95=0.5,
            kinship_mode="loco",
            scenarios=list(SCENARIOS),
        )
        path = self.root / "plan.json"
        path.write_text(json.dumps(plan))
        output = self.root / "simulation"
        simulate(path, output)
        report = json.loads((output / "calibration.json").read_text())
        self.assertEqual(report["commercial_status"], "not_validated")
        self.assertEqual(len(report["scenarios"]), 6)
        self.assertFalse(report["simulation_criteria_met"])
        self.assertIsNotNone(validate_provenance(output, required=True))
        self.assertGreater(wilson_interval(0, 3)[1], 0.5)

    def comparison(self):
        columns = (
            "chr",
            "pos",
            "ref",
            "alt",
            "effect_allele",
            "other_allele",
            "beta",
            "standard_error",
            "neg_log10_pvalue",
        )
        row = dict(zip(columns, ("a", 1, "A", "C", "C", "A", 0.2, 0.1, 1.4), strict=True))
        left, right = self.root / "left.tsv", self.root / "right.tsv"
        for path in (left, right):
            write_tsv(path, columns, [row])
        contract = dict(
            dataset_id="synthetic",
            assembly_id="v1",
            input_genotypes_sha256="a" * 64,
            input_phenotypes_sha256="b" * 64,
            covariates_sha256="c" * 64,
            model="quantitative-LMM",
            covariance_strategy="fixed-null-global",
            test="two-sided-t",
            effect_encoding="ALT_0_1_2",
            engine="internal",
            version="1",
            license="MIT",
            data_scope="synthetic",
        )
        plan = dict(
            schema_version=1,
            left_sha256=compute_sha256(left),
            right_sha256=compute_sha256(right),
            left_contract=contract,
            right_contract={**contract, "engine": "independent-fixture"},
            beta_atol=1e-7,
            se_atol=1e-7,
            logp_atol=1e-7,
        )
        path = self.root / "compare.json"
        path.write_text(json.dumps(plan))
        return left, right, path, plan

    def test_same_contract_comparison_and_alternate_test_are_distinct(self):
        left, right, path, plan = self.comparison()
        compare_results(left, right, path, self.root / "same")
        self.assertEqual(
            json.loads((self.root / "same/comparison.json").read_text())["status"], "passed"
        )
        plan["right_contract"]["test"] = "LRT"
        path.write_text(json.dumps(plan))
        compare_results(left, right, path, self.root / "different")
        result = json.loads((self.root / "different/comparison.json").read_text())
        self.assertEqual(result["status"], "not_comparable")
        self.assertEqual(result["commercial_status"], "not_validated")

    def test_modified_input_rejected_before_output(self):
        left, right, path, _ = self.comparison()
        right.write_text(right.read_text().replace("0.2", "0.8"))
        with self.assertRaisesRegex(ValueError, "predeclared plan"):
            compare_results(left, right, path, self.root / "bad")
        self.assertFalse((self.root / "bad").exists())

    def test_same_engine_different_version_is_not_independent_implementation(self):
        left, right, path, plan = self.comparison()
        plan["right_contract"].update(engine="INTERNAL", version="2")
        path.write_text(json.dumps(plan))
        compare_results(left, right, path, self.root / "same-engine")
        report = json.loads((self.root / "same-engine/comparison.json").read_text())
        self.assertFalse(report["independent_engine"])
        self.assertEqual(report["status"], "failed")
