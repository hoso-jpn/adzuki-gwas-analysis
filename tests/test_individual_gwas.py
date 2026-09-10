"""Synthetic individual data and independent regression checks; no customer fixtures."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import stats

from adzuki_gwas_analysis.individual_gwas import run_individual_gwas
from adzuki_gwas_analysis.individual_inputs import load_config, load_individual_inputs
from adzuki_gwas_analysis.mixed_model import fit_null, test_marker
from adzuki_gwas_analysis.provenance import validate_provenance
from adzuki_gwas_analysis.reference import load_reference_bundle
from adzuki_gwas_analysis.tables import read_tsv, write_tsv
from tests.reference_support import write_bundle


class MixedModelTests(unittest.TestCase):
    def test_identity_kinship_matches_independent_ordinary_regression(self):
        rng = np.random.default_rng(913)
        x = rng.binomial(2, 0.3, 40).astype(float)
        y = 1.7 * x + rng.normal(size=40)
        model = fit_null(y, np.ones((40, 1)), np.eye(40))
        beta, se, logp = test_marker(model, x)
        expected = stats.linregress(x, y)
        self.assertAlmostEqual(beta, expected.slope, places=10)
        self.assertAlmostEqual(se, expected.stderr, places=10)
        self.assertAlmostEqual(logp, -np.log10(expected.pvalue), places=9)

    def test_covariance_weighted_test_matches_direct_gls(self):
        rng = np.random.default_rng(131)
        z = rng.normal(size=(32, 10))
        k = z @ z.T / 10
        c = np.column_stack((np.ones(32), rng.normal(size=32)))
        x, y = rng.normal(size=(2, 32))
        model = fit_null(y, c, k)
        beta, se, logp = test_marker(model, x)
        design = np.column_stack((c, x))
        inverse = np.linalg.inv(k + model.delta * np.eye(32))
        information_inverse = np.linalg.inv(design.T @ inverse @ design)
        coefficients = information_inverse @ design.T @ inverse @ y
        residual = y - design @ coefficients
        expected_se = np.sqrt((residual @ inverse @ residual) / 29 * information_inverse[-1, -1])
        self.assertAlmostEqual(beta, coefficients[-1], places=9)
        self.assertAlmostEqual(se, expected_se, places=9)
        self.assertAlmostEqual(logp, -np.log10(2 * stats.t.sf(abs(beta / se), 29)), places=9)

    def test_rank_deficiency_and_invalid_kinship_fail(self):
        with self.assertRaisesRegex(ValueError, "rank deficient"):
            fit_null(np.arange(8, dtype=float), np.ones((8, 2)), np.eye(8))
        with self.assertRaisesRegex(ValueError, "positive semidefinite"):
            fit_null(np.arange(8, dtype=float), np.ones((8, 1)), -np.eye(8))


class IndividualGWASTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bundle_path = write_bundle(self.root / "reference")
        self.bundle = load_reference_bundle(self.bundle_path)
        self.config_path = self.root / "run.toml"
        self.config_values = dict(
            schema_version=1,
            dataset_id="synthetic_trait",
            cohort_id="synthetic-cohort",
            analysis_id="synthetic-analysis",
            reference="Synthetic",
            assembly_id="synthetic-v1",
            species="synthetic plant",
            trait="mass",
            trait_unit="g",
            trait_coding="untransformed mass",
            trait_type="quantitative",
            design_note="synthetic one observation per sample",
            genotype_encoding="ALT_0_1_2",
            data_scope="synthetic",
            data_use_confirmed=True,
            permission_reference="synthetic",
            transfer_method="generated in temporary directory",
            retention_policy="test duration",
            deletion_policy="TemporaryDirectory cleanup",
            output_ownership="synthetic fixture",
            clustering_distance=10,
            covariates=["block_numeric"],
        )
        self.write_config()
        rng = np.random.default_rng(819)
        self.ids = ["001", "NA", *(f"sample-{i}" for i in range(62))]
        self.genotype = rng.binomial(2, 0.35, size=(64, 26)).astype(float)
        self.genotype[:, -1] = 0
        self.genotype[3, 2] = np.nan
        self.y = 3 * self.genotype[:, 0] + rng.normal(size=64)
        self.paths = {
            name: self.root / (name + ".tsv")
            for name in ("genotypes", "phenotypes", "samples", "variants")
        }
        markers = [f"m{i}" for i in range(26)]
        write_tsv(
            self.paths["genotypes"],
            ("sample_id", *markers),
            [
                {
                    "sample_id": identifier,
                    **{
                        marker: "NA" if np.isnan(self.genotype[i, j]) else int(self.genotype[i, j])
                        for j, marker in enumerate(markers)
                    },
                }
                for i, identifier in enumerate(self.ids)
            ],
        )
        write_tsv(
            self.paths["phenotypes"],
            ("sample_id", "trait", "value", "unit"),
            [
                {"sample_id": identifier, "trait": "mass", "value": self.y[i], "unit": "g"}
                for i, identifier in reversed(list(enumerate(self.ids)))
            ],
        )
        write_tsv(
            self.paths["samples"],
            ("sample_id", "cohort_id", "block_numeric"),
            [
                {"sample_id": identifier, "cohort_id": "synthetic-cohort", "block_numeric": i % 3}
                for i, identifier in enumerate(self.ids)
            ],
        )
        write_tsv(
            self.paths["variants"],
            ("marker_id", "chr", "pos", "ref", "alt", "assembly_id"),
            [
                {
                    "marker_id": marker,
                    "chr": "chrA",
                    "pos": i * 4 + 1,
                    "ref": "A",
                    "alt": "C",
                    "assembly_id": "synthetic-v1",
                }
                for i, marker in enumerate(markers)
            ],
        )

    def write_config(self, **overrides):
        values = {**self.config_values, **overrides}
        self.config_path.write_text("\n".join(f"{k} = {json.dumps(v)}" for k, v in values.items()))

    def load(self):
        return load_individual_inputs(
            **self.paths, bundle=self.bundle, config=load_config(self.config_path, self.bundle)
        )

    def test_identity_alignment_and_missing_dosage(self):
        data = self.load()
        self.assertEqual(data.sample_ids[:2], ["001", "NA"])
        np.testing.assert_array_equal(data.phenotype, self.y)
        self.assertTrue(np.isnan(data.dosage[3, 2]))

    def test_duplicate_missing_sample_and_mixed_units_fail(self):
        path = self.paths["phenotypes"]
        original = path.read_text()
        for malformed in (
            original + original.splitlines()[1] + "\n",
            "\n".join(original.splitlines()[:-1]) + "\n",
            original.replace("\tg\n", "\tkg\n", 1),
        ):
            with self.subTest():
                path.write_text(malformed)
                with self.assertRaises(ValueError):
                    self.load()
        path.write_text(original)

    def test_wrong_reference_dosage_and_memory_limit_fail(self):
        variants = self.paths["variants"]
        original = variants.read_text()
        variants.write_text(original.replace("synthetic-v1", "wrong-assembly", 1))
        with self.assertRaisesRegex(ValueError, "assembly"):
            self.load()
        variants.write_text(original)
        genotypes = self.paths["genotypes"]
        original = genotypes.read_text()
        genotypes.write_text(original.replace("\tNA\t", "\t-9\t", 1))
        with self.assertRaisesRegex(ValueError, "0/1/2/NA"):
            self.load()
        genotypes.write_text(original)
        self.write_config(max_genotype_cells=5)
        with self.assertRaisesRegex(ValueError, "cell limit"):
            self.load()

    def test_synthetic_e2e_qc_plots_accounting_and_reproducibility(self):
        outputs = [self.root / "run1", self.root / "run2"]
        for output in outputs:
            run_individual_gwas(
                **self.paths,
                bundle_path=self.bundle_path,
                config_path=self.config_path,
                output_dir=output,
            )
            self.assertIsNotNone(validate_provenance(output, required=True))
        self.assertEqual(
            (outputs[0] / "association_results.tsv").read_bytes(),
            (outputs[1] / "association_results.tsv").read_bytes(),
        )
        result = list(read_tsv(outputs[0] / "normalized_summary.tsv"))
        self.assertEqual(len(result), 25)
        lead = max(result, key=lambda row: float(row["neg_log10_pvalue"]))
        self.assertEqual(lead["pos"], "1")
        self.assertGreater(float(lead["beta"]), 0)
        self.assertEqual(lead["effect_allele"], "C")
        self.assertTrue(list(read_tsv(outputs[0] / "candidate_snps.tsv")))
        diag = json.loads((outputs[0] / "statistical_diagnostics.json").read_text())
        self.assertEqual((diag["n_input"], diag["n_tests"], diag["n_excluded"]), (26, 25, 1))
        self.assertEqual(diag["analysis_origin"], "individual_gwas")
        self.assertIn(
            "Individual-level quantitative GWAS was run",
            (outputs[0] / "analysis_report.md").read_text(),
        )
        self.assertTrue((outputs[0] / "manhattan.png").is_file())
        self.assertTrue((outputs[0] / "qq.png").is_file())

    def test_pca_and_covariate_contract(self):
        self.write_config(n_pcs=2)
        output = self.root / "with-pcs"
        run_individual_gwas(
            **self.paths,
            bundle_path=self.bundle_path,
            config_path=self.config_path,
            output_dir=output,
        )
        self.assertEqual(len(list(read_tsv(output / "pc_scores.tsv"))), 64)
        model = json.loads((output / "model.json").read_text())
        self.assertEqual(model["n_pcs"], 2)
        self.assertIn("block_numeric", model["covariate_transform"])
