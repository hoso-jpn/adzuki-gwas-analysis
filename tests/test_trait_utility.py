"""Synthetic evidence contracts; no real trait performance is asserted by these tests."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from adzuki_gwas_analysis.assay_contract import IDENTITY_FIELDS
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    validate_provenance,
)
from adzuki_gwas_analysis.tables import read_tsv, write_tsv
from adzuki_gwas_analysis.trait_utility import (
    TEXT_FIELDS,
    _read_units,
    assess_evidence,
    run_trait_review,
)


class TraitEvidenceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.design = {key: key + "-synthetic" for key in IDENTITY_FIELDS}
        self.design.update(
            ref="A",
            alt="C",
            state="population_validated",
            data_scope="public",
            target_population="p1",
        )
        self.hashes = {
            key: character * 64
            for key, character in (("validation", "a"), ("discovery", "b"), ("analysis_run", "c"))
        }
        self.meta = {key: key + "-fixture" for key in TEXT_FIELDS}
        self.meta.update({key: self.design[key] for key in IDENTITY_FIELDS})
        self.meta.update(
            schema_version=1,
            data_scope="public",
            desired_direction="increase",
            effect_allele="C",
            analyst_approved=True,
            independence_basis="independent_families",
            minimum_samples=4,
            minimum_effect=0.1,
            minimum_utility=0.1,
            trait="mass",
            trait_unit="g",
            target_population="p1",
            environment="year1-site1",
            validation_data_sha256="a" * 64,
            discovery_data_sha256="b" * 64,
            analysis_run_sha256="c" * 64,
            association_effect=dict(estimate=2.0, lower95=1.0, upper95=3.0, unit="g"),
            selection_contrast=dict(estimate=1.0, lower95=0.5, upper95=1.5, unit="g"),
        )
        self.discovery = {"d1": dict(sample_id="d1", family_id="df1")}
        self.validation = {
            f"v{i}": dict(
                sample_id=f"v{i}",
                family_id=f"vf{i}",
                population="p1",
                environment="year1-site1",
                trait="mass",
                unit="g",
                value=str(i),
            )
            for i in range(4)
        }

    def assess(self, **overrides):
        return assess_evidence(
            {**self.meta, **overrides}, self.design, self.discovery, self.validation, self.hashes
        )

    def test_assay_accuracy_alone_never_becomes_trait_utility(self):
        source = self.root / "assay"
        source.mkdir()
        write_tsv(source / "marker_states.tsv", self.design.keys(), [self.design])
        finish_provenance(
            source,
            operation="assay_review",
            inputs={},
            initial={},
            parameters={},
            environment=generation_environment(),
            families=[],
        )
        target = self.root / "review"
        run_trait_review(assay_dir=source, output_dir=target)
        row = list(read_tsv(target / "evidence_axes.tsv"))[0]
        self.assertEqual(
            row["genotyping_validation"], "genotyping_validated_in_declared_population"
        )
        self.assertEqual(row["trait_utility"], "not_assessed")
        self.assertFalse(list(read_tsv(target / "trait_supported_panel.tsv")))
        self.assertIsNotNone(validate_provenance(target, required=True))

    def test_reviewed_numeric_evidence_and_synthetic_suppression(self):
        self.assertEqual(self.assess()["trait_utility"], "supported_within_declared_scope")
        result = self.assess(data_scope="synthetic")
        self.assertEqual(result["trait_utility"], "not_assessed")
        self.assertEqual(result["simulated_trait_utility"], "supported_within_declared_scope")

    def test_null_reversed_and_utility_failure(self):
        for effect in (
            dict(estimate=0.0, lower95=-1.0, upper95=1.0, unit="g"),
            dict(estimate=-2.0, lower95=-3.0, upper95=-1.0, unit="g"),
        ):
            with self.subTest(effect=effect):
                self.assertEqual(
                    self.assess(association_effect=effect)["trait_utility"], "not_supported"
                )
        result = self.assess(
            selection_contrast=dict(estimate=0.0, lower95=-0.5, upper95=0.5, unit="g")
        )
        self.assertEqual(result["association_replication"], "supported")
        self.assertEqual(result["trait_utility"], "not_supported")

    def test_reuse_scope_relatedness_and_unreviewed_fail(self):
        for override in (
            dict(validation_dataset_id=self.meta["discovery_dataset_id"]),
            dict(environment="different-year"),
            dict(analyst_approved=False),
        ):
            with self.subTest(override=override):
                result = self.assess(**override)
                self.assertEqual(result["trait_utility"], "not_supported")
                self.assertEqual(result["association_replication"], "not_supported")
        self.validation["v0"]["family_id"] = "df1"
        self.assertIn("independent_family_scope_not_met", self.assess()["reasons"])
        self.assertEqual(self.assess()["association_replication"], "not_supported")

    def test_identity_hash_unit_and_unknown_contract_rejected(self):
        cases = (
            dict(primer_version="different"),
            dict(validation_data_sha256="d" * 64),
            dict(effect_allele="A"),
            dict(extra=True),
            dict(association_effect=dict(estimate=1, lower95=2, upper95=3, unit="g")),
        )
        for override in cases:
            with self.subTest(override=override), self.assertRaises(ValueError):
                self.assess(**override)

    def test_duplicate_biological_units_are_not_independent_replicates(self):
        path = self.root / "units.tsv"
        row = self.validation["v0"]
        write_tsv(path, row.keys(), [row, row])
        with self.assertRaisesRegex(ValueError, "duplicate biological"):
            _read_units(path, validation=True)

    def test_complete_evidence_import_binds_all_sources(self):
        from adzuki_gwas_analysis.loader import compute_sha256

        source = self.root / "assay"
        source.mkdir()
        design = {**self.design, "data_scope": "synthetic"}
        write_tsv(source / "marker_states.tsv", design.keys(), [design])
        finish_provenance(
            source,
            operation="assay_review",
            inputs={},
            initial={},
            parameters={},
            environment=generation_environment(),
            families=[],
        )
        discovery, validation, external = (
            self.root / name for name in ("discovery.tsv", "validation.tsv", "external.json")
        )
        write_tsv(discovery, ("sample_id", "family_id"), self.discovery.values())
        write_tsv(validation, self.validation["v0"].keys(), self.validation.values())
        external.write_text(json.dumps({"description": "synthetic external analysis fixture"}))
        meta = copy.deepcopy(self.meta)
        meta.update(
            data_scope="synthetic",
            discovery_data_sha256=compute_sha256(discovery),
            validation_data_sha256=compute_sha256(validation),
            analysis_run_sha256=compute_sha256(external),
        )
        evidence = self.root / "evidence.json"
        evidence.write_text(json.dumps(meta))
        target = self.root / "complete"
        run_trait_review(
            assay_dir=source,
            output_dir=target,
            evidence_path=evidence,
            discovery_path=discovery,
            validation_path=validation,
            analysis_run_path=external,
        )
        report = json.loads((target / "trait_review.json").read_text())
        self.assertEqual(report["assessment"]["n_biological_samples"], 4)
        self.assertEqual(report["assessment"]["trait_utility"], "not_assessed")
        self.assertIsNotNone(validate_provenance(target, required=True))
