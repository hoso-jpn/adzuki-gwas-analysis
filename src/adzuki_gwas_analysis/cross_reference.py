"""Compare candidates only through a reviewed, checksum-bound assembly alignment."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.chain_mapping import load_chain_metadata, map_chain
from adzuki_gwas_analysis.gene_annotation import annotate_candidate, load_genes
from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    output_transaction,
    write_json,
)
from adzuki_gwas_analysis.reference import (
    ReferenceBundle,
    load_reference_bundle,
    reverse_complement,
)
from adzuki_gwas_analysis.tables import read_tsv, write_tsv

CANDIDATE_FIELDS = (
    "candidate_id",
    "dataset_id",
    "cohort_id",
    "analysis_id",
    "reference",
    "assembly_id",
    "chr",
    "pos",
    "ref",
    "alt",
    "trait",
    "trait_unit",
    "trait_coding",
    "test",
    "effect_scale",
    "effect_allele",
    "other_allele",
    "effect_orientation",
    "beta",
    "se",
    "neg_log10_pvalue",
)
ANNOTATION_FIELDS = (
    "side",
    "candidate_id",
    "dataset_id",
    "reference",
    "assembly_id",
    "chr",
    "pos",
    "ref",
    "alt",
    "trait",
    "annotation_status",
    "gene_id",
    "gene_start",
    "gene_end",
    "gene_strand",
    "distance_bp",
)
COMPARISON_FIELDS = tuple(
    f"{side}_{field}" for side in ("source", "target") for field in CANDIDATE_FIELDS
) + (
    "source_gene_ids",
    "target_gene_ids",
    "strand",
    "chain_id",
    "chain_score",
    "block_size",
    "allele_alignment",
    "effect_direction_comparison",
    "cohort_relationship",
    "interpretation",
)
MAPPING_FIELDS = (
    "candidate_id",
    "source_assembly",
    "source_chr",
    "source_pos",
    "source_ref",
    "source_alt",
    "status",
    "reason",
    "target_assembly",
    "target_chr",
    "target_pos",
    "target_reference_base",
    "mapped_ref",
    "mapped_alt",
    "strand",
    "chain_id",
    "chain_score",
    "block_size",
)


def load_candidates(
    path: Path, bundle: ReferenceBundle, trait: str, max_candidates: int
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    identifiers: set[str] = set()
    variants: set[tuple[str, int, str, str]] = set()
    families: set[tuple[str, ...]] = set()
    for row in read_tsv(path, CANDIDATE_FIELDS):
        if any(not row[key] for key in CANDIDATE_FIELDS):
            raise ValueError("candidate identity/statistical fields must be explicit")
        if len(rows) >= max_candidates:
            raise ValueError("candidate count exceeds declared operating limit")
        bundle.check_dataset(row["dataset_id"], row["reference"], row["assembly_id"])
        bundle.check_snp(row["chr"], int(row["pos"]), row["ref"], row["alt"])
        if row["trait"] != trait:
            raise ValueError("only the explicitly requested trait can be compared")
        if (
            any(not math.isfinite(float(row[key])) for key in ("beta", "se", "neg_log10_pvalue"))
            or float(row["se"]) <= 0
            or float(row["neg_log10_pvalue"]) < 0
        ):
            raise ValueError("invalid candidate effect/SE/p-value")
        if row["effect_orientation"] not in ("as_reported", "aligned_to_ALT", "unresolved_strand"):
            raise ValueError("unknown effect orientation")
        effect_pair = {row["effect_allele"], row["other_allele"]}
        allowed_pairs = [{row["ref"], row["alt"]}]
        if row["effect_orientation"] == "unresolved_strand":
            allowed_pairs.append({reverse_complement(row["ref"]), reverse_complement(row["alt"])})
        if effect_pair not in allowed_pairs:
            raise ValueError("effect alleles disagree with declared orientation")
        variant = (row["chr"], int(row["pos"]), row["ref"], row["alt"])
        if row["candidate_id"] in identifiers or variant in variants:
            raise ValueError("duplicate candidate or variant identity")
        identifiers.add(row["candidate_id"])
        variants.add(variant)
        families.add(
            tuple(
                row[key]
                for key in (
                    "dataset_id",
                    "cohort_id",
                    "analysis_id",
                    "trait",
                    "trait_unit",
                    "trait_coding",
                    "test",
                    "effect_scale",
                )
            )
        )
        rows.append(row)
    if len(families) > 1:
        raise ValueError("each candidate table must have exactly one analysis family")
    return rows


def _effect_comparison(source: dict[str, str], target: dict[str, str], strand: str) -> str:
    if any(
        row[key] == "unknown"
        for row in (source, target)
        for key in ("test", "effect_scale", "trait_coding")
    ):
        return "not_assessed_unknown_effect_or_test_contract"
    if any(row["effect_orientation"] == "unresolved_strand" for row in (source, target)):
        return "not_assessed_unresolved_strand"
    if any(
        source[key] != target[key] for key in ("trait_unit", "trait_coding", "effect_scale", "test")
    ):
        return "not_assessed_different_effect_or_test_contract"
    allele = (
        reverse_complement(source["effect_allele"]) if strand == "-" else source["effect_allele"]
    )
    effect = float(source["beta"])
    if allele == target["other_allele"]:
        effect *= -1
    elif allele != target["effect_allele"]:
        return "not_assessed_allele_unresolved"
    product = effect * float(target["beta"])
    return (
        "zero_effect"
        if product == 0
        else "same_sign_after_alignment"
        if product > 0
        else "opposite_sign_after_alignment"
    )


def run_cross_reference(
    *,
    source_table: Path,
    source_bundle: Path,
    trait: str,
    output_dir: Path,
    target_table: Path | None = None,
    target_bundle: Path | None = None,
    chain_metadata: Path | None = None,
    max_candidates: int = 10_000,
) -> None:
    if not trait or max_candidates <= 0 or (target_table is None) != (target_bundle is None):
        raise ValueError("trait, valid limit, and paired target table/bundle are required")
    if chain_metadata is not None and target_bundle is None:
        raise ValueError("chain comparison requires a target table and bundle")
    source = load_reference_bundle(source_bundle)
    target = load_reference_bundle(target_bundle) if target_bundle is not None else None
    inputs = {
        "source_candidates": source_table,
        "source_bundle": source_bundle,
        "source_fasta": source.fasta,
    }
    if source.annotation is not None:
        inputs["source_annotation"] = source.annotation
    if target is not None and target_table is not None and target_bundle is not None:
        if target.assembly_id == source.assembly_id or target.species != source.species:
            raise ValueError(
                "cross-reference comparison requires different assemblies of the same species"
            )
        inputs.update(
            target_candidates=target_table, target_bundle=target_bundle, target_fasta=target.fasta
        )
        if target.annotation is not None:
            inputs["target_annotation"] = target.annotation
    meta: dict[str, Any] | None = None
    chain: Path | None = None
    if chain_metadata is not None and target is not None:
        meta, chain = load_chain_metadata(chain_metadata, source, target)
        inputs.update(chain_metadata=chain_metadata, chain=chain)
    initial, environment = input_checksums(inputs), generation_environment()
    source_rows = load_candidates(source_table, source, trait, max_candidates)
    target_rows = (
        load_candidates(target_table, target, trait, max_candidates)
        if target is not None and target_table is not None
        else []
    )
    annotations: list[dict[str, Any]] = []
    gene_ids: dict[tuple[str, str], str] = {}
    for side, bundle, rows in (("source", source, source_rows), ("target", target, target_rows)):
        if bundle is None:
            continue
        genes = load_genes(bundle)
        for row in rows:
            entries = annotate_candidate(row, bundle, genes)
            annotations.extend({"side": side, **entry} for entry in entries)
            gene_ids[(side, row["candidate_id"])] = ";".join(
                sorted({entry["gene_id"] for entry in entries if entry["gene_id"]})
            )
    hits = (
        map_chain(chain, source, target, source_rows)
        if chain is not None and target is not None
        else {}
    )
    target_at: dict[tuple[str, int], list[dict[str, str]]] = {}
    for row in target_rows:
        target_at.setdefault((row["chr"], int(row["pos"])), []).append(row)
    mappings: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    accounting: list[dict[str, Any]] = []
    compared_targets: dict[str, list[str]] = {}
    for row in source_rows:
        available = hits.get(row["candidate_id"], [])
        base = {
            "candidate_id": row["candidate_id"],
            "source_assembly": source.assembly_id,
            "source_chr": row["chr"],
            "source_pos": row["pos"],
            "source_ref": row["ref"],
            "source_alt": row["alt"],
            "target_assembly": target.assembly_id if target is not None else "",
        }
        status, reason = "not_assessed", "no_correspondence_asset"
        if chain is not None:
            status, reason = "unmapped", "no_ungapped_block_for_source_position"
        if available and target is not None:
            status = "multimap" if len(available) > 1 else "mapped"
            reason = "multiple_chain_hits" if len(available) > 1 else "unique_in_supplied_chain"
            for hit in available:
                ref, alt = row["ref"], row["alt"]
                if hit["strand"] == "-":
                    ref, alt = reverse_complement(ref), reverse_complement(alt)
                target_base = target.sequence(
                    hit["target_chr"], hit["target_pos"], hit["target_pos"]
                )
                local_status, local_reason = status, reason
                if status == "mapped" and target_base not in (ref, alt):
                    local_status, local_reason = (
                        "ambiguous",
                        "target_reference_allele_not_reconciled",
                    )
                    status, reason = local_status, local_reason
                mappings.append(
                    {
                        **base,
                        **hit,
                        "target_reference_base": target_base,
                        "mapped_ref": ref,
                        "mapped_alt": alt,
                        "status": local_status,
                        "reason": local_reason,
                    }
                )
                if local_status != "mapped":
                    continue
                matching = [
                    other
                    for other in target_at.get((hit["target_chr"], hit["target_pos"]), [])
                    if {other["ref"], other["alt"]} == {ref, alt}
                ]
                if not matching:
                    reason = "mapped_coordinate_no_matching_target_candidate"
                for other in matching:
                    comparisons.append(
                        {
                            **{f"source_{k}": row[k] for k in CANDIDATE_FIELDS},
                            **{f"target_{k}": other[k] for k in CANDIDATE_FIELDS},
                            **hit,
                            "source_gene_ids": gene_ids[("source", row["candidate_id"])],
                            "target_gene_ids": gene_ids[("target", other["candidate_id"])],
                            "allele_alignment": "same_REF_ALT"
                            if target_base == ref
                            else "REF_ALT_swapped",
                            "effect_direction_comparison": _effect_comparison(
                                row, other, hit["strand"]
                            ),
                            "cohort_relationship": "same_cohort"
                            if row["cohort_id"] == other["cohort_id"]
                            else "different_ids_independence_unverified",
                            "interpretation": "candidate_context_not_replication_or_meta_analysis",
                        }
                    )
                    compared_targets.setdefault(other["candidate_id"], []).append(
                        row["candidate_id"]
                    )
        else:
            mappings.append({**base, "status": status, "reason": reason})
        accounting.append(
            {
                "side": "source",
                "candidate_id": row["candidate_id"],
                "status": status,
                "reason": reason,
                "n_chain_hits": len(available),
                "linked_candidates": "",
            }
        )
    for row in target_rows:
        linked = compared_targets.get(row["candidate_id"], [])
        accounting.append(
            {
                "side": "target",
                "candidate_id": row["candidate_id"],
                "status": "compared" if linked else "not_compared",
                "reason": "linked_by_verified_alleles"
                if linked
                else "no_source_comparison_not_evidence_of_absence",
                "n_chain_hits": "",
                "linked_candidates": ";".join(linked),
            }
        )
    with output_transaction(output_dir) as stage:
        write_tsv(stage / "candidate_annotations.tsv", ANNOTATION_FIELDS, annotations)
        write_tsv(stage / "coordinate_mappings.tsv", MAPPING_FIELDS, mappings)
        write_tsv(stage / "cross_reference_comparison.tsv", COMPARISON_FIELDS, comparisons)
        write_tsv(
            stage / "candidate_accounting.tsv",
            ("side", "candidate_id", "status", "reason", "n_chain_hits", "linked_candidates"),
            accounting,
        )
        summary = {
            "schema_version": 1,
            "trait": trait,
            "n_source": len(source_rows),
            "n_target": len(target_rows),
            "n_comparisons": len(comparisons),
            "cross_reference_status": "assessed_against_supplied_chain"
            if chain is not None
            else "not_assessed",
            "source_status_counts": {
                state: sum(
                    entry["side"] == "source" and entry["status"] == state for entry in accounting
                )
                for state in ("mapped", "unmapped", "multimap", "ambiguous", "not_assessed")
            },
        }
        write_json(stage / "comparison_summary.json", summary)
        (stage / "comparison_report.md").write_text(
            "# Cross-reference candidate context\n\n"
            f"Trait: {trait}. Source candidates: {len(source_rows)}; "
            f"target candidates: {len(target_rows)}. "
            f"Matched comparisons: {len(comparisons)}. "
            f"Status: {summary['cross_reference_status']}.\n\n"
            "See candidate_accounting.tsv for every candidate, "
            "coordinate_mappings.tsv for every hit, "
            "and candidate_annotations.tsv for reference-local gene overlaps/nearest distances. "
            "When supplied, the chain has external review evidence and a checksum "
            "recorded in analysis_run.json. "
            "A unique hit means unique within the supplied alignment, "
            "not proven genome-wide orthology.\n\n"
            "Unmapped, ambiguous, or unselected candidates are not evidence of absent association. "
            "No p-values are recalculated or pooled. Reference-specific families remain separate. "
            "The same cohort mapped to two references is not independent replication; "
            "different cohort IDs do not establish independence either. "
            "Raw p-value differences are not biological contrasts. "
            "Gene proximity does not identify causal genes, LD blocks or independent QTLs.\n",
            encoding="utf-8",
        )
        finish_provenance(
            stage,
            operation="cross_reference",
            inputs=inputs,
            initial=initial,
            environment=environment,
            parameters={
                "trait": trait,
                "max_candidates": max_candidates,
                "chain_metadata": meta,
                "source": source.metadata(),
                "target": target.metadata() if target is not None else None,
            },
            families=[
                {
                    key: rows[0][key]
                    for key in ("dataset_id", "cohort_id", "analysis_id", "trait", "test")
                }
                for rows in (source_rows, target_rows)
                if rows
            ],
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("source-table", "source-bundle", "output-dir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    for flag in ("target-table", "target-bundle", "chain-metadata"):
        parser.add_argument("--" + flag, type=Path)
    parser.add_argument("--trait", required=True)
    parser.add_argument("--max-candidates", type=int, default=10_000)
    try:
        run_cross_reference(**vars(parser.parse_args(argv)))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"cross-reference comparison failed: {type(exc).__name__}: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
