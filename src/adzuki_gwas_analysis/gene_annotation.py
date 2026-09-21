"""Reference-local gene overlap/nearest distance from explicit GFF3/GTF gene features."""

from __future__ import annotations

import re
from typing import Any
from urllib.parse import unquote

from adzuki_gwas_analysis.reference import ReferenceBundle


def load_genes(bundle: ReferenceBundle) -> dict[str, list[dict[str, Any]]]:
    genes: dict[str, list[dict[str, Any]]] = {}
    if bundle.annotation is None:
        return genes
    with bundle.annotation.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if fields[2] != "gene":
                continue
            if "=" in fields[8]:
                attrs = dict(part.split("=", 1) for part in fields[8].split(";") if "=" in part)
                identifier = unquote(attrs.get("ID", ""))
            else:
                attrs = dict(re.findall(r'(\w+)\s+"([^"\r\n]*)"\s*;?', fields[8]))
                identifier = attrs.get("gene_id", "")
            if not identifier or any(c in identifier for c in "\t\r\n"):
                raise ValueError("gene feature requires a safe GFF3 ID or GTF gene_id")
            genes.setdefault(fields[0], []).append(
                {
                    "gene_id": identifier,
                    "gene_start": int(fields[3]),
                    "gene_end": int(fields[4]),
                    "gene_strand": fields[6],
                }
            )
    return genes


def annotate_candidate(
    row: dict[str, str], bundle: ReferenceBundle, genes: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    identity = {
        key: row[key]
        for key in (
            "candidate_id",
            "dataset_id",
            "reference",
            "assembly_id",
            "chr",
            "pos",
            "ref",
            "alt",
            "trait",
        )
    }
    available = genes.get(row["chr"], [])
    if not available:
        return [
            {
                **identity,
                "annotation_status": "not_assessed"
                if bundle.annotation is None
                else "no_gene_features_on_contig",
                "gene_id": "",
                "gene_start": "",
                "gene_end": "",
                "gene_strand": "",
                "distance_bp": "",
            }
        ]
    position = int(row["pos"])
    distances = [
        (max(gene["gene_start"] - position, position - gene["gene_end"], 0), gene)
        for gene in available
    ]
    nearest = min(distance for distance, _ in distances)
    return [
        {
            **identity,
            **gene,
            "distance_bp": distance,
            "annotation_status": "overlap" if distance == 0 else "nearest",
        }
        for distance, gene in distances
        if distance == nearest
    ]
