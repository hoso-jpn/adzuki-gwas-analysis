"""Reference-checked candidate flanks and source-separated neighboring variants."""

from __future__ import annotations

import argparse
import gzip
import hashlib
from bisect import bisect_left, bisect_right
from collections.abc import Iterator
from pathlib import Path
from typing import Any
from urllib.parse import quote

from adzuki_gwas_analysis.provenance import (
    finish_provenance,
    generation_environment,
    input_checksums,
    output_transaction,
)
from adzuki_gwas_analysis.reference import ReferenceBundle, load_reference_bundle
from adzuki_gwas_analysis.tables import read_tsv, write_tsv

CONTEXT_FIELDS = (
    "schema_version",
    "candidate_id",
    "dataset_id",
    "trait",
    "reference",
    "assembly_id",
    "chr",
    "pos",
    "ref",
    "alt",
    "strand",
    "interval_start",
    "interval_end",
    "snp_offset_0based",
    "sequence_sha256",
    "flanks_complete",
    "summary_neighbor_count",
    "vcf_neighbor_count",
    "mask_overlap",
    "annotation_status",
    "allele_source",
    "candidate_input_sha256",
)
NEIGHBOR_FIELDS = (
    "candidate_id",
    "dataset_id",
    "reference",
    "assembly_id",
    "chr",
    "pos",
    "ref",
    "alt",
    "source",
    "source_sha256",
    "filter",
    "distance_bp",
)


def _variant(row: dict[str, str], bundle: ReferenceBundle, dataset_id: str) -> dict[str, Any]:
    bundle.check_dataset(
        row.get("dataset_id", dataset_id),
        row.get("reference", bundle.reference),
        row.get("assembly_id", bundle.assembly_id),
    )
    if row.get("dataset_id", dataset_id) != dataset_id:
        raise ValueError("mixed datasets are not allowed in one context run")
    chrom, pos = row["chr"], int(row["pos"])
    if "ref" in row and "alt" in row:
        ref, alt, source = row["ref"], row["alt"], "explicit_REF_ALT"
    elif "allele1" in row and "allele0" in row:
        base = bundle.sequence(chrom, pos, pos)
        first, other = row["allele1"], row["allele0"]
        if first == other or base not in (first, other):
            raise ValueError("neither unique GWAS allele matches FASTA")
        ref, alt = (first, other) if base == first else (other, first)
        source = "FASTA_match_of_allele1_allele0;effect_orientation_not_inferred"
    else:
        raise ValueError("explicit REF/ALT or both GWAS alleles are required")
    bundle.check_snp(chrom, pos, ref, alt)
    return {"chr": chrom, "pos": pos, "ref": ref, "alt": alt, "allele_source": source}


def _vcf_rows(path: Path, bundle: ReferenceBundle) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    reference_seen = header_seen = False
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##reference="):
                if reference_seen or line.strip().partition("=")[2] != bundle.assembly_id:
                    raise ValueError("VCF reference does not uniquely match assembly_id")
                reference_seen = True
            elif line.startswith("#CHROM\t"):
                if not reference_seen or header_seen:
                    raise ValueError("VCF needs one matching ##reference and one column header")
                header_seen = True
            elif not line.startswith("#"):
                fields = line.rstrip("\r\n").split("\t")
                if not header_seen or len(fields) < 8:
                    raise ValueError("invalid VCF header/record")
                for alt in fields[4].split(","):
                    yield {
                        "chr": fields[0],
                        "pos": fields[1],
                        "ref": fields[3],
                        "alt": alt,
                        "filter": fields[6],
                    }
    if not header_seen:
        raise ValueError("VCF column header is missing")


def _masked(bundle: ReferenceBundle, chrom: str, start: int, end: int) -> str:
    if bundle.mask is None:
        return "unavailable"
    with bundle.mask.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if fields[0] == chrom and int(fields[1]) < end and int(fields[2]) >= start:
                return "true"
    return "false"


def run_context(
    *,
    candidate_table: Path,
    bundle_path: Path,
    dataset_id: str,
    output_dir: Path,
    flank_bp: int,
    neighbor_window_bp: int,
    summary_table: Path | None = None,
    cohort_vcf: Path | None = None,
) -> None:
    if flank_bp < 0 or neighbor_window_bp < 0:
        raise ValueError("flank and neighbor windows must be nonnegative")
    bundle = load_reference_bundle(bundle_path)
    bundle.check_dataset(dataset_id, bundle.reference, bundle.assembly_id)
    inputs = {
        "candidates": candidate_table,
        "reference_bundle": bundle_path,
        "reference_fasta": bundle.fasta,
    }
    for role, path in (
        ("summary", summary_table),
        ("cohort_vcf", cohort_vcf),
        ("annotation", bundle.annotation),
        ("mask", bundle.mask),
    ):
        if path is not None:
            inputs[role] = path
    initial, environment = input_checksums(inputs), generation_environment()
    candidates: list[dict[str, Any]] = []
    identifiers: set[str] = set()
    loci: set[tuple[str, int, str, str]] = set()
    for row in read_tsv(candidate_table, ("chr", "pos", "trait")):
        variant = _variant(row, bundle, dataset_id)
        key = (variant["chr"], variant["pos"], variant["ref"], variant["alt"])
        identity = hashlib.sha256(
            (dataset_id + "\0" + bundle.assembly_id + "\0" + repr(key)).encode()
        ).hexdigest()[:20]
        candidate_id = row.get("candidate_id", "candidate-" + identity)
        if not candidate_id or candidate_id in identifiers or key in loci or not row["trait"]:
            raise ValueError("duplicate/empty candidate identity or missing trait")
        identifiers.add(candidate_id)
        loci.add(key)
        chrom, pos = variant["chr"], variant["pos"]
        start, end = max(1, pos - flank_bp), min(bundle.contig(chrom).length, pos + flank_bp)
        strand = row.get("strand", "+")
        sequence = bundle.sequence(chrom, start, end, strand)
        candidates.append(
            {
                "schema_version": 1,
                "candidate_id": candidate_id,
                "dataset_id": dataset_id,
                "trait": row["trait"],
                "reference": bundle.reference,
                "assembly_id": bundle.assembly_id,
                **variant,
                "strand": strand,
                "interval_start": start,
                "interval_end": end,
                "snp_offset_0based": pos - start if strand == "+" else end - pos,
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "flanks_complete": start == pos - flank_bp and end == pos + flank_bp,
                "summary_neighbor_count": 0,
                "vcf_neighbor_count": 0,
                "mask_overlap": _masked(bundle, chrom, start, end),
                "annotation_status": "available" if bundle.annotation else "unavailable",
                "candidate_input_sha256": initial["candidates"],
            }
        )
    by_chrom = {
        chrom.name: sorted(
            [c for c in candidates if c["chr"] == chrom.name],
            key=lambda candidate: candidate["pos"],
        )
        for chrom in bundle.contigs
    }
    positions = {chrom: [c["pos"] for c in values] for chrom, values in by_chrom.items()}
    neighbors: list[dict[str, Any]] = []
    seen: set[tuple[str, str, int, str, str]] = set()
    for source, path in (("summary", summary_table), ("cohort_vcf", cohort_vcf)):
        if path is None:
            continue
        rows = _vcf_rows(path, bundle) if source == "cohort_vcf" else read_tsv(path, ("chr", "pos"))
        for row in rows:
            chrom, pos = row["chr"], int(row["pos"])
            if not 1 <= pos <= bundle.contig(chrom).length:
                raise ValueError("neighbor coordinate out of range")
            end = pos + len(row.get("ref", "N")) - 1
            nearby = by_chrom[chrom][
                bisect_left(positions[chrom], pos - neighbor_window_bp) : bisect_right(
                    positions[chrom], end + neighbor_window_bp
                )
            ]
            if not nearby:
                continue
            if source == "summary":
                variant = _variant(row, bundle, dataset_id)
            else:
                ref, alt = row["ref"], row["alt"]
                if (
                    not ref
                    or not alt
                    or any(base not in "ACGT" for base in ref + alt)
                    or ref == alt
                    or bundle.sequence(chrom, pos, end) != ref
                ):
                    raise ValueError("VCF neighbor REF/ALT is unsupported or inconsistent")
                variant = {"chr": chrom, "pos": pos, "ref": ref, "alt": alt}
            for candidate in nearby:
                if all(variant[k] == candidate[k] for k in ("chr", "pos", "ref", "alt")):
                    continue
                key_neighbor = (
                    candidate["candidate_id"],
                    source,
                    pos,
                    variant["ref"],
                    variant["alt"],
                )
                if key_neighbor in seen:
                    raise ValueError("duplicate neighboring variant record")
                seen.add(key_neighbor)
                candidate[
                    "vcf_neighbor_count" if source == "cohort_vcf" else "summary_neighbor_count"
                ] += 1
                neighbors.append(
                    {
                        "candidate_id": candidate["candidate_id"],
                        "dataset_id": dataset_id,
                        "reference": bundle.reference,
                        "assembly_id": bundle.assembly_id,
                        **variant,
                        "source": source,
                        "source_sha256": initial[source],
                        "filter": row.get("filter", "not_provided"),
                        "distance_bp": pos - candidate["pos"],
                    }
                )
    with output_transaction(output_dir) as stage:
        write_tsv(stage / "candidate_context.tsv", CONTEXT_FIELDS, candidates)
        write_tsv(stage / "neighboring_variants.tsv", NEIGHBOR_FIELDS, neighbors)
        with (stage / "flanking_sequences.fasta").open("w") as handle:
            for candidate in candidates:
                header = " ".join(
                    f"{key}={quote(str(candidate[key]), safe='')}"
                    for key in ("dataset_id", "reference", "assembly_id", "chr", "pos", "strand")
                )
                handle.write(
                    f">{quote(candidate['candidate_id'], safe='')} {header} "
                    f"sha256={candidate['sequence_sha256']}\n{candidate['sequence']}\n"
                )
        finish_provenance(
            stage,
            operation="candidate_context",
            inputs=inputs,
            initial=initial,
            parameters={
                "flank_bp": flank_bp,
                "neighbor_window_bp": neighbor_window_bp,
                "reference": bundle.metadata(),
            },
            environment=environment,
            families=[{"dataset_id": dataset_id, "scope": "context only; no LD inference"}],
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("candidate-table", "bundle-path", "output-dir"):
        parser.add_argument("--" + option, type=Path, required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--flank-bp", type=int, required=True)
    parser.add_argument("--neighbor-window-bp", type=int, required=True)
    parser.add_argument("--summary-table", type=Path)
    parser.add_argument("--cohort-vcf", type=Path)
    args = parser.parse_args(argv)
    try:
        run_context(**vars(args))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"context generation failed: {type(exc).__name__}: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
