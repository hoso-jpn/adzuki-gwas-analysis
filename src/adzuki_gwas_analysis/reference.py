"""Offline, checksum-bound reference assets and indexed 1-based sequence access."""

from __future__ import annotations

import argparse
import json
import re
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.loader import compute_sha256


class ReferenceError(ValueError):
    """A reference asset or coordinate violates the declared contract."""


def required_text(table: dict[str, Any], key: str) -> str:
    value = table.get(key)
    if not isinstance(value, str) or not value or any(c in value for c in "\t\r\n"):
        raise ReferenceError(f"{key}: expected nonempty single-line text")
    return value


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


@dataclass(frozen=True)
class Contig:
    name: str
    length: int
    offset: int
    line_bases: int
    line_bytes: int


@dataclass(frozen=True)
class ReferenceBundle:
    reference: str
    assembly_id: str
    species: str
    source: str
    license: str
    data_scope: str
    dataset_ids: tuple[str, ...]
    fasta: Path
    contigs: tuple[Contig, ...]
    checksums: dict[str, str]
    annotation: Path | None = None
    mask: Path | None = None

    def check_dataset(self, dataset_id: str, reference: str, assembly_id: str) -> None:
        if dataset_id not in self.dataset_ids:
            raise ReferenceError("dataset_id is not bound to this reference bundle")
        if (reference, assembly_id) != (self.reference, self.assembly_id):
            raise ReferenceError("reference/assembly_id mismatch")

    def contig(self, name: str) -> Contig:
        for contig in self.contigs:
            if contig.name == name:
                return contig
        raise ReferenceError(f"unknown contig: {name}")

    def sequence(self, chrom: str, start: int, end: int, strand: str = "+") -> str:
        contig = self.contig(chrom)
        if strand not in ("+", "-") or not 1 <= start <= end <= contig.length:
            raise ReferenceError("invalid strand or 1-based inclusive interval")
        chunks: list[bytes] = []
        offset = start - 1
        with self.fasta.open("rb") as handle:
            while offset < end:
                row, column = divmod(offset, contig.line_bases)
                length = min(contig.line_bases - column, end - offset)
                handle.seek(contig.offset + row * contig.line_bytes + column)
                chunks.append(handle.read(length))
                offset += length
        sequence = b"".join(chunks).decode("ascii").upper()
        if len(sequence) != end - start + 1 or not re.fullmatch("[ACGTRYSWKMBDHVN]+", sequence):
            raise ReferenceError("FASTA changed or index no longer matches sequence")
        return reverse_complement(sequence) if strand == "-" else sequence

    def check_snp(self, chrom: str, pos: int, ref: str, alt: str) -> None:
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
            raise ReferenceError("only explicit uppercase biallelic A/C/G/T SNPs are supported")
        if ref == alt or self.sequence(chrom, pos, pos) != ref:
            raise ReferenceError("reference allele mismatch or identical REF/ALT")

    def metadata(self) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "reference": self.reference,
            "assembly_id": self.assembly_id,
            "species": self.species,
            "data_scope": self.data_scope,
            "dataset_ids": list(self.dataset_ids),
            "asset_sha256": self.checksums,
            "contigs": [{"name": c.name, "length": c.length} for c in self.contigs],
            "coordinate_system": "1-based-inclusive",
        }


def _asset(root: Path, table: dict[str, Any], assembly_id: str) -> tuple[Path, str]:
    if required_text(table, "assembly_id") != assembly_id:
        raise ReferenceError("asset assembly_id mismatch")
    relative = Path(required_text(table, "path"))
    if relative.is_absolute() or ".." in relative.parts:
        raise ReferenceError("asset path must be relative and stay inside the bundle")
    path = root / relative
    if any(p.is_symlink() for p in (path, *path.parents)) or not path.is_file():
        raise ReferenceError("asset is missing, not a file, or uses a symlink")
    digest = required_text(table, "sha256")
    if not re.fullmatch("[0-9a-f]{64}", digest) or compute_sha256(path) != digest:
        raise ReferenceError("asset checksum mismatch")
    return path, digest


def _scan_fasta(path: Path) -> tuple[Contig, ...]:
    """Independently compute FAI without retaining genome sequences in memory."""
    result: list[Contig] = []
    name = ""
    length = offset = line_bases = line_bytes = previous_bases = previous_bytes = 0
    with path.open("rb") as handle:
        while raw := handle.readline():
            if raw.startswith(b">"):
                if name:
                    if not length:
                        raise ReferenceError("empty FASTA contig")
                    result.append(Contig(name, length, offset, line_bases, line_bytes))
                header = raw[1:].strip().split()
                if not header:
                    raise ReferenceError("empty FASTA header")
                name = header[0].decode("ascii")
                if name in {c.name for c in result}:
                    raise ReferenceError("duplicate FASTA contig")
                length = line_bases = line_bytes = previous_bases = previous_bytes = 0
                offset = handle.tell()
                continue
            bases = raw.rstrip(b"\r\n")
            if not name or not re.fullmatch(b"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", bases):
                raise ReferenceError("invalid FASTA sequence")
            if not line_bases:
                line_bases, line_bytes = len(bases), len(raw)
            elif previous_bases != line_bases or previous_bytes != line_bytes:
                raise ReferenceError("FASTA wrapping changed before the last sequence line")
            if len(bases) > line_bases:
                raise ReferenceError("FASTA line exceeds indexed width")
            previous_bases, previous_bytes = len(bases), len(raw)
            length += len(bases)
    if not name or not length:
        raise ReferenceError("empty FASTA or contig")
    result.append(Contig(name, length, offset, line_bases, line_bytes))
    return tuple(result)


def _read_fai(path: Path) -> tuple[Contig, ...]:
    entries: list[Contig] = []
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split("\t")
        if len(fields) != 5:
            raise ReferenceError("FAI must contain five tab-separated fields")
        try:
            entries.append(Contig(fields[0], *(int(value) for value in fields[1:])))
        except ValueError as exc:
            raise ReferenceError("invalid FAI integer") from exc
    return tuple(entries)


def _check_intervals(path: Path, contigs: tuple[Contig, ...], *, bed: bool) -> None:
    lengths = {c.name: c.length for c in contigs}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                raise ReferenceError("embedded FASTA in GFF is not supported")
            if line.startswith("#") or not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            if (bed and len(columns) < 3) or (not bed and len(columns) != 9):
                raise ReferenceError("invalid BED/GFF/GTF record")
            chrom = columns[0]
            try:
                start, end = (
                    (int(columns[1]) + 1, int(columns[2]))
                    if bed
                    else (int(columns[3]), int(columns[4]))
                )
            except ValueError as exc:
                raise ReferenceError("invalid annotation coordinate") from exc
            if chrom not in lengths or not 1 <= start <= end <= lengths[chrom]:
                raise ReferenceError("annotation/mask contig or interval mismatch")
            if not bed and columns[6] not in ("+", "-", ".", "?"):
                raise ReferenceError("invalid annotation strand")


def load_reference_bundle(path: Path) -> ReferenceBundle:
    try:
        table = tomllib.loads(path.read_text(encoding="utf-8"))
        if type(table.get("schema_version")) is not int or table["schema_version"] != 1:
            raise ReferenceError("unsupported reference bundle schema_version")
        assembly = required_text(table, "assembly_id")
        datasets = table.get("dataset_ids")
        if (
            not isinstance(datasets, list)
            or not datasets
            or any(
                not isinstance(d, str) or not d or any(c in d for c in "\t\r\n") for d in datasets
            )
            or len(set(datasets)) != len(datasets)
        ):
            raise ReferenceError("dataset_ids must be a nonempty list of unique identifiers")
        scope = required_text(table, "data_scope")
        if scope not in ("public", "customer", "synthetic"):
            raise ReferenceError("unknown data_scope")
        assets: dict[str, Path] = {}
        checksums = {"bundle": compute_sha256(path)}
        for role in ("fasta", "fai", "annotation", "mask"):
            if role not in table and role in ("annotation", "mask"):
                continue
            if not isinstance(table.get(role), dict):
                raise ReferenceError(f"missing/invalid {role} asset")
            assets[role], checksums[role] = _asset(path.parent, table[role], assembly)
        contigs = _scan_fasta(assets["fasta"])
        if _read_fai(assets["fai"]) != contigs:
            raise ReferenceError("FAI names, order, length or byte offsets mismatch FASTA")
        for role in ("annotation", "mask"):
            if role in assets:
                _check_intervals(assets[role], contigs, bed=role == "mask")
        return ReferenceBundle(
            reference=required_text(table, "reference"),
            assembly_id=assembly,
            species=required_text(table, "species"),
            source=required_text(table, "source"),
            license=required_text(table, "license"),
            data_scope=scope,
            dataset_ids=tuple(datasets),
            fasta=assets["fasta"],
            contigs=contigs,
            checksums=checksums,
            annotation=assets.get("annotation"),
            mask=assets.get("mask"),
        )
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as exc:
        raise ReferenceError("cannot read a valid reference bundle or asset") from exc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--assembly-id", required=True)
    args = parser.parse_args(argv)
    try:
        bundle = load_reference_bundle(args.bundle)
        bundle.check_dataset(args.dataset_id, args.reference, args.assembly_id)
        print(json.dumps(bundle.metadata(), sort_keys=True, indent=2))
    except ReferenceError as exc:
        parser.exit(1, f"reference validation failed: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
