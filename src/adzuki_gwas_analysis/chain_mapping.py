"""Stream UCSC chain blocks: source is q/query, destination is t/target."""

from __future__ import annotations

import bisect
import itertools
import math
import re
import tomllib
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.loader import compute_sha256
from adzuki_gwas_analysis.reference import ReferenceBundle, required_text


def load_chain_metadata(
    path: Path, source: ReferenceBundle, target: ReferenceBundle
) -> tuple[dict[str, Any], Path]:
    meta = tomllib.loads(path.read_text(encoding="utf-8"))
    if type(meta.get("schema_version")) is not int or meta["schema_version"] != 1:
        raise ValueError("unsupported chain metadata schema")
    fields = (
        "source_assembly",
        "target_assembly",
        "source_uri",
        "version",
        "license",
        "validation_reference",
        "reviewer_id",
        "chain_path",
        "chain_sha256",
    )
    for key in fields:
        required_text(meta, key)
    if set(meta) != set(fields) | {"schema_version", "reviewed", "direction"}:
        raise ValueError("missing/unknown chain metadata field")
    if meta["reviewed"] is not True or meta["direction"] != "query_to_target":
        raise ValueError("chain must have reviewed evidence and explicit query_to_target direction")
    if (meta["source_assembly"], meta["target_assembly"]) != (
        source.assembly_id,
        target.assembly_id,
    ):
        raise ValueError("chain assembly mismatch")
    relative = Path(meta["chain_path"])
    chain = path.parent / relative
    if (
        relative.is_absolute()
        or ".." in relative.parts
        or any(p.is_symlink() for p in (chain, *chain.parents))
    ):
        raise ValueError("chain path must stay inside the metadata directory without symlinks")
    if (
        not re.fullmatch("[0-9a-f]{64}", meta["chain_sha256"])
        or compute_sha256(chain) != meta["chain_sha256"]
    ):
        raise ValueError("chain checksum mismatch")
    return meta, chain


def map_chain(
    chain: Path,
    source: ReferenceBundle,
    target: ReferenceBundle,
    candidates: list[dict[str, str]],
    *,
    max_blocks: int = 1_000_000,
    max_total_hits: int = 100_000,
) -> dict[str, list[dict[str, Any]]]:
    """Map SNP points, preserving every hit; gaps are never linearly interpolated."""
    indexed: dict[str, list[tuple[int, str]]] = {}
    hits: dict[str, list[dict[str, Any]]] = {row["candidate_id"]: [] for row in candidates}
    for row in candidates:
        indexed.setdefault(row["chr"], []).append((int(row["pos"]), row["candidate_id"]))
    for points in indexed.values():
        points.sort()
    chain_ids: set[str] = set()
    active: dict[str, Any] | None = None
    blocks = total_hits = 0
    with chain.open(encoding="ascii") as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if fields[0] == "chain":
                if active is not None or len(fields) != 13:
                    raise ValueError("incomplete or invalid chain header")
                (
                    _,
                    score_text,
                    tname,
                    tsize_text,
                    tstrand,
                    ts,
                    te,
                    qname,
                    qsize_text,
                    qstrand,
                    qs,
                    qe,
                    identifier,
                ) = fields
                tsize, qsize = int(tsize_text), int(qsize_text)
                tstart, tend, qstart, qend = map(int, (ts, te, qs, qe))
                score = float(score_text)
                if identifier in chain_ids or not math.isfinite(score) or score < 0:
                    raise ValueError("duplicate chain ID or invalid score")
                chain_ids.add(identifier)
                if tstrand != "+" or qstrand not in ("+", "-"):
                    raise ValueError("v1 requires target strand + and query strand +/-")
                if target.contig(tname).length != tsize or source.contig(qname).length != qsize:
                    raise ValueError("chain contig/assembly length mismatch")
                if not (0 <= tstart < tend <= tsize and 0 <= qstart < qend <= qsize):
                    raise ValueError("chain header interval is outside reference")
                active = dict(
                    tname=tname,
                    qname=qname,
                    qsize=qsize,
                    strand=qstrand,
                    tcursor=tstart,
                    qcursor=qstart,
                    tend=tend,
                    qend=qend,
                    chain_id=identifier,
                    score=score,
                )
                continue
            if active is None or len(fields) not in (1, 3):
                raise ValueError("unexpected chain block")
            size, *gaps = map(int, fields)
            blocks += 1
            if size <= 0 or any(gap < 0 for gap in gaps) or blocks > max_blocks:
                raise ValueError("invalid chain block or declared block limit exceeded")
            a = active
            if a["tcursor"] + size > a["tend"] or a["qcursor"] + size > a["qend"]:
                raise ValueError("chain blocks exceed header endpoints")
            start = a["qcursor"] if a["strand"] == "+" else a["qsize"] - a["qcursor"] - size
            points = indexed.get(a["qname"], [])
            left = bisect.bisect_left(points, (start + 1, ""))
            for position, candidate_id in itertools.islice(points, left, None):
                if position > start + size:
                    break
                query_offset = position - 1 if a["strand"] == "+" else a["qsize"] - position
                target_position = a["tcursor"] + query_offset - a["qcursor"] + 1
                total_hits += 1
                if total_hits > max_total_hits:
                    raise ValueError("declared mapping hit limit exceeded")
                hits[candidate_id].append(
                    {
                        "target_chr": a["tname"],
                        "target_pos": target_position,
                        "strand": a["strand"],
                        "chain_id": a["chain_id"],
                        "chain_score": a["score"],
                        "block_size": size,
                    }
                )
            a["tcursor"] += size
            a["qcursor"] += size
            if gaps:
                a["tcursor"] += gaps[0]
                a["qcursor"] += gaps[1]
            else:
                if (a["tcursor"], a["qcursor"]) != (a["tend"], a["qend"]):
                    raise ValueError("chain terminal block does not reach declared endpoints")
                active = None
    if active is not None or not chain_ids:
        raise ValueError("chain is empty or missing its terminal block")
    return hits
