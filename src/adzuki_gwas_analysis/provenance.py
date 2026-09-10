"""Deterministic analysis-generation provenance and atomic output publication."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from importlib import metadata
from pathlib import Path
from typing import Any

from adzuki_gwas_analysis.analysis.output_safety import check_output_dir_is_safe
from adzuki_gwas_analysis.loader import compute_sha256

PROVENANCE_FILE = "analysis_run.json"
CONTRACT_FILE = "bundle_contract.json"


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def write_json(path: Path, value: object) -> None:
    path.write_bytes(json_bytes(value))


@contextmanager
def output_transaction(output: Path) -> Iterator[Path]:
    check_output_dir_is_safe(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".analysis-stage-", dir=output.parent))
    try:
        yield stage
        check_output_dir_is_safe(output)
        os.replace(stage, output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def generation_environment() -> dict[str, Any]:
    package = Path(__file__).resolve().parent
    source = hashlib.sha256()
    for path in sorted(package.rglob("*.py")):
        source.update(path.relative_to(package).as_posix().encode() + b"\0")
        source.update(bytes.fromhex(compute_sha256(path)))
    root = package.parent.parent
    lock_hash: str | None = None
    commit: str | None = None
    if (root / "uv.lock").is_file():
        lock_hash = compute_sha256(root / "uv.lock")
        try:
            value = subprocess.check_output(
                ["git", "-C", str(root), "rev-parse", "HEAD"],
                stderr=subprocess.DEVNULL,
                timeout=3,
                text=True,
            ).strip()
            if re.fullmatch("[0-9a-f]{40}", value):
                commit = value
        except (OSError, subprocess.SubprocessError):
            pass
    versions: dict[str, str | None] = {}
    for name in ("adzuki-gwas-analysis", "numpy", "pandas", "scipy", "matplotlib"):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return {
        "python": platform.python_version(),
        "packages": versions,
        "source_sha256": source.hexdigest(),
        "git_commit": commit,
        "lock_sha256": lock_hash,
        "code_identity": "source_sha256 is authoritative",
    }


def input_checksums(inputs: dict[str, Path]) -> dict[str, str]:
    """Persist role names and hashes only, never private local file paths."""
    return {role: compute_sha256(path) for role, path in sorted(inputs.items())}


def artifact_checksums(root: Path) -> dict[str, str]:
    result = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError("artifact symlink is not allowed")
        if path.is_file() and path not in (root / PROVENANCE_FILE, root / CONTRACT_FILE):
            result[path.relative_to(root).as_posix()] = compute_sha256(path)
    return result


def finish_provenance(
    stage: Path,
    *,
    operation: str,
    inputs: dict[str, Path],
    initial: dict[str, str],
    parameters: dict[str, Any],
    environment: dict[str, Any],
    families: list[dict[str, Any]],
) -> None:
    if input_checksums(inputs) != initial:
        raise ValueError("analysis input changed during execution")
    identity = {
        "operation": operation,
        "inputs": initial,
        "parameters": parameters,
        "environment": environment,
        "families": families,
    }
    document = {
        "schema_version": 1,
        "run_id": hashlib.sha256(json_bytes(identity)).hexdigest(),
        **identity,
        "artifacts": artifact_checksums(stage),
    }
    write_json(stage / PROVENANCE_FILE, document)
    write_json(
        stage / CONTRACT_FILE,
        {
            "schema_version": 1,
            "requires_provenance": True,
            "analysis_run_sha256": compute_sha256(stage / PROVENANCE_FILE),
        },
    )


def validate_provenance(root: Path, *, required: bool = False) -> dict[str, Any] | None:
    manifest, contract = root / PROVENANCE_FILE, root / CONTRACT_FILE
    if (
        not manifest.exists()
        and not contract.exists()
        and not required
        and not manifest.is_symlink()
        and not contract.is_symlink()
    ):
        return None
    if (
        manifest.is_symlink()
        or contract.is_symlink()
        or not manifest.is_file()
        or not contract.is_file()
    ):
        raise ValueError("incomplete analysis-generation provenance")
    try:
        record, marker = json.loads(manifest.read_bytes()), json.loads(contract.read_bytes())
        if (
            not isinstance(record, dict)
            or not isinstance(marker, dict)
            or type(record.get("schema_version")) is not int
            or record.get("schema_version") != 1
            or type(marker.get("schema_version")) is not int
            or marker.get("schema_version") != 1
            or marker.get("requires_provenance") is not True
        ):
            raise ValueError("unsupported analysis provenance schema")
        if marker.get("analysis_run_sha256") != compute_sha256(manifest):
            raise ValueError("analysis-generation manifest checksum mismatch")
        keys = ("operation", "inputs", "parameters", "environment", "families")
        identity = {key: record[key] for key in keys}
        if record["run_id"] != hashlib.sha256(json_bytes(identity)).hexdigest():
            raise ValueError("analysis run identity mismatch")
        if not isinstance(record["environment"], dict) or not isinstance(record["inputs"], dict):
            raise ValueError("analysis-generation metadata is incomplete")
        if record["artifacts"] != artifact_checksums(root):
            raise ValueError("analysis artifact checksum/inventory mismatch")
        return record
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("invalid analysis-generation metadata") from exc


def run_audited_analysis(
    *,
    command: str,
    manifest_path: Path,
    data_dir: Path,
    output_dir: Path,
    dataset_id: str | None,
    alpha: float,
    fdr_level: float,
    clustering_distance: int | None,
    threshold: float = 1e-5,
) -> None:
    from adzuki_gwas_analysis.analysis.batch import run_batch
    from adzuki_gwas_analysis.analysis.pipeline import run_candidates
    from adzuki_gwas_analysis.manifest import load_manifest

    manifest = load_manifest(manifest_path)
    if command not in ("batch", "candidates"):
        raise ValueError("unsupported audited analysis command")
    entries = manifest.datasets if command == "batch" else (manifest.get(str(dataset_id)),)
    inputs = {
        "manifest": manifest_path,
        **{entry.dataset_id: data_dir / entry.member_filename for entry in entries},
    }
    initial, environment = input_checksums(inputs), generation_environment()
    parameters = {
        "alpha": alpha,
        "fdr_level": fdr_level,
        "clustering_distance": clustering_distance,
        "threshold": threshold,
        "clustering_rule": "physical distance; not LD or independently defined QTL",
    }
    families = [
        {
            "dataset_id": e.dataset_id,
            "reference": e.reference,
            "trait": e.trait,
            "test": "LRT",
            "pvalue_column": manifest.pvalue_columns.primary,
            "source_sha256": e.member_sha256,
        }
        for e in entries
    ]
    with output_transaction(output_dir) as stage:
        if command == "batch":
            run_batch(
                manifest_path=manifest_path,
                data_dir=data_dir,
                output_dir=stage,
                alpha=alpha,
                fdr_level=fdr_level,
                threshold=threshold,
                clustering_distance=clustering_distance,
            )
        else:
            if clustering_distance is None:
                raise ValueError("clustering_distance is required")
            run_candidates(
                manifest_path=manifest_path,
                data_dir=data_dir,
                dataset_id=str(dataset_id),
                output_dir=stage,
                alpha=alpha,
                fdr_level=fdr_level,
                clustering_distance=clustering_distance,
            )
        finish_provenance(
            stage,
            operation="dryad." + command,
            inputs=inputs,
            initial=initial,
            parameters=parameters,
            environment=environment,
            families=families,
        )
        validate_provenance(stage, required=True)
