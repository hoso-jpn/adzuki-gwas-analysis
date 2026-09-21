"""Orchestration for the customer-facing report/audit delivery package (Issue #11).

This module is a **consumer**, not a producer: it never re-validates the original
``.assoc.txt`` files, never recomputes Bonferroni/BH/lambda_GC, and never re-clusters
candidates. It reads a candidate-enabled ``batch`` output directory (validated by
:func:`adzuki_gwas_analysis.analysis.report_validation.validate_analysis_dir`, which is
also this function's sole gate for schema/consistency checks) and publishes a
self-contained delivery package: two Markdown reports, a copy of every dataset's derived
artifacts (never the raw ``.assoc.txt`` files -- an explicit allowlist of filenames is
copied, never a recursive directory copy), and 3 machine-readable reproducibility files.

No network access, no external service, and no LLM call happens anywhere in this module;
every piece of report text is built by
:mod:`adzuki_gwas_analysis.analysis.report_content`'s plain string templates.

Determinism: nothing here embeds a wall-clock generation timestamp. The same
``--analysis-dir`` (and the same report-generation environment) always produces the same
delivery package content, so a checksum comparison between two runs is meaningful and
tests never need to inject or normalize a clock.

Confidentiality: the delivery package never contains ``os.environ``, a hostname (never
``platform.node()``), a username, or any absolute filesystem path -- see
:func:`_report_generation_environment` and
:mod:`adzuki_gwas_analysis.analysis.report_content`'s path-building helpers, which only
ever emit paths relative to the delivery package root.

Transaction semantics mirror :mod:`adzuki_gwas_analysis.analysis.batch`: the full package
is built in a staging directory next to ``--output-dir`` and only moved into place, via one
``os.replace``, after every step succeeds; any failure removes the staging directory and
leaves ``--output-dir`` untouched.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import cast

import pandas as pd

from adzuki_gwas_analysis.analysis.output_safety import check_output_dir_is_safe
from adzuki_gwas_analysis.analysis.report_content import (
    build_analysis_report_markdown,
    build_executive_summary_markdown,
    build_run_manifest,
    build_software_versions,
)
from adzuki_gwas_analysis.analysis.report_validation import (
    ValidatedAnalysisDir,
    validate_analysis_dir,
)
from adzuki_gwas_analysis.loader import compute_sha256
from adzuki_gwas_analysis.provenance import CONTRACT_FILE, PROVENANCE_FILE, validate_provenance

#: Packages whose installed version is recorded in software_versions.json. Deliberately a
#: fixed, curated list -- never a dump of every installed package or of os.environ.
_PACKAGES_TO_REPORT: tuple[str, ...] = (
    "adzuki-gwas-analysis",
    "pandas",
    "numpy",
    "scipy",
    "matplotlib",
)

#: (attribute on DatasetArtifactPaths, artifact role) for each of a dataset's 7 files,
#: copied into artifacts/<dataset_id>/ under its own original filename.
_PER_DATASET_ARTIFACT_ROLES: tuple[tuple[str, str], ...] = (
    ("manhattan", "manhattan_plot"),
    ("qq", "qq_plot"),
    ("diagnostics", "statistical_diagnostics"),
    ("significant_variants", "significant_variants"),
    ("association_peaks", "association_peaks"),
    ("candidate_snps", "candidate_snps"),
    ("candidate_ranking", "candidate_ranking"),
)


@dataclass(frozen=True, slots=True)
class ReportOutcome:
    """Everything produced by one :func:`run_report` call."""

    output_dir: Path
    analysis_dir: Path
    n_datasets: int
    n_signals: int
    n_candidates: int
    executive_summary_path: Path
    analysis_report_path: Path
    run_manifest_path: Path


def _atomic_write_text(content: str, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=output_path.parent, suffix=".tmp")
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        tmp_path.write_text(content, encoding="utf-8")
        os.replace(tmp_path, output_path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


def _atomic_write_json(document: dict[str, object], output_path: Path) -> None:
    _atomic_write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", output_path)


def _copy_artifact(src: Path, dst: Path, *, relative_path: str, role: str) -> dict[str, object]:
    """Copy one allowlisted file into the staging tree and record its checksum.

    ``shutil.copyfile`` (content only, no metadata) rather than ``copy2`` -- this
    repository's confidentiality design never carries source-filesystem metadata
    (timestamps included) into a customer-facing package.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    return {"relative_path": relative_path, "sha256": compute_sha256(dst), "role": role}


def _package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for name in _PACKAGES_TO_REPORT:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def _detect_report_generation_git_commit() -> str | None:
    """Best-effort ``git rev-parse HEAD`` for the checkout running this report, or ``None``.

    Never raises: a wheel install with no ``.git`` directory, a missing ``git`` binary, or
    any subprocess failure all fall back to ``None`` rather than aborting report
    generation over a provenance nicety. This is the commit of the checkout *generating
    this report* -- it is never claimed to be the commit that generated the
    ``--analysis-dir`` artifacts being consumed (see
    :mod:`adzuki_gwas_analysis.analysis.report_content`'s
    ``analysis_generation_environment`` sentinel).
    """
    repo_root: Path | None = None
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            repo_root = parent
            break
    if repo_root is None:
        return None
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    commit = result.stdout.strip()
    return commit or None


def _report_generation_environment() -> dict[str, object]:
    """The environment generating *this report* -- never the environment's hostname/user.

    Deliberately uses ``platform.platform()``/``platform.python_version()`` (OS/arch/
    interpreter facts) and never ``platform.node()`` or any environment-variable dump --
    see this module's docstring's Confidentiality section.
    """
    return {
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "packages": _package_versions(),
        "git_commit": _detect_report_generation_git_commit(),
    }


def _copy_all_artifacts(
    validated: ValidatedAnalysisDir, staging_dir: Path
) -> list[dict[str, object]]:
    artifacts: list[dict[str, object]] = []
    batch_summary_relative = "artifacts/batch_summary.tsv"
    artifacts.append(
        _copy_artifact(
            validated.batch_summary_path,
            staging_dir / batch_summary_relative,
            relative_path=batch_summary_relative,
            role="batch_summary",
        )
    )
    for dataset in validated.datasets:
        for attribute, role in _PER_DATASET_ARTIFACT_ROLES:
            src: Path = getattr(dataset.paths, attribute)
            relative_path = f"artifacts/{dataset.dataset_id}/{src.name}"
            artifacts.append(
                _copy_artifact(
                    src, staging_dir / relative_path, relative_path=relative_path, role=role
                )
            )
    if validated.analysis_generation is not None:
        for name in (PROVENANCE_FILE, CONTRACT_FILE):
            relative_path = "artifacts/" + name
            artifacts.append(
                _copy_artifact(
                    validated.analysis_dir / name,
                    staging_dir / relative_path,
                    relative_path=relative_path,
                    role="analysis_generation_provenance",
                )
            )
        validate_provenance(staging_dir / "artifacts", required=True)
    return artifacts


def _write_dataset_artifact_bundle(
    validated: ValidatedAnalysisDir, staging_dir: Path
) -> list[dict[str, object]]:
    artifacts = _copy_all_artifacts(validated, staging_dir)

    executive_summary_relative = "executive_summary.md"
    analysis_report_relative = "analysis_report.md"
    _atomic_write_text(
        build_executive_summary_markdown(validated), staging_dir / executive_summary_relative
    )
    _atomic_write_text(
        build_analysis_report_markdown(validated), staging_dir / analysis_report_relative
    )
    artifacts.append(
        {
            "relative_path": executive_summary_relative,
            "sha256": compute_sha256(staging_dir / executive_summary_relative),
            "role": "executive_summary",
        }
    )
    artifacts.append(
        {
            "relative_path": analysis_report_relative,
            "sha256": compute_sha256(staging_dir / analysis_report_relative),
            "role": "analysis_report",
        }
    )

    checksums_relative = "reproducibility/input_checksums.tsv"
    checksums_table = pd.DataFrame(
        [
            {
                "dataset_id": d.dataset_id,
                "reference": d.reference,
                "trait": d.trait,
                "source_sha256": d.source_sha256,
            }
            for d in validated.datasets
        ]
    )
    checksums_path = staging_dir / checksums_relative
    checksums_path.parent.mkdir(parents=True, exist_ok=True)
    checksums_path.write_text(checksums_table.to_csv(sep="\t", index=False), encoding="utf-8")
    artifacts.append(
        {
            "relative_path": checksums_relative,
            "sha256": compute_sha256(checksums_path),
            "role": "input_checksums",
        }
    )

    software_versions_relative = "reproducibility/software_versions.json"
    software_versions_doc = build_software_versions(
        report_generation_environment=_report_generation_environment(),
        analysis_generation_environment=(
            cast(dict[str, object], validated.analysis_generation["environment"])
            if validated.analysis_generation is not None
            else None
        ),
    )
    _atomic_write_json(software_versions_doc, staging_dir / software_versions_relative)
    artifacts.append(
        {
            "relative_path": software_versions_relative,
            "sha256": compute_sha256(staging_dir / software_versions_relative),
            "role": "software_versions",
        }
    )
    return artifacts


def run_report(
    *, analysis_dir: Path, output_dir: Path, require_provenance: bool = False
) -> ReportOutcome:
    """Validate ``--analysis-dir``, then publish the full delivery package to ``--output-dir``.

    Every read of ``--analysis-dir`` happens inside
    :func:`~adzuki_gwas_analysis.analysis.report_validation.validate_analysis_dir`, before
    any staging directory is created -- an invalid or inconsistent analysis directory
    raises and creates nothing. The delivery package itself is built entirely inside a
    staging directory next to ``--output-dir`` and published with one final ``os.replace``;
    any failure at any step removes the staging directory and leaves ``--output-dir``
    untouched.
    """
    validated = validate_analysis_dir(Path(analysis_dir), require_provenance=require_provenance)

    output_dir = Path(output_dir)
    check_output_dir_is_safe(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    staging_dir = Path(tempfile.mkdtemp(dir=output_dir.parent, prefix=".report-staging-"))
    try:
        artifacts = _write_dataset_artifact_bundle(validated, staging_dir)

        run_manifest_relative = "reproducibility/run_manifest.json"
        run_manifest_doc = build_run_manifest(
            validated,
            provenance={"software_versions_file": "reproducibility/software_versions.json"},
            artifacts=sorted(artifacts, key=lambda a: str(a["relative_path"])),
        )
        _atomic_write_json(run_manifest_doc, staging_dir / run_manifest_relative)

        os.replace(staging_dir, output_dir)
    except BaseException:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise

    return ReportOutcome(
        output_dir=output_dir,
        analysis_dir=Path(analysis_dir),
        n_datasets=len(validated.datasets),
        n_signals=sum(d.n_signals for d in validated.datasets),
        n_candidates=sum(d.n_candidates for d in validated.datasets),
        executive_summary_path=output_dir / "executive_summary.md",
        analysis_report_path=output_dir / "analysis_report.md",
        run_manifest_path=output_dir / run_manifest_relative,
    )
