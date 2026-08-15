"""Command-line entry point for validating GWAS summary-statistics files.

Usage::

    adzuki-gwas-validate --manifest manifest.toml --data-dir data/raw
    adzuki-gwas-validate --manifest manifest.toml --data-dir data/raw --output result.json

Exits non-zero if the manifest fails to load, or if any dataset fails
validation.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from adzuki_gwas_analysis.errors import ManifestError
from adzuki_gwas_analysis.manifest import load_manifest
from adzuki_gwas_analysis.validate import ValidationResult, validate_all


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="adzuki-gwas-validate",
        description=("Validate GWAS summary-statistics files against the input-contract manifest."),
    )
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="Path to the manifest TOML file",
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        type=Path,
        help="Directory containing the extracted .assoc.txt files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional path to write a machine-readable JSON validation summary",
    )
    return parser


def _write_summary(output_path: Path, results: list[ValidationResult]) -> None:
    payload = [asdict(result) for result in results]
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    """Run the validation CLI; returns the process exit code (0 = success)."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        manifest = load_manifest(args.manifest)
    except ManifestError as exc:
        print(f"ERROR: failed to load manifest {args.manifest}: {exc}", file=sys.stderr)
        return 1

    results = validate_all(manifest.datasets, args.data_dir)

    for result in results:
        status = "OK" if result.success else "FAIL"
        print(
            f"[{status}] dataset={result.dataset_id} reference={result.reference} "
            f"trait={result.trait} filename={result.filename} "
            f"row_count={result.row_count} sha256={result.sha256}"
        )
        if not result.success:
            print(f"        error: {result.error}", file=sys.stderr)

    if args.output is not None:
        _write_summary(args.output, results)

    failed = [r for r in results if not r.success]
    if failed:
        print(
            f"\n{len(failed)} of {len(results)} dataset(s) failed validation.",
            file=sys.stderr,
        )
        return 1

    print(f"\nAll {len(results)} dataset(s) passed validation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
