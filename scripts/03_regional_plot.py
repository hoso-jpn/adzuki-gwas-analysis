"""Backward-compatible wrapper: regional association plot from GWAS summary statistics.

Thin wrapper around
:func:`adzuki_gwas_analysis.analysis.pipeline.run_single_regional` (see
Issue #3). Keeps the exact same ``--input``/``--chrom``/``--start``/
``--end``/``--output``/``--title`` arguments as before, but now looks up
``--input`` in the schema v1 manifest and validates it before plotting;
produces no output if validation fails.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from adzuki_gwas_analysis.analysis.pipeline import run_single_regional
from adzuki_gwas_analysis.errors import GwasContractError
from adzuki_gwas_analysis.manifest import load_manifest

MANIFEST_PATH = Path("manifest.toml")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a regional association plot from GWAS summary statistics."
    )
    parser.add_argument("--input", required=True, help="GWAS association file")
    parser.add_argument("--chrom", required=True, help="Chromosome name, e.g. Chr07")
    parser.add_argument("--start", type=int, required=True, help="Start position")
    parser.add_argument("--end", type=int, required=True, help="End position")
    parser.add_argument("--output", required=True, help="Output PNG file")
    parser.add_argument("--title", default=None, help="Plot title")
    args = parser.parse_args()

    input_path = Path(args.input)
    manifest = load_manifest(MANIFEST_PATH)
    try:
        entry = manifest.get_by_filename(input_path.name)
    except KeyError:
        print(
            f"ERROR: {input_path.name!r} is not a member_filename in {MANIFEST_PATH}; "
            f"expected one of the schema v1 dataset files",
            file=sys.stderr,
        )
        return 1

    try:
        outcome = run_single_regional(
            manifest_path=MANIFEST_PATH,
            data_dir=input_path.parent,
            dataset_id=entry.dataset_id,
            chrom=args.chrom,
            start=args.start,
            end=args.end,
            output_path=Path(args.output),
            title=args.title,
            region_id=f"{args.chrom}_{args.start}_{args.end}",
        )
    except GwasContractError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Saved: {outcome.plot.output_path}")
    print("Top variant in region:")
    top = outcome.top_variant
    print(
        pd.Series(
            {
                "chr": top.chrom,
                "pos": top.pos,
                "allele1": top.allele1,
                "allele0": top.allele0,
                "af": top.af,
                "beta": top.beta,
                "pval": top.pval,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
