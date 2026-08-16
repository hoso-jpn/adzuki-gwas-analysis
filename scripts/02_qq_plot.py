"""Backward-compatible wrapper: QQ plot for Miyagi water permeability.

Thin wrapper around :func:`adzuki_gwas_analysis.analysis.pipeline.run_qq`
(see Issue #3). Validates the dataset against the schema v1 manifest before
plotting; produces no output if validation fails -- unlike the original
version of this script, out-of-range p-values are never silently dropped.
Equivalent to::

    uv run adzuki-gwas-analyze qq --output-dir plots
"""

from __future__ import annotations

from pathlib import Path

from adzuki_gwas_analysis.analysis.pipeline import run_qq

DATASET_ID = "miyagi_water_permeability"
MANIFEST_PATH = Path("manifest.toml")
DATA_DIR = Path("data/raw")
OUTPUT_PATH = Path("plots/water_permeability_qq.png")


def main() -> None:
    run_qq(
        manifest_path=MANIFEST_PATH,
        data_dir=DATA_DIR,
        dataset_id=DATASET_ID,
        output_path=OUTPUT_PATH,
    )
    print("done")


if __name__ == "__main__":
    main()
