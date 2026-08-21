"""Backward-compatible wrapper: top-variant-per-region TSV for Miyagi water permeability.

Thin wrapper around
:func:`adzuki_gwas_analysis.analysis.pipeline.run_top_variants` (see
Issue #3). Regions now come from
``config/water_permeability_regions.toml`` rather than being hardcoded here.
Validates the dataset against the schema v1 manifest before extracting;
produces no output if validation fails -- unlike the original version of
this script, which checked ``pval > 0`` but never ``pval <= 1``.
Equivalent to::

    uv run adzuki-gwas-analyze top-variants \\
        --output results/water_permeability/top_variants_by_region.tsv
"""

from __future__ import annotations

from pathlib import Path

from adzuki_gwas_analysis.analysis.pipeline import run_top_variants

DATASET_ID = "miyagi_water_permeability"
MANIFEST_PATH = Path("manifest.toml")
DATA_DIR = Path("data/raw")
REGIONS_CONFIG_PATH = Path("config/water_permeability_regions.toml")
OUTPUT_PATH = Path("results/water_permeability/top_variants_by_region.tsv")


def main() -> None:
    table = run_top_variants(
        manifest_path=MANIFEST_PATH,
        data_dir=DATA_DIR,
        dataset_id=DATASET_ID,
        regions_config_path=REGIONS_CONFIG_PATH,
        output_path=OUTPUT_PATH,
    )
    pd_display = table.copy()
    pd_display["pval"] = pd_display["pval"].map("{:.3e}".format)
    pd_display["beta"] = pd_display["beta"].map("{:.3e}".format)
    print(pd_display.to_string(index=False))
    print(f"Saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
