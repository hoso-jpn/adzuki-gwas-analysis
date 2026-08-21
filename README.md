# adzuki-gwas-analysis

Re-analysis and visualization of publicly available GWAS summary statistics for adzuki bean (*Vigna angularis*) water permeability.

This repository provides a reproducible workflow for:
- Visualizing GWAS summary statistics
- Generating Manhattan and QQ plots
- Extracting top associated variants
- Creating regional association plots around major signals

The analysis is based on publicly available data from Chien et al. 2025 and Dryad.

---

## Scope

This repository focuses on reproducible re-analysis and visualization of public GWAS summary statistics.

It does **not** claim novel QTL discovery.

---

## Data Source

| Item | Detail |
|---|---|
| Publication | Chien et al. 2025, *Science* 388: eads2871 |
| DOI | https://doi.org/10.1126/science.ads2871 |
| Dryad | https://datadryad.org/dataset/doi:10.5061/dryad.8w9ghx3xv |
| Trait | Water permeability |

---

## Results

### Manhattan Plot

![Water permeability Manhattan plot](plots/water_permeability_manhattan.png)

### QQ Plot

![Water permeability QQ plot](plots/water_permeability_qq.png)


### Top Variants by Region
> Note: Regional windows were defined post-hoc based on visual inspection of the Manhattan plot and are intended for visualization purposes only. They should not be interpreted as independently defined QTL intervals.

| Region | Chr | Position | Beta | p-value |
|---|---|---|---|---|
| Chr07 5-7 Mb | Chr07 | 6,112,438 | 0.205 | 1.95e-14 |
| Chr05 0.5-1.5 Mb | Chr05 | 773,719 | 0.224 | 6.71e-13 |
| Chr11 7-17 Mb | Chr11 | 16,588,878 | 0.148 | 3.82e-10 |
| Chr07 32.0-33.5 Mb | Chr07 | 32,805,358 | 0.131 | 6.93e-09 |
| Chr09 27-30 Mb | Chr09 | 28,421,501 | 0.098 | 3.03e-08 |

Note: Beta values represent effect size estimates from the linear mixed model as reported in the original GWAS summary statistics. The phenotypic scale and exact interpretation follow the original dataset definition.

### Regional Plots

### Chr07: 5–7 Mb

![Chr07 5-7 Mb regional plot](plots/water_permeability_Chr07_5_7Mb_regional.png)

### Chr07: 32.0–33.5 Mb

![Chr07 32.0-33.5 Mb regional plot](plots/water_permeability_Chr07_32_33_5Mb_regional.png)

### Chr09: 27–30 Mb

![Chr09 27-30 Mb regional plot](plots/water_permeability_Chr09_27_30Mb_regional.png)

### Chr05: 0.5–1.5 Mb

![Chr05 0.5-1.5 Mb regional plot](plots/water_permeability_Chr05_0_5_1_5Mb_regional.png)

### Chr11: 7–17 Mb

![Chr11 7-17 Mb regional plot](plots/water_permeability_Chr11_7_17Mb_regional.png)

---

## Reproducibility

All commands below reproduce the plots and tables shown in this README for the
`miyagi_water_permeability` dataset only. They do **not** cover the other 5 datasets in the
Dryad archive (Miyagi/Shumari x water_permeability/red_seedcoat/mottled_black_seedcoat) --
those are out of scope for this repository's current plots/results (see
[Issue #3](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/3)).

### Validation is mandatory, not optional

Every analysis command below -- both the unified CLI and the legacy per-script wrappers --
runs schema v1 validation (`validate_dataset()`, [Issue #1](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/1)/[#2](https://github.com/hoso-jpn/adzuki-gwas-analysis/pull/2))
against the input file **before** producing any plot or TSV. If validation fails, no output is
written -- out-of-range or malformed rows are never silently dropped or ignored.

### Environment (uv only)

```bash
uv sync --locked
```

`pandas`, `numpy`, and `matplotlib` are declared as regular runtime dependencies in
`pyproject.toml`/`uv.lock`; no separate conda environment is needed.

### Input Data

Download from Dryad and place the extracted `.assoc.txt` files under `data/raw/` (never
committed to this repository -- see `.gitignore`):
```
data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt
```

### Validate first

```bash
uv run adzuki-gwas-validate --manifest manifest.toml --data-dir data/raw
```

### Unified CLI: `adzuki-gwas-analyze`

Each of `manhattan`, `qq`, `regions`, and `top-variants` runs schema v1 validation and then
writes one kind of output under `--output-dir` (or to the exact path passed to `--output`):

```bash
uv run adzuki-gwas-analyze manhattan --output-dir plots
uv run adzuki-gwas-analyze qq --output-dir plots
uv run adzuki-gwas-analyze regions --output-dir plots
uv run adzuki-gwas-analyze top-variants --output results/water_permeability/top_variants_by_region.tsv
```

`adzuki-gwas-analyze all` is a **single-output-directory bundle command**: it validates the
input once and then writes all of the following under one `--output-dir`, in one pass --

- `<dataset_id>_manhattan.png`
- `<dataset_id>_qq.png`
- one regional PNG per region in `--regions-config`, named by that region's
  `output_filename`
- `top_variants_by_region.tsv`

It is **not** a reproducer of this README's tracked `plots/`/`results/` layout above -- see
"Legacy scripts" below for that. Use it by pointing `--output-dir` at a scratch directory
outside anything tracked by git -- e.g. for a one-pass smoke test of all four outputs
together:

```bash
OUTPUT_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze all --output-dir "$OUTPUT_DIR"
find "$OUTPUT_DIR" -maxdepth 1 -type f -print
```

Common flags: `--manifest` (default `manifest.toml`), `--data-dir` (default `data/raw`),
`--dataset-id` (default `miyagi_water_permeability`), `--regions-config` (default
`config/water_permeability_regions.toml`), `--threshold` (default `1e-5`, see below).

**Pointing `--output-dir` directly at `plots` (with the default `--regions-config`):**
schema validation still runs first, as with every command above, but no tracked-artifact
equivalence check or overwrite confirmation runs before anything is written. Concretely, for
`all` or `regions`:

- the 5 regional PNGs are named identically to the tracked files in
  [`config/water_permeability_regions.toml`](config/water_permeability_regions.toml) and are
  replaced in place;
- `all`'s `<dataset_id>_manhattan.png` / `<dataset_id>_qq.png` are added as new,
  differently-named files -- the tracked `water_permeability_manhattan.png` /
  `water_permeability_qq.png` are left untouched;
- `all`'s `top_variants_by_region.tsv` is added directly under `plots/`, not written to
  `results/water_permeability/top_variants_by_region.tsv`.

### Legacy scripts (kept as thin backward-compatible wrappers)

| Script | Description | Input | Output |
|---|---|---|---|
| `01_manhattan_plot.py` | Genome-wide Manhattan plot for water permeability GWAS. | `data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt` | `plots/water_permeability_manhattan.png` |
| `02_qq_plot.py` | QQ plot for GWAS quality-control visualization. | (same) | `plots/water_permeability_qq.png` |
| `03_regional_plot.py` | Regional association plot for a user-specified chromosome interval. | GWAS summary statistics, chromosome, start/end positions | Regional plot PNG |
| `04_extract_top_variants_by_region.py` | Extracts the top variant within each configured visualization window. | GWAS summary statistics + `config/water_permeability_regions.toml` | `results/water_permeability/top_variants_by_region.tsv` |

These scripts are thin wrappers around `src/adzuki_gwas_analysis/analysis/` -- they run the
same mandatory validation and produce the same output as the unified CLI, and take the exact
same command-line arguments they always have. Unlike `adzuki-gwas-analyze all`/`regions` (a
single-output-directory bundle command meant for a scratch directory -- see above), these
wrappers are the intended path for regenerating the tracked
`plots`/`results/water_permeability/` artifacts under their existing names and layout.

Before running them for that purpose, confirm a clean working tree (`git status`); after
running, review what changed before committing -- e.g. `git diff --stat -- plots results`
for the PNGs, and a byte-for-byte check such as `cmp` against the pre-existing file for
`results/water_permeability/top_variants_by_region.tsv` -- rather than committing
regenerated output unreviewed:

```bash
python scripts/01_manhattan_plot.py
python scripts/02_qq_plot.py
python scripts/03_regional_plot.py \
  --input data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt \
  --chrom Chr07 \
  --start 5000000 \
  --end 7000000 \
  --output plots/water_permeability_Chr07_5_7Mb_regional.png \
  --title "Water Permeability GWAS: Chr07 5-7 Mb"

python scripts/03_regional_plot.py \
  --input data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt \
  --chrom Chr07 \
  --start 32000000 \
  --end 33500000 \
  --output plots/water_permeability_Chr07_32_33_5Mb_regional.png \
  --title "Water Permeability GWAS: Chr07 32.0-33.5 Mb"

python scripts/03_regional_plot.py \
  --input data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt \
  --chrom Chr09 \
  --start 27000000 \
  --end 30000000 \
  --output plots/water_permeability_Chr09_27_30Mb_regional.png \
  --title "Water Permeability GWAS: Chr09 27-30 Mb"

python scripts/03_regional_plot.py \
  --input data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt \
  --chrom Chr05 \
  --start 500000 \
  --end 1500000 \
  --output plots/water_permeability_Chr05_0_5_1_5Mb_regional.png \
  --title "Water Permeability GWAS: Chr05 0.5-1.5 Mb"

python scripts/03_regional_plot.py \
  --input data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt \
  --chrom Chr11 \
  --start 7000000 \
  --end 17000000 \
  --output plots/water_permeability_Chr11_7_17Mb_regional.png \
  --title "Water Permeability GWAS: Chr11 7-17 Mb"

python scripts/04_extract_top_variants_by_region.py
```

### Post-hoc visualization regions

The 5 chromosome windows used above (and defined in
[`config/water_permeability_regions.toml`](config/water_permeability_regions.toml)) were chosen
post-hoc by visual inspection of the genome-wide Manhattan plot. They are **not** independently
defined QTL intervals, **not** LD blocks, and not the output of any linkage or fine-mapping
analysis -- treat them as visualization conveniences only.

### The `1e-5` threshold line

The dashed line drawn on Manhattan and regional plots at `1e-5` (configurable via
`--threshold`) is a **legacy visualization threshold**, not a Bonferroni-corrected or
genome-wide significance level. No multiple-testing correction is computed by this repository
(see Issue #3's out-of-scope list); a statistically corrected threshold is a candidate for a
future Issue.

### What this repository is (and is not)

This repository visualizes and re-analyzes existing, publicly deposited GWAS summary
statistics -- it does **not** re-run the underlying GWAS. `pval` is always the
likelihood-ratio-test (LRT) p-value (see [`docs/gwas_input_contract.md`](docs/gwas_input_contract.md)).
Miyagi, Shumari, and Longxiaodou 4 (used by
[adzuki-snp-pipeline](https://github.com/hoso-jpn/adzuki-snp-pipeline)) are three distinct,
non-interchangeable reference coordinate systems -- positions from one must never be
interpreted against another.

### `results/water_permeability/top_snps.tsv`

This file's generation rule is **not verified**: no script producing it exists anywhere in
this repository's git history, and this Issue does not attempt to reconstruct or regenerate
it. It is out of scope here and tracked as a separate, unresolved question.

### Smoke-testing against real data (developers)

Smoke-test `adzuki-gwas-analyze all` against the real `data/raw/` file by writing to a
scratch directory -- never directly to `plots/` or `results/` -- and then compare the
result against the tracked files:

```bash
OUTPUT_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze all --output-dir "$OUTPUT_DIR"
find "$OUTPUT_DIR" -maxdepth 1 -type f -print
cmp "$OUTPUT_DIR/top_variants_by_region.tsv" \
    results/water_permeability/top_variants_by_region.tsv
```

Compare the regional/Manhattan/QQ PNGs under `$OUTPUT_DIR` against the tracked `plots/`
files by inspection (dimensions and content); this repository does not assert byte-for-byte
PNG equivalence, only TSV equivalence and PNG dimension/content equivalence.

This is a correctness check, not a performance benchmark. Any wall-time or peak-memory
numbers reported for this repository were measured on a single Apple Silicon Mac and are
recorded only to confirm the analysis fits comfortably in memory (one dataset processed at
a time) -- they say nothing about Linux/production performance.

---

## Input Contract and Validation

[Issue #1](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/1) /
[PR #2](https://github.com/hoso-jpn/adzuki-gwas-analysis/pull/2) added a machine-readable
manifest (`manifest.toml`) and a schema validator (`src/adzuki_gwas_analysis/`) covering all
6 GWAS summary-statistics files in the Dryad dataset above (3 traits x 2 reference genomes:
Miyagi, Shumari). See [`docs/gwas_input_contract.md`](docs/gwas_input_contract.md) for the
full contract: the 6-dataset list, why Miyagi/Shumari (and Longxiaodou 4) are not
interchangeable coordinate systems, why `pval` is the likelihood-ratio-test p-value used as
this repository's primary statistic, how to obtain and checksum-verify the raw data, and how
to run the validator. Raw data is never committed to this repository.
[Issue #3](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/3) migrated the
`miyagi_water_permeability` analysis scripts above onto this validator/loader; the other 5
datasets are not yet migrated.

```bash
uv sync --locked
uv run adzuki-gwas-validate --manifest manifest.toml --data-dir data/raw
```

---

## Related Repositories

- [adzuki-snp-pipeline](https://github.com/hoso-jpn/adzuki-snp-pipeline)
- [genomic-prediction-resnet-hybrid](https://github.com/hoso-jpn/genomic-prediction-resnet-hybrid)

---

## Author

**Hoso**
Plant Genetics x Bioinformatics x Physical AI

- GitHub: https://github.com/hoso-jpn
- Researchmap: https://researchmap.jp/hosokawa-yusuke

---

## License

MIT License
