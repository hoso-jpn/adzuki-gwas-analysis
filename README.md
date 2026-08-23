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
writes one kind of output under `--output-dir` (or to the exact path passed to `--output`).
None of the examples below write directly into the tracked `plots/`/`results/` directories --
see "Legacy scripts" below for the wrappers that do:

```bash
INDIVIDUAL_OUTPUT_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze manhattan --output-dir "$INDIVIDUAL_OUTPUT_DIR"
uv run adzuki-gwas-analyze qq --output-dir "$INDIVIDUAL_OUTPUT_DIR"
uv run adzuki-gwas-analyze regions --output-dir "$INDIVIDUAL_OUTPUT_DIR"
uv run adzuki-gwas-analyze top-variants \
  --output "$INDIVIDUAL_OUTPUT_DIR/top_variants_by_region.tsv"

find "$INDIVIDUAL_OUTPUT_DIR" -maxdepth 1 -type f -print
```

`adzuki-gwas-analyze all` is a **single-output-directory bundle command**: it validates the
input once and then writes all of the following under one `--output-dir`, in one pass --

- `<dataset_id>_manhattan.png`
- `<dataset_id>_qq.png`
- one regional PNG per region in `--regions-config`, named by that region's
  `output_filename`
- `top_variants_by_region.tsv`

It is **not** a reproducer of this README's tracked `plots/`/`results/` layout above -- see
"Legacy scripts" below for that. Use a separate scratch directory from the individual
commands above -- reusing the same directory would let `all`'s bundle output overwrite what
those just wrote -- e.g. for a one-pass smoke test of all four outputs together:

```bash
BUNDLE_OUTPUT_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze all --output-dir "$BUNDLE_OUTPUT_DIR"
find "$BUNDLE_OUTPUT_DIR" -maxdepth 1 -type f -print
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
genome-wide significance level, and this default is unchanged by
[Issue #7](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/7)'s `diagnostics`
subcommand below. Bonferroni/BH-FDR/lambda_GC are computed by `adzuki-gwas-analyze
diagnostics` as a separate command with its own `--alpha`/`--fdr-level` flags -- see
[Statistical Diagnostics](#statistical-diagnostics-bonferroni--bh-fdr--genomic-inflation-factor)
below.

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

This is a correctness check, not a performance benchmark. Wall-time and peak-memory
figures reported anywhere in this repository's history (Issue/PR descriptions) are
measured on whichever specific machine ran that specific verification -- recorded only to
confirm the analysis fits comfortably in memory (one dataset processed at a time) -- and
are never a performance guarantee or SLA for any other environment.

---

## Statistical Diagnostics (Bonferroni / BH-FDR / genomic inflation factor)

[Issue #7](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/7) adds
`adzuki-gwas-analyze diagnostics`: Bonferroni family-wise error rate (FWER) correction,
Benjamini-Hochberg false discovery rate (FDR) correction, and the genomic inflation factor
(lambda_GC), computed against one dataset's validated `pval` column. This is a **post-hoc
diagnostic pass over already-published GWAS summary statistics -- not a GWAS re-run**, and
it does not re-correct population structure, kinship, or batch effects.

```bash
DIAGNOSTICS_OUTPUT_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze diagnostics --output-dir "$DIAGNOSTICS_OUTPUT_DIR"
find "$DIAGNOSTICS_OUTPUT_DIR" -maxdepth 1 -type f -print
```

Writes exactly two files, always: `statistical_diagnostics.tsv` (one summary row) and
`significant_variants.tsv` (the union of Bonferroni/BH discoveries, in original input row
order -- header-only with zero rows when there are none). `--alpha` (default `0.05`) and
`--fdr-level` (default `0.05`) are independent from each other and from the legacy
`--threshold` above. This subcommand does not change what `manhattan`/`qq`/`regional`/
`regions`/`top-variants`/`all` produce.

### The multiple-testing family

The family is fixed to **one `dataset_id` x the manifest's declared
`pvalue_columns.primary` column (`pval`) x every variant that passed schema v1 validation
for that one file**. The 6 Dryad files are never combined; Miyagi and Shumari are never
combined; the 3 traits are never combined; a post-hoc visualization region is never treated
as its own family; and `p_wald`/`p_score` are never substituted for `pval`. The
implementation reads the corrected column's name from the loaded manifest object, not a
hardcoded `"pval"` literal, and asserts that the number of p-values it loads equals the row
count schema v1 validation actually counted for that file (not `manifest.toml`'s declared
`row_count` taken on faith).

### Bonferroni

Family-wise alpha default `0.05`; threshold `alpha / m` where `m` is the family's test
count; a variant is significant iff `p <= alpha / m` (equivalently
`min(p * m, 1.0) <= alpha`). This is a threshold **for this dataset's own multiple-testing
family only** -- not a universal genome-wide-significance level -- and it is not adjusted
for linkage disequilibrium among SNPs (`m` is the raw variant count, not an
LD-pruned effective-test count), so it may be conservative.

### Benjamini-Hochberg FDR

Computed via `scipy.stats.false_discovery_control(pvalues, method="bh")`. The adjusted
p-value column is named `pval_bh`; it is a BH-adjusted p-value, **not a Storey-style
q-value** (a different estimator this repository does not implement), and this repository
never calls it a "q-value" unqualified. `bh_raw_p_cutoff` is the largest raw `pval` among
rejected variants, or empty when there are zero discoveries.

**Dependency condition:** the original Benjamini & Hochberg (1995) FDR guarantee holds
under independence; Benjamini & Yekutieli (2001) extended it to positive regression
dependency on a subset (PRDS) -- a specific, narrower condition than "any dependence
structure." SNPs in this dataset are correlated through linkage disequilibrium, and this
repository has **not** demonstrated that `pval` here satisfies PRDS. This repository does
not claim FDR control under arbitrary dependence, does not conflate the 1995 (independence)
and 2001 (PRDS) results, and does not implement the Benjamini-Yekutieli correction itself.

### Genomic inflation factor (lambda_GC)

```
chi2_values     = scipy.stats.chi2.isf(pvalues, df=1)
expected_median = scipy.stats.chi2.ppf(0.5, df=1)
lambda_gc       = median(chi2_values) / expected_median
```

Uses the survival-function side (`isf`) to avoid catastrophic cancellation at small
p-values, never substitutes `-2 * log(p)`, never rounds `expected_median` to a fixed
literal, and maps `p == 1` to `chi2 == 0`. lambda_GC is reported **as a diagnostic value
only**: nothing here divides a test statistic or p-value by it, and no
genomic-control-adjusted p-value is produced. lambda_GC alone does not establish or rule out
population stratification, kinship, or batch effects as the cause of any inflation --
polygenicity also inflates GWAS test statistics and is not distinguishable from confounding
by lambda_GC alone (Bulik-Sullivan et al. 2015, "LD Score Regression Distinguishes
Confounding from Polygenicity"). LD Score regression, which could make that distinction, is
out of scope here because this repository has neither individual-level genotypes nor a
validated LD reference.

**Evidence for the `df=1` assumption, stated precisely:** Dryad's own dataset metadata for
this archive states the GWAS was run with GEMMA; Dryad's data dictionary defines `pval` as
the likelihood-ratio-test (LRT) p-value; GEMMA's univariate linear mixed model (whose output
columns match this dataset's schema exactly) tests a single scalar marker effect `beta` via
`H0: beta = 0` against `H1: beta != 0` for each SNP, and `df=1` follows from that
single-free-parameter test formulation. **The publication's own Methods/Supplementary
Methods text has not been checked** (`science.org` returned `403 Forbidden`; no accessible
preprint was found), so `df=1` is an explicit analysis assumption grounded in the confirmed
output contract and the GEMMA model documentation -- not a fact confirmed by reading the
paper itself. If information contradicting `df=1` for this dataset is found later, lambda_GC
must be re-evaluated before being relied upon.

### Real-data measurements are this machine's reference values only

Any wall-time or peak-RSS figures recorded for `diagnostics` (in Issue #7's PR description)
were measured on one specific Linux x86_64 machine at one point in time and are recorded
only to confirm the computation fits comfortably in memory for one dataset -- they are not a
performance benchmark and say nothing about other hardware.

---

## Batch Analysis (all 6 datasets)

[Issue #9](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/9) adds
`adzuki-gwas-analyze batch`: Manhattan + QQ + Bonferroni/BH-FDR/lambda_GC diagnostics for
**every** dataset declared in `manifest.toml`, in that file's own order, in one command.
This is the standard, auditable artifact set this repository produces per public GWAS
dataset. Candidate-SNP extraction on top of it is implemented (opt-in via
`--clustering-distance`; see
[Candidate SNP Extraction](#candidate-snp-extraction-association-peaks-and-downstream-validation-priority)
below); flanking-SNP reporting and customer-facing summaries remain future work, **not**
implemented here.

```bash
BATCH_OUTPUT_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze batch --output-dir "$BATCH_OUTPUT_DIR"
find "$BATCH_OUTPUT_DIR" -type f | sort
```

`--output-dir` is required and must not already exist as a non-empty directory (batch
publishes all 25 files as one all-or-nothing unit; see
[Output structure and transaction semantics](#output-structure-and-transaction-semantics)
below). There is no `--dataset-id` -- `batch` always processes all 6 canonical datasets.
`--alpha`/`--fdr-level` (both default `0.05`) and `--threshold` (default `1e-5`, the same
legacy visualization line described above) are applied independently to each dataset.
`batch` does not render regional plots or a top-variant-by-region TSV (both specific to
`miyagi_water_permeability`'s post-hoc regions) and does not change what `manhattan`/`qq`/
`regional`/`regions`/`top-variants`/`all`/`diagnostics` produce.

### One validation, one load, per dataset

For each of the 6 datasets, `batch` validates against schema v1 exactly once, loads its
analysis DataFrame exactly once, and reuses that one DataFrame for its Manhattan plot, QQ
plot, and diagnostics -- never re-validating or re-loading the same file for a second
artifact (unlike calling the standalone `manhattan`/`qq`/`diagnostics` subcommands three
times per dataset, which would do exactly that). No two datasets' DataFrames, p-value
arrays, or adjusted-p-value arrays are ever held in memory at once: only a small
per-dataset result (scalars and output-relative paths) survives after that dataset's 4
files are written, before the next dataset starts.

### Each dataset is its own multiple-testing family

Exactly as for the standalone `diagnostics` subcommand above: each dataset's family is
**that one dataset_id x `pval` x its own schema-v1-validated variants**, computed and
recorded independently of the other 5. Running `batch` never combines the 6 files, never
combines Miyagi and Shumari, never combines the 3 traits, and never treats the sum of all 6
datasets' row counts as a shared correction family. That sum -- **8,187,994** rows across
`miyagi_water_permeability` (1,741,385), `miyagi_red_seedcoat` (1,255,203),
`miyagi_mottled_black_seedcoat` (1,255,203), `shumari_water_permeability` (1,471,837),
`shumari_red_seedcoat` (1,232,183), and `shumari_mottled_black_seedcoat` (1,232,183) -- is
reported in `batch_summary.tsv` only as a total processed-row count, never as `m` for any
Bonferroni/BH computation. Miyagi and Shumari coordinates are never compared against each
other anywhere in `batch`'s output.

### Plot titles reflect trait and reference genome

Each dataset's Manhattan/QQ title is built from its manifest `trait`/`reference` fields
(never by string-splitting `dataset_id`), using this fixed trait-label mapping:

| `trait` | Label |
|---|---|
| `water_permeability` | Water Permeability |
| `red_seedcoat` | Red Seed Coat Color |
| `mottled_black_seedcoat` | Mottled Black Seed Coat Color |

e.g. `miyagi_water_permeability` renders `Water Permeability GWAS (Miyagi reference)` on
its Manhattan plot and `QQ Plot: Water Permeability GWAS (Miyagi reference)` on its QQ
plot. The standalone `manhattan`/`qq` subcommands and the legacy `scripts/01`-`04`
wrappers keep their own existing default titles, unchanged by `batch`.

### Output structure and transaction semantics

```text
<output-dir>/
├── batch_summary.tsv
├── miyagi_water_permeability/       (4 files: _manhattan.png, _qq.png, statistical_diagnostics.tsv, significant_variants.tsv)
├── miyagi_red_seedcoat/             (4 files)
├── miyagi_mottled_black_seedcoat/   (4 files)
├── shumari_water_permeability/      (4 files)
├── shumari_red_seedcoat/            (4 files)
└── shumari_mottled_black_seedcoat/  (4 files)
```

25 files total (6 x 4 + 1). The full tree is built in a staging directory next to
`--output-dir` and only moved into place after all 6 datasets and `batch_summary.tsv` have
succeeded; any failure (a bad `--alpha`/`--fdr-level`/`--threshold`, a failed validation, a
row-count mismatch, a plotting or TSV-writing error) removes the staging directory and
leaves `--output-dir` untouched -- no half-finished batch is ever left looking complete.

`batch_summary.tsv` has one row per dataset, in manifest order, with (at least) the same
per-dataset diagnostic columns as `statistical_diagnostics.tsv` above, plus `reference`,
`trait`, `visualization_threshold`, and the 4 artifact paths (`manhattan_path`, `qq_path`,
`diagnostics_path`, `significant_variants_path`) -- always relative, POSIX-style paths
under `--output-dir`, never an absolute or host-specific path.

### Real-data smoke test: always a scratch directory

```bash
SMOKE_OUTPUT_DIR="$(mktemp -d)"
MPLBACKEND=Agg uv run adzuki-gwas-analyze batch --output-dir "$SMOKE_OUTPUT_DIR"
```

Never point `--output-dir` at `plots/` or `results/` directly. Any wall-time/peak-RSS
figures recorded for a `batch` run are that one run's reference values on that one
machine, not a performance guarantee or SLA for any other environment.

---

## Candidate SNP Extraction (association peaks and downstream validation priority)

[Issue #10](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/10) adds
`adzuki-gwas-analyze candidates`: it clusters one dataset's own Bonferroni-or-BH-significant
variants (the same population `diagnostics`/`batch` write to `significant_variants.tsv`)
into physical-distance "signals" and produces a deterministic, auditable priority ordering
for which candidates a breeder should look at first. This is a **post-hoc consumer of
already-corrected summary statistics -- not a new statistical test, not linkage-disequilibrium
(LD) analysis, and not GWAS re-analysis.**

```bash
CANDIDATES_OUTPUT_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze candidates \
    --output-dir "$CANDIDATES_OUTPUT_DIR" \
    --clustering-distance 50000
find "$CANDIDATES_OUTPUT_DIR" -maxdepth 1 -type f -print
```

`candidates` is self-sufficient like `diagnostics`: it validates and loads its dataset
exactly once, computes the same Bonferroni/BH/lambda_GC diagnostics, and writes 5 files --
`statistical_diagnostics.tsv`, `significant_variants.tsv`, `association_peaks.tsv`,
`candidate_snps.tsv`, and `candidate_ranking.tsv` -- always consistent with the same
invocation's own `--alpha`/`--fdr-level`, never a stale file read back from a previous,
possibly differently-parameterized `diagnostics` run.

`adzuki-gwas-analyze batch --clustering-distance N` additively writes the same 3 candidate
files into every one of the 6 dataset directories (43 files total instead of 25).
**Omitting `--clustering-distance` from `batch` leaves its original 25-file,
`schema_version=1` output completely unchanged** -- candidate extraction is opt-in and does
not alter any existing subcommand's default behavior.

### Scientific scope -- read this before interpreting any output

This repository has no individual-level genotype data, so **no linkage disequilibrium can
be computed here**. Given that constraint:

- A "signal" is a **physical-distance cluster only** -- adjacent significant variants on
  the same chromosome, in the same dataset, within `--clustering-distance` base pairs of
  some other already-clustered variant in the same signal (a chained/interval merge, not
  merely within that distance of the signal's first variant). It is **not** an LD block and
  **not** an independently defined QTL interval.
- A signal's **lead variant** is whichever variant this repository's own deterministic
  tie-break rule (below) selects first. It is **not** asserted to be the causal variant, and
  no other variant in the same signal is asserted to be merely a proxy for it.
- `priority_tier`/`priority_reasons` describe **downstream validation priority** --
  which candidates to look at first -- built only from already-computed, already-explained
  quantities (this dataset's own Bonferroni/BH significance flags). This is **not** a
  biological-importance ranking, **not** a probability of being a true positive, and
  **not** a validated or experimentally confirmed breeding marker.
- `dataset_id`/`reference`/`trait` are carried on every output row. A signal never spans
  more than one chromosome or more than one dataset, and Miyagi/Shumari coordinates are
  never combined, compared, or jointly clustered -- exactly as for `diagnostics`/`batch`
  above.

### `--clustering-distance` has no default

There is no scientifically justified default clustering window anywhere in this
repository. The existing `config/*_regions.toml` windows
([Post-hoc visualization regions](#post-hoc-visualization-regions) above) are human-picked,
post-hoc visualization windows chosen by eye from a Manhattan plot -- not an LD or
QTL-interval estimate -- and this repository has no individual-level genotypes from which
an LD-based window could be derived. `--clustering-distance` is therefore **required on
both `candidates` and, when opted into, `batch`**, and this repository never silently
assumes a value. Any value used in this README's examples (e.g. `50000`) is illustrative
only, not a recommendation.

### Clustering rule

Within one dataset's significant-variant population, grouped by chromosome (in natural
numeric order -- `Chr2` before `Chr10` -- never lexicographic order) and sorted by
position: a new signal starts whenever the gap between a variant's position and the running
maximum position of the currently open signal exceeds `--clustering-distance`. This is a
**chained merge** -- a run of variants each within the distance of some other
already-clustered neighbor can span more than `--clustering-distance` end-to-end -- recorded
verbatim in every `association_peaks.tsv` row's `clustering_method` column so the file is
self-describing without cross-referencing this section.

### Lead variant and its tie-break

The lead variant of a signal is chosen by, in order: (1) smallest primary `pval`; (2) on a
tie, the larger-magnitude `beta` (a real, always-finite effect-size estimate schema v1
already guarantees for every row -- see the effect-size note above); (3) on a further tie,
the smallest `pos`; (4) on a full tie (e.g. two multi-allelic records at the same position
with identical `pval` and `|beta|`), the row that appeared first in the source
`.assoc.txt` file. This deliberately differs from the pre-existing
[regional top-variant selection](#regional-plots) (`select_top_variant`), which breaks ties
purely by first-occurrence file order with no candidate-prioritization claim attached; a
signal's lead variant here is presented as the headline candidate for that signal, so an
effect-size-based tie-break is used instead.

### Priority tier and reasons

`priority_tier` is `1` for Bonferroni-significant candidates and `2` for candidates
significant under Benjamini-Hochberg FDR only -- derived purely from the two boolean flags
`diagnostics` already computes, with no additional weighting, scoring formula, or hidden
parameter. `candidate_rank` orders the whole dataset's candidate population by
`priority_tier` first, then the same tie-break used for lead-variant selection.
`priority_reasons` is a semicolon-joined, fully explainable list (e.g.
`bonferroni_significant;lead_variant_of_signal;member_of_multi_variant_signal`) --
never a single opaque score presented as if it were a probability of biological
importance.

### Output structure

```text
<output-dir>/
├── statistical_diagnostics.tsv
├── significant_variants.tsv
├── association_peaks.tsv     (1 row per signal)
├── candidate_snps.tsv        (1 row per significant variant, genome-ordered)
└── candidate_ranking.tsv     (1 row per significant variant, priority-ordered)
```

`association_peaks.tsv` columns: `schema_version`, `dataset_id`, `reference`, `trait`,
`signal_id`, `chromosome`, `start`, `end`, `n_significant_variants`, `lead_pos`,
`lead_allele1`, `lead_allele0`, `lead_pval`, `lead_pval_bonferroni`, `lead_pval_bh`,
`clustering_method`, `clustering_distance`.

`candidate_snps.tsv` columns: `schema_version`, `dataset_id`, `reference`, `trait`,
`signal_id`, `chr`, `pos`, `allele1`, `allele0`, `af`, `beta`, `pval`, `pval_bonferroni`,
`pval_bh`, `bonferroni_significant`, `bh_significant`, `is_lead_variant`.

`candidate_ranking.tsv` columns: `schema_version`, `dataset_id`, `reference`, `trait`,
`signal_id`, `chr`, `pos`, `allele1`, `allele0`, `candidate_rank`, `priority_tier`,
`priority_reasons`.

All 3 files are always produced, header-only with zero rows (never an omitted file) when a
dataset has zero significant variants -- matching `significant_variants.tsv`'s existing
convention.

### `batch`'s opt-in schema

`run_batch(..., clustering_distance=None)` (the default) produces exactly the same 25-file
tree and `batch_summary.tsv` `schema_version=1` this repository has always produced.
Supplying `clustering_distance` bumps `batch_summary.tsv` to `schema_version=2`, which adds
`clustering_distance`, `n_signals`, `n_candidates`, `association_peaks_path`,
`candidate_snps_path`, and `candidate_ranking_path` columns -- computed and reported
independently per dataset, never as a shared 6-dataset ranking population.

### What is explicitly out of scope here

Reference-genome coordinate/sequence-asset contracts, flanking-sequence and
neighboring-variant extraction, and ARMS marker candidate and primer design are separate,
already-tracked follow-up issues (see [Related Repositories](#related-repositories) and
this repository's issue tracker) -- **none of them are implemented by this section.**
Customer-facing report generation on top of this section's output is implemented; see
[Customer Report and Audit Package](#customer-report-and-audit-package) below.

---

## Customer Report and Audit Package

[Issue #11](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/11) adds
`adzuki-gwas-analyze report`: a **consumer**, not a new analysis, that turns an existing
candidate-enabled `batch` output into a self-contained delivery package -- a short,
non-specialist executive summary, a detailed technical report, and a machine-readable
audit trail. It never re-validates `.assoc.txt` files, never recomputes Bonferroni/BH/
lambda_GC, and never re-clusters candidates; it only reads, cross-checks, and repackages
artifacts `batch` already produced.

```bash
BATCH_OUTPUT_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze batch --output-dir "$BATCH_OUTPUT_DIR" --clustering-distance 50000

DELIVERY_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze report --analysis-dir "$BATCH_OUTPUT_DIR" --output-dir "$DELIVERY_DIR"
find "$DELIVERY_DIR" -maxdepth 2 -print | sort
```

### Input contract: a candidate-enabled `batch` output only

`--analysis-dir` must be a `batch` output whose `batch_summary.tsv` has
`schema_version=2` (i.e. that `batch` run was given an explicit `--clustering-distance` --
see [Candidate SNP Extraction](#candidate-snp-extraction-association-peaks-and-downstream-validation-priority)
above). A `schema_version=1` (no-candidates) `batch` output is **rejected** with an
actionable error ("rerun batch with an explicit --clustering-distance") rather than
silently proceeding with an implicit clustering distance -- `report` never invents a
parameter on the caller's behalf.

### Output structure

```text
<output-dir>/
├── executive_summary.md
├── analysis_report.md
├── artifacts/
│   ├── batch_summary.tsv
│   └── <dataset_id>/   (the same 7 files batch wrote for that dataset, copied verbatim)
└── reproducibility/
    ├── input_checksums.tsv       (dataset_id, reference, trait, source_sha256)
    ├── software_versions.json    (report-generation environment only -- see below)
    └── run_manifest.json         (parameters, per-dataset counts, scientific scope,
                                    and a checksummed inventory of every delivered file)
```

Only an explicit allowlist of derived filenames is ever copied -- never a recursive
directory copy, and **never** the raw `.assoc.txt` files. Every reference inside
`executive_summary.md`/`analysis_report.md` is a path relative to the delivery package
root, so the package remains self-contained if copied to another machine.

### Cross-artifact consistency, checked before anything is written

Before generating any report content, `report` reads `batch_summary.tsv` and every
dataset's `statistical_diagnostics.tsv`/`significant_variants.tsv`/
`association_peaks.tsv`/`candidate_snps.tsv`/`candidate_ranking.tsv`, and confirms they
still agree with each other: `n_tests`/Bonferroni and BH discovery counts/lambda_GC match
between `batch_summary.tsv` and that dataset's own diagnostics file; `n_signals`/
`n_candidates` match the corresponding files' row counts; every file's `dataset_id`/
`reference`/`trait` match the row it came from; `candidate_rank` is a dense `1..n`
sequence; and every `signal_id` a candidate references actually exists in
`association_peaks.tsv`. A dataset with zero candidates (header-only files,
`n_candidates=0`) is a normal, fully-checked state, not an error. Any disagreement aborts
report generation with no output written -- `report` never generates a report that would
silently contradict its own source data.

### `software_versions.json`: what it does and does not claim

The batch/candidate artifacts `report` consumes do not themselves record the software
versions or Git commit that generated them. `software_versions.json` therefore records
**only the environment generating the report itself** (`report_generation_environment`:
Python version, platform, this package's and its runtime dependencies' versions, and a
best-effort local Git commit) and sets `analysis_generation_environment` to the literal
string `"unavailable_from_source_artifacts"` -- it never substitutes the report-generation
environment for the unknown analysis-generation one.

### Confidentiality

The delivery package never contains `os.environ`, a hostname (`platform.node()` is never
called), a username, or an absolute filesystem path. No network access, external
service, or LLM call happens anywhere in `report` -- every piece of report text is
generated from fixed, offline string templates.

### Scientific boundary (unchanged from the sections above)

`report` restates, rather than reinterprets, the constraints already established by
[Statistical Diagnostics](#statistical-diagnostics-bonferroni--bh-fdr--genomic-inflation-factor)
and [Candidate SNP Extraction](#candidate-snp-extraction-association-peaks-and-downstream-validation-priority):
this is a post-hoc re-analysis of already-published summary statistics (the GWAS itself
was not re-run); association signals are physical-distance clusters only, never LD blocks
or independently established QTL intervals; lead variants and candidates are never
asserted to be causal or validated breeding markers; `priority_tier` is a downstream
validation priority, never a biological-importance ranking; and Miyagi/Shumari
coordinates, or any two datasets' candidates, are never combined into one shared ranking.

### Atomic publication

Same staging-directory-then-`os.replace` transaction pattern as `batch` (the check itself
is shared, in `analysis/output_safety.py`): the full package is built in a staging
directory next to `--output-dir` and only published after every step succeeds. Any
failure -- an invalid `--analysis-dir`, a path-safety violation, a cross-artifact
inconsistency, or an I/O error partway through -- removes the staging directory and
leaves `--output-dir` untouched.

### Real-data smoke test: always a scratch directory

```bash
SMOKE_ANALYSIS_DIR="$(mktemp -d)"
MPLBACKEND=Agg uv run adzuki-gwas-analyze batch \
    --output-dir "$SMOKE_ANALYSIS_DIR" --clustering-distance 50000

SMOKE_DELIVERY_DIR="$(mktemp -d)"
uv run adzuki-gwas-analyze report \
    --analysis-dir "$SMOKE_ANALYSIS_DIR" --output-dir "$SMOKE_DELIVERY_DIR"
```

As with `batch`, `50000` above is this smoke test's own illustrative parameter, not a
recommended default -- there is none in this repository. Never point either
`--output-dir` at `plots/` or `results/` directly.

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
