# GWAS summary-statistics input contract

This document describes the input-contract and validation infrastructure added in
[Issue #1](https://github.com/hoso-jpn/adzuki-gwas-analysis/issues/1): a machine-readable
manifest plus a schema validator for the 6 public GWAS summary-statistics files this
repository re-analyzes. It does **not** describe a GWAS re-analysis pipeline, LD analysis,
population-structure correction, or genomic-selection training -- see
[Scope and limitations](#scope-and-limitations) below for why.

## Scope

This repository re-analyzes and visualizes **already-computed public GWAS summary
statistics** (per-variant effect size, standard error, and p-values). It is not a
genotype/phenotype analysis pipeline, and the `src/adzuki_gwas_analysis` package added by
Issue #1 only validates that a summary-statistics file matches its declared input
contract -- it does not perform or re-run any statistical test.

## Scope and limitations

The Dryad dataset backing this repository (see [Data source](#data-source) below) contains
only summary statistics: one row per variant, with an already-fitted effect size,
standard error, and p-value. It does **not** include:

- per-individual genotype calls
- per-individual phenotype measurements
- a sample-to-genotype/phenotype correspondence table

Without these, this repository **cannot** perform: GWAS re-execution, LD (linkage
disequilibrium) analysis, population-structure/kinship correction, or genomic-selection
(GS) model training. Any of those would require the underlying per-individual data, which
is not part of this input contract.

## The 6 datasets

The Dryad archive `adzuki_GWAS_data.zip` contains exactly 6 files: 3 traits (water
permeability, red seedcoat color, mottled black seedcoat color) each mapped to 2
reference genome assemblies (Miyagi, Shumari).

| dataset_id | reference | trait | source filename |
|---|---|---|---|
| `miyagi_water_permeability` | Miyagi | water_permeability | `mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt` |
| `miyagi_red_seedcoat` | Miyagi | red_seedcoat | `mapped_to_Miyagi_red_seedcoat.maf_0.05.assoc.txt` |
| `miyagi_mottled_black_seedcoat` | Miyagi | mottled_black_seedcoat | `mapped_to_Miyagi_mottled_black_seedcoat.maf_0.05.assoc.txt` |
| `shumari_water_permeability` | Shumari | water_permeability | `mapped_to_Shumari_water_permeability.maf_0.05.assoc.txt` |
| `shumari_red_seedcoat` | Shumari | red_seedcoat | `mapped_to_Shumari_red_seedcoat.maf_0.05.assoc.txt` |
| `shumari_mottled_black_seedcoat` | Shumari | mottled_black_seedcoat | `mapped_to_Shumari_mottled_black_seedcoat.maf_0.05.assoc.txt` |

This exact list (reference, trait, and filename) is machine-readable in
[`manifest.toml`](../manifest.toml) at the repository root.

## Miyagi and Shumari are different reference genomes

**Miyagi** (wild-type) and **Shumari** (cultivated) are two distinct adzuki bean reference
genome assemblies -- coordinates in a Miyagi-mapped file and a Shumari-mapped file are
**not interchangeable** without an explicit liftover. Neither is the same reference as
**Longxiaodou 4** (`GCF_016808095.1`), the reference genome used by the sibling repository
[adzuki-snp-pipeline](https://github.com/hoso-jpn/adzuki-snp-pipeline).

**Do not annotate a Longxiaodou 4 GFF (or any Longxiaodou-4-coordinate feature) directly
against `chr`/`pos` values from these files.** A coordinate transformation (liftover)
between assemblies is required first; none is implemented in this repository. Comparing
results across the two references within this dataset (Miyagi vs. Shumari) also requires
care for the same reason, and is explicitly out of scope for Issue #1.

## The `pval` column is the likelihood ratio test (LRT) p-value

Every file has three p-value columns, all confirmed against Dryad's own data dictionary:

| column | meaning |
|---|---|
| `p_wald` | Wald test p-value |
| `pval` | **Likelihood ratio test (LRT) p-value** |
| `p_score` | Score test p-value |

`pval` is the primary p-value used by this repository's existing plotting scripts
(`scripts/01_manhattan_plot.py`, `02_qq_plot.py`, `03_regional_plot.py`,
`04_extract_top_variants_by_region.py`) and remains the primary column for the validator,
for compatibility with that existing code. `p_wald` and `p_score` are validated as required
columns with the same `0 < p <= 1` contract, but are not used as the primary result by
anything in this repository yet.

## Data source

- Publication: Chien et al. 2025, *Science* 388: eads2871
- Dryad dataset: <https://doi.org/10.5061/dryad.8w9ghx3xv> (dataset ID 149675, version
  356599, license CC0-1.0)
- Archive: `adzuki_GWAS_data.zip` (257,167,472 bytes)

The archive is not distributed by this repository (see
[Raw data is never committed](#raw-data-is-never-committed) below). To obtain it:

1. Download `adzuki_GWAS_data.zip` from the Dryad dataset page above (`GET
   /api/v2/versions/356599/files` finds the current per-file download link). Dryad's file
   download API requires an `Authorization: Bearer <token>` header, per its own published
   API documentation. This repository does not use any Dryad API credential; the archive
   used here was obtained through the Dryad web UI instead, without inspecting how the UI
   itself authorizes its downloads (session cookie, signed URL, or otherwise). A direct,
   unauthenticated API/curl download attempt returns `401`/`403`, consistent with the
   documented Bearer-token requirement, not evidence of anything beyond that.
2. Verify the archive's SHA-256 against `[archive].sha256` in `manifest.toml` (see
   [Checksum verification](#checksum-verification) below) before extracting anything.
3. Extract the 6 `.assoc.txt` members into `data/raw/` (already `.gitignore`d).

## Checksum verification

```bash
# Archive checksum (compare against [archive].sha256 in manifest.toml)
shasum -a 256 /path/to/adzuki_GWAS_data.zip

# Per-file checksum after extraction (compare against each
# [[datasets]] entry's member_sha256 in manifest.toml)
shasum -a 256 data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt
```

The validator (see below) performs this same SHA-256 check automatically for all 6 files,
in addition to the full schema/row validation.

## Running the validator

```bash
uv sync --locked
uv run adzuki-gwas-validate --manifest manifest.toml --data-dir data/raw
```

Add `--output result.json` to also write a machine-readable JSON summary (dataset ID,
reference, trait, filename, row count, checksum, and pass/fail per dataset). The CLI exits
non-zero if the manifest fails to load or if any dataset fails validation; every error
message includes the dataset ID, file name, and (for row-level errors) the row number and
column involved.

## Raw data is never committed

`adzuki_GWAS_data.zip`, its 6 extracted `.assoc.txt` members, and
`selection_sweep_scanning.xlsx` (a separate, unrelated Dryad file, also out of scope) are
never committed to this repository. `data/raw/` and `*.assoc.txt` are excluded via
`.gitignore`; verify with:

```bash
git check-ignore -v data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt
```

## What Issue #1 does not change

Issue #1 adds the manifest, validator, and CLI described above. It does **not** modify or
migrate the existing analysis scripts (`scripts/01_manhattan_plot.py` through
`04_extract_top_variants_by_region.py`), and does not change any existing `plots/` or
`results/water_permeability/` artifact. Migrating those scripts to use this validator/loader
is tracked as a follow-up issue. `results/water_permeability/top_snps.tsv` also has a known,
unresolved reproducibility gap (it is not produced by any of the 4 existing scripts) that
is recorded, not fixed, by Issue #1.
