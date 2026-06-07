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
### Scripts

| Script | Description | Input | Output |
|---|---|---|---|
| `01_manhattan_plot.py` | Generates a genome-wide Manhattan plot for water permeability GWAS. | GWAS summary statistics | `plots/water_permeability_manhattan.png` |
| `02_qq_plot.py` | Generates a QQ plot for GWAS quality-control visualization. | GWAS summary statistics | `plots/water_permeability_qq.png` |
| `03_regional_plot.py` | Generates regional association plots for user-specified chromosome intervals. | GWAS summary statistics, chromosome, start/end positions | Regional plot PNG files |
| `04_extract_top_variants_by_region.py` | Extracts the top associated variant within each visualization window. | GWAS summary statistics | `results/water_permeability/top_variants_by_region.tsv` |

### Environment

```bash
conda activate bioinfo
conda install pandas numpy matplotlib -y
```

### Input Data

Download from Dryad and place under data/raw/:
```
data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt
```

### Run

The following commands reproduce all plots and summary tables shown in this README.

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
