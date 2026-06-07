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


### QQ Plot


### Top Variants by Region

| Region | Chr | Position | Beta | p-value |
|---|---|---|---|---|
| Chr07 5-7 Mb | Chr07 | 6,112,438 | 0.205 | 1.95e-14 |
| Chr05 0.5-1.5 Mb | Chr05 | 773,719 | 0.224 | 6.71e-13 |
| Chr11 7-17 Mb | Chr11 | 16,588,878 | 0.148 | 3.82e-10 |
| Chr07 32.0-33.5 Mb | Chr07 | 32,805,358 | 0.131 | 6.93e-09 |
| Chr09 27-30 Mb | Chr09 | 28,421,501 | 0.098 | 3.03e-08 |

### Regional Plots

**Chr07: 5-7 Mb**

**Chr07: 32.0-33.5 Mb**

**Chr09: 27-30 Mb**

**Chr05: 0.5-1.5 Mb**

**Chr11: 7-17 Mb**

---

## Reproducibility

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

```bash
python scripts/01_manhattan_plot.py
python scripts/02_qq_plot.py
python scripts/03_regional_plot.py --input data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt --chrom Chr07 --start 5000000 --end 7000000 --output plots/water_permeability_Chr07_5_7Mb_regional.png
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
