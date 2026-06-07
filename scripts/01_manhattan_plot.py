import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILE = "data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"

df = pd.read_csv(FILE, sep="\t")

df["minuslog10p"] = -np.log10(df["pval"])

chroms = sorted(df["chr"].unique())

x = []
xticks = []
xticklabels = []

offset = 0

for chrom in chroms:
    sub = df[df["chr"] == chrom].copy()

    sub["x"] = sub["pos"] + offset

    plt.scatter(
        sub["x"],
        sub["minuslog10p"],
        s=2,
        alpha=0.6
    )

    xticks.append(
        sub["x"].median()
    )

    xticklabels.append(
        chrom.replace("Chr","")
    )

    offset += sub["pos"].max()

plt.axhline(
    -np.log10(1e-5),
    linestyle="--"
)

plt.xticks(
    xticks,
    xticklabels,
    rotation=0
)

plt.xlabel("Chromosome")
plt.ylabel("-log10(p)")
plt.title("Water Permeability GWAS")

plt.tight_layout()

plt.savefig(
    "plots/water_permeability_manhattan.png",
    dpi=300
)

print("done")
