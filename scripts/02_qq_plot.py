import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILE = "data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"

df = pd.read_csv(FILE, sep="\t")

p = df["pval"].dropna()
p = p[(p > 0) & (p <= 1)]

observed = -np.log10(np.sort(p))
expected = -np.log10(np.arange(1, len(p) + 1) / (len(p) + 1))

plt.figure(figsize=(6, 6))
plt.scatter(expected, observed, s=3, alpha=0.5)

max_val = max(expected.max(), observed.max())
plt.plot([0, max_val], [0, max_val], linestyle="--")

plt.xlabel("Expected -log10(p)")
plt.ylabel("Observed -log10(p)")
plt.title("QQ Plot: Water Permeability GWAS")

plt.tight_layout()
plt.savefig("plots/water_permeability_qq.png", dpi=300)

print("done")
