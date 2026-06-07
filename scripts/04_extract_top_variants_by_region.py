import pandas as pd

INPUT_FILE = "data/raw/mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"
OUTPUT_FILE = "results/water_permeability/top_variants_by_region.tsv"

regions = [
    ("Chr07_5_7Mb", "Chr07", 5_000_000, 7_000_000),
    ("Chr07_32_33_5Mb", "Chr07", 32_000_000, 33_500_000),
    ("Chr09_27_30Mb", "Chr09", 27_000_000, 30_000_000),
    ("Chr05_0_5_1_5Mb", "Chr05", 500_000, 1_500_000),
    ("Chr11_7_17Mb", "Chr11", 7_000_000, 17_000_000),
]

df = pd.read_csv(INPUT_FILE, sep="\t")

rows = []

for region_name, chrom, start, end in regions:
    sub = df[
        (df["chr"] == chrom)
        & (df["pos"] >= start)
        & (df["pos"] <= end)
        & (df["pval"] > 0)
    ].copy()

    if sub.empty:
        raise ValueError(f"No variants found in {chrom}:{start}-{end}")

    top = sub.loc[sub["pval"].idxmin()]

    rows.append({
        "region": region_name,
        "chr": top["chr"],
        "pos": int(top["pos"]),
        "allele1": top["allele1"],
        "allele0": top["allele0"],
        "af": top["af"],
        "beta": top["beta"],
        "pval": top["pval"],
    })

out = pd.DataFrame(rows)
out.to_csv(OUTPUT_FILE, sep="\t", index=False)

pd.set_option("display.float_format", "{:.3e}".format)
print(out)
print(f"Saved: {OUTPUT_FILE}")
