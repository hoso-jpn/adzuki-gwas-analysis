import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Create a regional association plot from GWAS summary statistics."
    )
    parser.add_argument("--input", required=True, help="GWAS association file")
    parser.add_argument("--chrom", required=True, help="Chromosome name, e.g. Chr07")
    parser.add_argument("--start", type=int, required=True, help="Start position")
    parser.add_argument("--end", type=int, required=True, help="End position")
    parser.add_argument("--output", required=True, help="Output PNG file")
    parser.add_argument("--title", default=None, help="Plot title")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    sub = df[
        (df["chr"] == args.chrom)
        & (df["pos"] >= args.start)
        & (df["pos"] <= args.end)
    ].copy()

    if sub.empty:
        raise ValueError(
            f"No variants found in {args.chrom}:{args.start}-{args.end}"
        )

    sub = sub[(sub["pval"] > 0) & (sub["pval"] <= 1)]
    sub["minuslog10p"] = -np.log10(sub["pval"])

    top = sub.loc[sub["pval"].idxmin()]

    plt.figure(figsize=(8, 4))
    plt.scatter(
        sub["pos"],
        sub["minuslog10p"],
        s=8,
        alpha=0.7,
    )

    plt.axhline(-np.log10(1e-5), linestyle="--")
    plt.scatter(
        [top["pos"]],
        [-np.log10(top["pval"])],
        s=40,
        marker="*",
    )

    plt.xlabel(f"Position on {args.chrom} (bp)")
    plt.ylabel("-log10(p)")
    plt.title(
        args.title
        or f"Regional Plot: {args.chrom}:{args.start}-{args.end}"
    )

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)

    print(f"Saved: {args.output}")
    print("Top variant in region:")
    print(top[["chr", "pos", "allele1", "allele0", "af", "beta", "pval"]])


if __name__ == "__main__":
    main()
