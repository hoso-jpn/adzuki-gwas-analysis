"""Generate small synthetic assets without copying any real reference sequence."""

from pathlib import Path

from adzuki_gwas_analysis.loader import compute_sha256


def write_bundle(
    root: Path,
    *,
    sequence: str = "ACGT" * 80,
    assembly: str = "synthetic-v1",
    dataset: str = "synthetic_trait",
    reference: str = "Synthetic",
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference.fasta"
    lines = [sequence[i : i + 40] for i in range(0, len(sequence), 40)]
    fasta.write_text(">chrA\n" + "\n".join(lines) + "\n", encoding="ascii")
    fai = root / "reference.fasta.fai"
    width = min(40, len(sequence))
    fai.write_text(f"chrA\t{len(sequence)}\t6\t{width}\t{width + 1}\n", encoding="ascii")
    annotation = root / "genes.gff3"
    annotation.write_text(
        f"##gff-version 3\nchrA\tsynthetic\tgene\t1\t{min(20, len(sequence))}"
        "\t.\t+\t.\tID=gene1;Name=SyntheticGene\n",
        encoding="ascii",
    )
    text = (
        f'schema_version = 1\nreference = "{reference}"\nassembly_id = "{assembly}"\n'
        'species = "synthetic plant"\nsource = "synthetic fixture"\nlicense = "CC0"\n'
        f'data_scope = "synthetic"\ndataset_ids = ["{dataset}"]\n'
    )
    for role, path in (("fasta", fasta), ("fai", fai), ("annotation", annotation)):
        text += (
            f'\n[{role}]\npath = "{path.name}"\nassembly_id = "{assembly}"\n'
            f'sha256 = "{compute_sha256(path)}"\n'
        )
    bundle = root / "reference_bundle.toml"
    bundle.write_text(text, encoding="utf-8")
    return bundle
