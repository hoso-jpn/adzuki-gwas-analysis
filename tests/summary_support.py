"""Synthetic customer summary metadata with noncanonical source column names."""

from pathlib import Path

SUMMARY_HEADER = "CHROM\tBP\tREF\tALT\tEA\tOA\tB\tSE\tP\tN\tEAF\n"


def write_summary_metadata(path: Path, **overrides) -> Path:
    values = dict(
        schema_version=1,
        dataset_id="synthetic_trait",
        cohort_id="synthetic-cohort",
        analysis_id="analysis-1",
        reference="Synthetic",
        assembly_id="synthetic-v1",
        species="synthetic plant",
        trait="trait",
        trait_unit="g",
        trait_coding="higher mass",
        test="wald",
        effect_type="beta",
        effect_scale="g",
        standard_error_scale="g",
        effect_strand="forward",
        target_effect_allele="ALT",
        sample_size_definition="analyzed_individuals",
    )
    values.update(overrides)
    import json

    text = "\n".join(f"{key} = {json.dumps(value)}" for key, value in values.items())
    text += (
        '\n[columns]\nchr = "CHROM"\npos = "BP"\nref = "REF"\nalt = "ALT"\n'
        'effect_allele = "EA"\nother_allele = "OA"\neffect = "B"\nstandard_error = "SE"\n'
        'pvalue = "P"\nsample_size = "N"\neffect_allele_frequency = "EAF"\n'
    )
    path.write_text(text)
    return path
