"""Shared test helpers for the analysis package's tests.

Builds small, synthetic, schema-v1-valid manifests and region configs
pointing at tiny fixture files under ``tests/fixtures/`` -- no test in this
package's suite reads anything under ``data/raw/`` (the real, un-tracked
Dryad files).
"""

from __future__ import annotations

from pathlib import Path

from adzuki_gwas_analysis.loader import compute_sha256, count_data_rows
from adzuki_gwas_analysis.manifest import SCHEMA_V1_DATASETS

HEADER = (
    "chr\trs\tpos\tn_miss\tallele1\tallele0\taf\tbeta\tse\t"
    "logl_H1\tl_remle\tl_mle\tp_wald\tpval\tp_score"
)

#: Canonical member_filename for the one dataset these tests actually read.
MIYAGI_WATER_PERMEABILITY_FILENAME = "mapped_to_Miyagi_water_permeability.maf_0.05.assoc.txt"

#: The other 5 schema-v1 datasets, never read by these tests -- only present
#: so the manifest satisfies "exactly 6 canonical datasets".
_OTHER_CANONICAL_ENTRIES = """
[[datasets]]
dataset_id = "miyagi_red_seedcoat"
reference = "Miyagi"
trait = "red_seedcoat"
member_filename = "mapped_to_Miyagi_red_seedcoat.maf_0.05.assoc.txt"
member_sha256 = "{fake}"
row_count = 1

[[datasets]]
dataset_id = "miyagi_mottled_black_seedcoat"
reference = "Miyagi"
trait = "mottled_black_seedcoat"
member_filename = "mapped_to_Miyagi_mottled_black_seedcoat.maf_0.05.assoc.txt"
member_sha256 = "{fake}"
row_count = 1

[[datasets]]
dataset_id = "shumari_water_permeability"
reference = "Shumari"
trait = "water_permeability"
member_filename = "mapped_to_Shumari_water_permeability.maf_0.05.assoc.txt"
member_sha256 = "{fake}"
row_count = 1

[[datasets]]
dataset_id = "shumari_red_seedcoat"
reference = "Shumari"
trait = "red_seedcoat"
member_filename = "mapped_to_Shumari_red_seedcoat.maf_0.05.assoc.txt"
member_sha256 = "{fake}"
row_count = 1

[[datasets]]
dataset_id = "shumari_mottled_black_seedcoat"
reference = "Shumari"
trait = "mottled_black_seedcoat"
member_filename = "mapped_to_Shumari_mottled_black_seedcoat.maf_0.05.assoc.txt"
member_sha256 = "{fake}"
row_count = 1
""".format(fake="b" * 64)

_MANIFEST_TEMPLATE = (
    """
schema_version = 1

[dryad]
doi = "10.5061/dryad.8w9ghx3xv"
dataset_id = 149675
version_id = 356599
version_number = 6
publication_doi = "10.1126/science.ads2871"

[archive]
filename = "adzuki_GWAS_data.zip"
size_bytes = 1
sha256 = "{archive_sha}"

[pvalue_columns]
p_wald = "Wald test p-value"
pval = "Likelihood ratio test (LRT) p-value"
p_score = "Score test p-value"
primary = "pval"

[[datasets]]
dataset_id = "miyagi_water_permeability"
reference = "Miyagi"
trait = "water_permeability"
member_filename = "{filename}"
member_sha256 = "{sha256}"
row_count = {row_count}
"""
    + _OTHER_CANONICAL_ENTRIES
)


def write_dataset_file(data_dir: Path, rows: list[str]) -> Path:
    """Write ``rows`` (tab-separated data lines, no header) as the Miyagi
    water-permeability fixture file under ``data_dir``; returns its path."""
    path = data_dir / MIYAGI_WATER_PERMEABILITY_FILENAME
    path.write_text(HEADER + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_manifest(manifest_path: Path, dataset_path: Path) -> Path:
    """Write a schema-v1 manifest whose ``miyagi_water_permeability`` entry
    matches ``dataset_path``'s real checksum/row-count; the other 5 entries
    are structurally valid placeholders that are never read."""
    manifest_path.write_text(
        _MANIFEST_TEMPLATE.format(
            archive_sha="a" * 64,
            filename=dataset_path.name,
            sha256=compute_sha256(dataset_path),
            row_count=count_data_rows(dataset_path),
        ),
        encoding="utf-8",
    )
    return manifest_path


#: The 6 (reference, trait) pairs, in the exact order manifest.toml itself declares them
#: (Miyagi x 3 traits, then Shumari x 3 traits) -- used by the batch ("Issue #9") fixture
#: helpers below, which (unlike write_manifest/write_dataset_file above) write a *real*,
#: independently checksummed fixture file for all 6 canonical datasets, not just
#: miyagi_water_permeability.
CANONICAL_REFERENCE_TRAIT_ORDER: tuple[tuple[str, str], ...] = (
    ("Miyagi", "water_permeability"),
    ("Miyagi", "red_seedcoat"),
    ("Miyagi", "mottled_black_seedcoat"),
    ("Shumari", "water_permeability"),
    ("Shumari", "red_seedcoat"),
    ("Shumari", "mottled_black_seedcoat"),
)


def write_six_dataset_files(
    data_dir: Path, rows_by_dataset_id: dict[str, list[str]]
) -> dict[str, Path]:
    """Write one small, schema-v1-valid fixture file per canonical dataset_id.

    ``rows_by_dataset_id`` must have exactly the 6 canonical dataset_ids as keys (see
    :data:`CANONICAL_REFERENCE_TRAIT_ORDER` /
    :data:`adzuki_gwas_analysis.manifest.SCHEMA_V1_DATASETS`). Returns
    ``{dataset_id: path}``, one real file per dataset -- unlike ``write_dataset_file``
    above, none of these 6 are placeholders.
    """
    paths: dict[str, Path] = {}
    for reference, trait in CANONICAL_REFERENCE_TRAIT_ORDER:
        spec = SCHEMA_V1_DATASETS[(reference, trait)]
        rows = rows_by_dataset_id[spec.dataset_id]
        path = data_dir / spec.member_filename
        path.write_text(HEADER + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
        paths[spec.dataset_id] = path
    return paths


def write_full_manifest(manifest_path: Path, dataset_paths: dict[str, Path]) -> Path:
    """Write a schema-v1 manifest whose all 6 entries match real fixture files.

    ``dataset_paths`` must be the return value of :func:`write_six_dataset_files` (or an
    equivalent ``{dataset_id: path}`` mapping covering all 6 canonical dataset_ids); each
    entry's ``member_sha256``/``row_count`` is computed directly from that dataset's real
    file, in :data:`CANONICAL_REFERENCE_TRAIT_ORDER`'s order, matching
    ``manifest.toml``'s own declared dataset order (Miyagi x 3 traits, then Shumari x 3
    traits).
    """
    entries = []
    for reference, trait in CANONICAL_REFERENCE_TRAIT_ORDER:
        spec = SCHEMA_V1_DATASETS[(reference, trait)]
        path = dataset_paths[spec.dataset_id]
        entries.append(
            f"""
[[datasets]]
dataset_id = "{spec.dataset_id}"
reference = "{reference}"
trait = "{trait}"
member_filename = "{spec.member_filename}"
member_sha256 = "{compute_sha256(path)}"
row_count = {count_data_rows(path)}
"""
        )
    manifest_path.write_text(
        """
schema_version = 1

[dryad]
doi = "10.5061/dryad.8w9ghx3xv"
dataset_id = 149675
version_id = 356599
version_number = 6
publication_doi = "10.1126/science.ads2871"

[archive]
filename = "adzuki_GWAS_data.zip"
size_bytes = 1
sha256 = "{archive_sha}"

[pvalue_columns]
p_wald = "Wald test p-value"
pval = "Likelihood ratio test (LRT) p-value"
p_score = "Score test p-value"
primary = "pval"
""".format(archive_sha="a" * 64)
        + "".join(entries),
        encoding="utf-8",
    )
    return manifest_path


#: Row template used by build_candidate_enabled_batch_fixture below: chr, pos, and pval are
#: filled in per row; every other column is a fixed, schema-v1-valid placeholder value.
_CANDIDATE_FIXTURE_ROW_TEMPLATE = (
    "{chrom}\t.\t{pos}\t0\tA\tG\t0.30\t0.05\t0.01\t70.0\t30.0\t28.0\t1e-3\t{pval}\t1e-3"
)

#: (chrom, pos, pval) rows per dataset_id, deliberately varied so report tests exercise:
#: multiple references, multiple traits, a dataset with zero candidates
#: (miyagi_mottled_black_seedcoat), a single-signal dataset, and a multi-signal dataset.
_CANDIDATE_FIXTURE_ROWS_BY_DATASET_ID: dict[str, list[tuple[str, int, str]]] = {
    "miyagi_water_permeability": [
        ("Chr01", 1_000_000, "1e-9"),
        ("Chr01", 1_000_200, "2e-9"),
        ("Chr01", 5_000_000, "9e-1"),
    ],
    "miyagi_red_seedcoat": [
        ("Chr01", 1_000_000, "1e-9"),
        ("Chr02", 2_000_000, "2e-9"),
        ("Chr02", 9_000_000, "9e-1"),
    ],
    "miyagi_mottled_black_seedcoat": [
        ("Chr01", 1_000_000, "8e-1"),
        ("Chr01", 2_000_000, "9e-1"),
        ("Chr02", 3_000_000, "7e-1"),
    ],
    "shumari_water_permeability": [
        ("Chr01", 1_000_000, "1e-9"),
        ("Chr01", 4_000_000, "8e-1"),
        ("Chr02", 5_000_000, "9e-1"),
    ],
    "shumari_red_seedcoat": [
        ("Chr01", 1_000_000, "1e-9"),
        ("Chr01", 1_000_300, "2e-9"),
        ("Chr02", 6_000_000, "3e-9"),
    ],
    "shumari_mottled_black_seedcoat": [
        ("Chr01", 1_000_000, "1e-9"),
        ("Chr02", 2_000_000, "9e-1"),
        ("Chr02", 3_000_000, "8e-1"),
    ],
}


def build_candidate_enabled_batch_fixture(
    tmp_path: Path, *, clustering_distance: int = 1000
) -> Path:
    """Build a real, schema-v2 candidate-enabled ``batch`` output for report tests.

    Runs the actual :func:`adzuki_gwas_analysis.analysis.batch.run_batch` against small
    synthetic fixtures rather than hand-authoring a ``batch_summary.tsv`` and per-dataset
    TSVs by hand -- this guarantees the fixture always matches the real schema/consistency
    contract :mod:`adzuki_gwas_analysis.analysis.report_validation` checks, with no risk of
    drifting from it. ``miyagi_mottled_black_seedcoat`` deliberately has zero candidates
    (every p-value is large), so report tests exercise that state for real rather than by
    construction. Returns the published ``analysis_dir`` (a fresh subdirectory of
    ``tmp_path``).
    """
    from adzuki_gwas_analysis.analysis.batch import run_batch

    data_dir = tmp_path / "fixture_data"
    data_dir.mkdir()
    manifest_path = tmp_path / "fixture_manifest.toml"
    analysis_dir = tmp_path / "fixture_analysis_dir"

    rows_by_dataset_id = {
        dataset_id: [
            _CANDIDATE_FIXTURE_ROW_TEMPLATE.format(chrom=chrom, pos=pos, pval=pval)
            for chrom, pos, pval in rows
        ]
        for dataset_id, rows in _CANDIDATE_FIXTURE_ROWS_BY_DATASET_ID.items()
    }
    dataset_paths = write_six_dataset_files(data_dir, rows_by_dataset_id)
    write_full_manifest(manifest_path, dataset_paths)

    run_batch(
        manifest_path=manifest_path,
        data_dir=data_dir,
        output_dir=analysis_dir,
        clustering_distance=clustering_distance,
    )
    return analysis_dir


def write_region_config(config_path: Path, regions: list[dict[str, object]]) -> Path:
    """Write a region config TOML for ``miyagi_water_permeability`` with ``regions``."""
    lines = [
        "config_schema_version = 1",
        'dataset_id = "miyagi_water_permeability"',
        "",
    ]
    for region in regions:
        lines.append("[[regions]]")
        for key, value in region.items():
            if isinstance(value, str):
                lines.append(f'{key} = "{value}"')
            else:
                lines.append(f"{key} = {value}")
        lines.append("")
    config_path.write_text("\n".join(lines), encoding="utf-8")
    return config_path
