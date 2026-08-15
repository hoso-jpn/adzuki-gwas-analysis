"""The column-level input contract for a GWAS summary-statistics file.

``REQUIRED_COLUMNS`` is the exact, ordered header observed in all 6 real
Dryad files (verified directly against the downloaded, extracted archive --
not assumed from documentation). ``rs`` is deliberately excluded from the
variant-identity key: it is legitimately absent or ``"."`` for this species
(no dbSNP-equivalent registry), so it must never be relied on as a unique
key.
"""

from __future__ import annotations

import math
from collections import Counter

from adzuki_gwas_analysis.errors import RowValidationError, SchemaError

REQUIRED_COLUMNS: tuple[str, ...] = (
    "chr",
    "rs",
    "pos",
    "n_miss",
    "allele1",
    "allele0",
    "af",
    "beta",
    "se",
    "logl_H1",
    "l_remle",
    "l_mle",
    "p_wald",
    "pval",
    "p_score",
)

# Columns whose value must be a finite float strictly greater than 0 and
# less than or equal to 1 (a p-value contract, not merely "any finite
# float"). All three GEMMA test statistics share this contract.
_PVALUE_COLUMNS: tuple[str, ...] = ("p_wald", "pval", "p_score")

# Columns whose value must merely be a finite float (no range restriction
# beyond being a real, non-NaN, non-infinite number).
_FINITE_FLOAT_COLUMNS: tuple[str, ...] = ("beta", "se", "logl_H1", "l_remle", "l_mle")

VariantKey = tuple[str, str, str, str]


def validate_header(*, dataset_id: str, path: str, header: list[str]) -> None:
    """Validate that ``header`` matches :data:`REQUIRED_COLUMNS` exactly.

    Checks both that no column name repeats and that the observed header is
    exactly (not just a superset or same-set-different-order of)
    ``REQUIRED_COLUMNS``, since a silently reordered or renamed column would
    otherwise be misread downstream without any error.
    """
    counts = Counter(header)
    duplicates = sorted(name for name, count in counts.items() if count > 1)
    if duplicates:
        raise SchemaError(
            dataset_id=dataset_id,
            path=path,
            reason=f"duplicate column name(s) in header: {duplicates}",
        )

    if tuple(header) != REQUIRED_COLUMNS:
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        unexpected = [c for c in header if c not in REQUIRED_COLUMNS]
        raise SchemaError(
            dataset_id=dataset_id,
            path=path,
            reason=(
                f"header does not match the required contract exactly. "
                f"missing={missing}, unexpected={unexpected}, "
                f"observed_order={list(header)}, required_order={list(REQUIRED_COLUMNS)}"
            ),
        )


def _parse_finite_float(
    *, dataset_id: str, path: str, row_number: int, column: str, raw_value: str
) -> float:
    try:
        value = float(raw_value)
    except ValueError as exc:
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column=column,
            value=raw_value,
            reason="not parseable as a float",
        ) from exc
    if not math.isfinite(value):
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column=column,
            value=raw_value,
            reason="must be finite (not NaN or infinite)",
        )
    return value


def _parse_int(*, dataset_id: str, path: str, row_number: int, column: str, raw_value: str) -> int:
    # int() itself would silently accept "3.0"-style strings via float
    # round-tripping in some codebases; require a strict base-10 integer
    # literal so "3.5" or "3e2" are correctly rejected, not truncated.
    stripped = raw_value.strip()
    if not (stripped.lstrip("-").isdigit()):
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column=column,
            value=raw_value,
            reason="not a base-10 integer literal",
        )
    return int(stripped)


def validate_row(*, dataset_id: str, path: str, row_number: int, row: dict[str, str]) -> VariantKey:
    """Validate one data row against the full column contract.

    ``row_number`` is 1-indexed and counts data rows only (matching
    :class:`~adzuki_gwas_analysis.errors.RowValidationError`'s convention).
    Returns the ``(chr, pos, allele0, allele1)`` variant-identity key so the
    caller can track duplicates across the whole file without re-parsing.

    Raises :class:`~adzuki_gwas_analysis.errors.RowValidationError` on the
    first violation found; does not clamp, drop, or silently coerce any
    value.
    """
    chrom = row["chr"].strip()
    if not chrom:
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column="chr",
            value=row["chr"],
            reason="must not be empty",
        )

    allele1 = row["allele1"].strip()
    if not allele1:
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column="allele1",
            value=row["allele1"],
            reason="must not be empty",
        )

    allele0 = row["allele0"].strip()
    if not allele0:
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column="allele0",
            value=row["allele0"],
            reason="must not be empty",
        )

    pos = _parse_int(
        dataset_id=dataset_id, path=path, row_number=row_number, column="pos", raw_value=row["pos"]
    )
    if pos <= 0:
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column="pos",
            value=row["pos"],
            reason="must be a positive integer",
        )

    n_miss = _parse_int(
        dataset_id=dataset_id,
        path=path,
        row_number=row_number,
        column="n_miss",
        raw_value=row["n_miss"],
    )
    if n_miss < 0:
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column="n_miss",
            value=row["n_miss"],
            reason="must be a non-negative integer",
        )

    af = _parse_finite_float(
        dataset_id=dataset_id, path=path, row_number=row_number, column="af", raw_value=row["af"]
    )
    if not (0.0 <= af <= 1.0):
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column="af",
            value=row["af"],
            reason="must satisfy 0 <= af <= 1",
        )

    for column in _PVALUE_COLUMNS:
        p_value = _parse_finite_float(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column=column,
            raw_value=row[column],
        )
        if not (0.0 < p_value <= 1.0):
            raise RowValidationError(
                dataset_id=dataset_id,
                path=path,
                row_number=row_number,
                column=column,
                value=row[column],
                reason="must satisfy 0 < p <= 1",
            )

    for column in _FINITE_FLOAT_COLUMNS:
        _parse_finite_float(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column=column,
            raw_value=row[column],
        )

    se = float(row["se"])
    if se < 0.0:
        raise RowValidationError(
            dataset_id=dataset_id,
            path=path,
            row_number=row_number,
            column="se",
            value=row["se"],
            reason="must satisfy se >= 0",
        )

    return (chrom, row["pos"].strip(), allele0, allele1)
