"""Small, explicit DNA screening models; no PCR or specificity guarantees."""

from __future__ import annotations

import math
from itertools import groupby

from adzuki_gwas_analysis.reference import reverse_complement

# SantaLucia (1998), DNA/DNA nearest-neighbor parameters (kcal/mol, cal/K/mol).
# https://doi.org/10.1073/pnas.95.4.1460 ; the DNA_NN3 convention documented by Biopython.
_NN = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}


def tm_nn(sequence: str, *, sodium_mM: float = 50, strand_concentration_nM: float = 25) -> float:
    """Perfect-match DNA/DNA Tm with entropy salt correction; Mg/DMSO not modeled."""
    if len(sequence) < 2 or any(base not in "ACGT" for base in sequence):
        raise ValueError("Tm requires an uppercase unambiguous DNA oligo")
    if not all(math.isfinite(x) and x > 0 for x in (sodium_mM, strand_concentration_nM)):
        raise ValueError("salt and strand concentrations must be finite and positive")
    enthalpy = entropy = 0.0
    for base in (sequence[0], sequence[-1]):
        dh, ds = (2.3, 4.1) if base in "AT" else (0.1, -2.8)
        enthalpy += dh
        entropy += ds
    for index in range(len(sequence) - 1):
        dh, ds = _NN[sequence[index : index + 2]]
        enthalpy += dh
        entropy += ds
    self_complementary = sequence == reverse_complement(sequence)
    if self_complementary:
        entropy -= 1.4
    concentration = strand_concentration_nM * 1e-9 / (1 if self_complementary else 2)
    entropy += 0.368 * (len(sequence) - 1) * math.log(sodium_mM * 1e-3)
    return 1000 * enthalpy / (entropy + 1.987 * math.log(concentration)) - 273.15


def complementary_run(first: str, second: str) -> int:
    """Longest contiguous complementary tract; a screen, not hairpin/dimer thermodynamics."""
    reverse = reverse_complement(second)
    previous = [0] * (len(reverse) + 1)
    longest = 0
    for base in first:
        current = [0]
        for index, other in enumerate(reverse):
            value = previous[index] + 1 if base == other else 0
            current.append(value)
            longest = max(longest, value)
        previous = current
    return longest


def primer_metrics(
    sequence: str, *, sodium_mM: float, strand_concentration_nM: float
) -> dict[str, float]:
    return {
        "length": len(sequence),
        "gc_fraction": (sequence.count("G") + sequence.count("C")) / len(sequence),
        "tm_C": tm_nn(
            sequence, sodium_mM=sodium_mM, strand_concentration_nM=strand_concentration_nM
        ),
        "max_homopolymer": max(len(list(bases)) for _, bases in groupby(sequence)),
        "self_complementary_run": complementary_run(sequence, sequence),
    }
