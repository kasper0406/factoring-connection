"""Check the residue-table identities used by the modular summary bound."""

import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

from sympy import primerange

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.dense_count_prediction import quarter_sign


def cell_counts(pairs, primary, shifts, previous):
    """Group product lifts by their public residue and joint count label."""
    base = primary * previous
    cells = defaultdict(list)
    for left, right in pairs:
        modulus = left * right
        labels = tuple(quarter_sign((left + shift) * (right + shift), primary) for shift in shifts)
        residue = modulus % base
        cells[residue, labels].append((modulus - residue) // base)
    return cells


def check_pool(cells, total, moduli):
    baseline = defaultdict(list)
    for (residue, _), lifts in cells.items():
        baseline[residue].append(len(lifts))
    without_extra = sum(max(counts) for counts in baseline.values())
    lower, upper = min(min(row) for row in cells.values()), max(max(row) for row in cells.values())
    powers = 0
    best_gain = Fraction(0)
    for extra in moduli:
        distance = 0
        energy = 0
        observation = defaultdict(list)
        for (residue, _), lifts in cells.items():
            counts = Counter(lift % extra for lift in lifts)
            size = len(lifts)
            distance += (extra - len(counts)) * size
            distance += sum(abs(extra * value - size) for value in counts.values())
            energy += extra * sum(value**2 for value in counts.values()) - size**2
            for low, count in counts.items():
                observation[residue, low].append(count)
        tv = Fraction(distance, 2 * extra * total)
        gain = Fraction(sum(max(counts) for counts in observation.values()) - without_extra, total)
        assert 0 <= gain <= tv <= 1
        assert 4 * total**2 * tv**2 <= len(cells) * energy
        best_gain = max(best_gain, gain)
        powers += energy
    # Nonzero additive frequencies for each prime are all reduced fractions.
    assert powers <= total * (upper - lower + 1 + max(moduli)**2)
    return best_gain


def check():
    cases = ((1024, 67, (1,), (1, 11)), (4096, 131, (1, 9, 25), (1, 11, 11 * 13)))
    for size, primary, shifts, previous_moduli in cases:
        moduli = list(map(int, primerange(17, 47)))
        primes = list(map(int, primerange(size, 5 * size // 4 + 1)))
        pairs = list(combinations(primes, 2))
        for previous in previous_moduli:
            cells = cell_counts(pairs, primary, shifts, previous)
            assert sum(map(len, cells.values())) == len(pairs)
            gain = check_pool(cells, len(pairs), moduli)
            combined = primary * previous * max(moduli)
            print(f'{len(pairs)} products, {len(shifts)} signs, CRT modulus up to {combined}, '
                  f'largest finite-sample optimal success gain {float(gain):.4f}')


if __name__ == '__main__':
    check()
