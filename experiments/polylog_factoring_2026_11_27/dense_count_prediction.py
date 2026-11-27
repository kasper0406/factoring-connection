"""Recover a component from noisy signs of shifted count residues."""

import sys
from collections import Counter
from fractions import Fraction
from math import cos, pi, sin
from pathlib import Path
from random import Random

from mpmath import mp
from sympy import sqrt_mod

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.large_mod_residue_evaluator import (
    find_component_in_residue,
)
from experiments.polylog_factoring_2026_11_27.quarter_fourier_recovery import (
    estimator,
    heavy_frequencies,
)


def quarter_sign(value, prime):
    residue = value % prime
    return 1 if 4 * min(residue, prime - residue) < prime else -1


def sum_candidates(modulus, prime, predict, rng, harmonics, threshold, samples, precise=False):
    """Return candidate residues for the component sum."""
    public_residue = modulus % prime
    trig_functions = (mp.cospi, mp.sinpi) if precise else (cos, sin)
    tested = set()
    for harmonic in harmonics:
        inverse = pow(harmonic, -1, prime)
        for trig in trig_functions:
            def demodulate(shift):
                phase = harmonic * (shift**2 + public_residue) % prime
                angle = mp.mpf(2 * phase) / prime if precise else 2 * pi * phase / prime
                return predict(shift) * trig(angle)

            energy = estimator(prime, demodulate, rng, samples, precise)
            for frequency in heavy_frequencies(prime, threshold, energy, precise):
                candidate = frequency * inverse % prime
                if candidate not in tested:
                    tested.add(candidate)
                    yield candidate


def component_from_sums(modulus, prime, candidates):
    """Translate sum candidates into residue classes and verify the components."""
    inverse_two = pow(2, -1, prime)
    searched = set()
    for total in candidates:
        discriminant = (total**2 - 4 * modulus) % prime
        for root in sqrt_mod(discriminant, prime, all_roots=True):
            residue = (total + root) * inverse_two % prime
            if residue in searched:
                continue
            searched.add(residue)
            component = find_component_in_residue(modulus, prime, residue, Counter())
            if component is not None and 1 < component < modulus and modulus % component == 0:
                return component
    return None


def check():
    rng = Random(261127)
    fixtures = ((101, 109, 17), (101**2, 109, 37), (1000003, 1000033, 1009),
                (1000000007, 1000000009, 32771))
    for left, right, prime in fixtures:
        modulus = left * right
        for harmonic in (1, 5):
            def planted(shift):
                return quarter_sign(harmonic * (left + shift) * (right + shift), prime)

            # Measure the planted predictor's actual agreement separately.
            correlation = sum(planted(shift) * quarter_sign((left + shift) * (right + shift), prime)
                              for shift in range(prime)) / prime
            assert correlation > 0
            candidates = list(sum_candidates(modulus, prime, planted, rng, (harmonic,),
                                             Fraction(1, 5), 16_384, precise=prime == 17))
            assert (left + right) % prime in candidates
            assert component_from_sums(modulus, prime, candidates) in (left, right)
            print(f'{modulus.bit_length()} bits, harmonic {harmonic}, correlation {correlation:.4f}, split verified')

    for prime in (17, 37, 67, 131):
        inverses = [pow(unit, -1, prime) for unit in range(1, prime)]
        # Enumerate all u at fixed public product and nonzero shift.
        max_bias = max(
            abs(sum(quarter_sign((unit + shift) * (product * inverse + shift), prime)
                    for unit, inverse in enumerate(inverses, 1))) / (prime - 1)
            for product in range(1, prime) for shift in range(1, prime)
        )
        print(f'auxiliary prime {prime}, maximum residue-only conditional bias {max_bias:.4f}')


if __name__ == '__main__':
    check()
