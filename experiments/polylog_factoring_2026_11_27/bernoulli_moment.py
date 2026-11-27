"""Benchmark the degree-eight moment on lower-input windows modulo N squared.

The moment itself uses exact integer arithmetic. The dense Fourier transform
only diagnoses the resulting bit signal on these small, synthetic fixtures.
"""

import argparse
import json
from cmath import exp
from math import cos, gcd, pi
from random import Random
from sys import float_info
from time import perf_counter
from typing import NamedTuple


class Block(NamedTuple):
    upper_digit: int
    velocity: int
    divisor: int


def weight(modulus, digit):
    """Return modulus**8 times B_8(digit / modulus) / B_8(0)."""
    modulus_squared = modulus * modulus
    product = digit * (modulus - digit)
    product_squared = product * product
    return modulus_squared**4 - product_squared * (
        20 * modulus_squared**2 + 40 * modulus_squared * product + 30 * product_squared
    )


def window_units(n, c0, offset, width):
    """Visit the cyclic lower-input window, omitting public stationary blocks."""
    for index in range(width):
        x0 = (index - offset) % n
        if x0 not in (c0, n - c0) and gcd(x0, n) == 1:
            yield x0


def progression_blocks(n, c0, offset, width):
    square = n * n
    blocks = []
    for x0 in window_units(n, c0, offset, width):
        inverse = pow(x0, -1, square)
        blocks.append(Block(
            ((x0 + c0 * c0 * inverse) % square) // n,
            (2 * c0 * inverse) % n,
            gcd(x0 * x0 - c0 * c0, n),
        ))
    return blocks


def progression_moment(blocks, j):
    """Numerator of the corrected moment with denominator N**8."""
    # The average weight of the N lifts of a block equals W_d(s), where
    # the upper digits are s, s+d, ..., s+N-d with multiplicity d.
    return sum(
        weight(d, (upper + velocity * j) % d) for upper, velocity, d in blocks
    )


def direct_moment(n, c0, j, offset, width):
    """Recompute the same numerator by explicit inversion for every lift."""
    square = n * n
    c_squared = (c0 + n * j)**2
    total = 0
    for x0 in window_units(n, c0, offset, width):
        for x in range(x0, square, n):
            upper = ((x + c_squared * pow(x, -1, square)) % square) // n
            total += weight(n, upper)
    numerator, remainder = divmod(total, n)
    assert remainder == 0
    return numerator


def gap_bit(numerator, scale, j):
    """Sign above a 1/64 gap; exercise one legal answer within the gap."""
    if 64 * abs(numerator) > scale:
        return 1 if numerator > 0 else -1
    return 1 if j % 2 == 0 else -1


def unit_peak(bits):
    """Find the strongest unit frequency by dense DFT, up to conjugation."""
    n = len(bits)
    phases = [exp(-2j * pi * digit / n) for digit in range(n)]
    frequencies = [k for k in range(1, n // 2 + 1) if gcd(k, n) == 1]
    peak, frequency = max(
        (abs(sum(bit * phases[(k * j) % n] for j, bit in enumerate(bits))) / n, k)
        for k in frequencies
    )
    return peak, frequency, len(frequencies)


def fixture_window(n, rng):
    # Fixture construction uses the hidden roots to select a singleton. The
    # moment and Fourier routines take only the public window and its answers.
    mixed_roots = [u for u in range(2, n - 1) if u * u % n == 1]
    width = max(1, n // (4 * (len(mixed_roots) + 2)))
    for attempt in range(1, 10_001):
        c0, offset = rng.randrange(1, n), rng.randrange(n)
        if gcd(c0, n) != 1:
            continue
        selected = [u for u in mixed_roots if (c0 * u + offset) % n < width]
        if len(selected) == 1:
            return c0, offset, width, selected[0], attempt
    raise AssertionError("No singleton found in fixture")


def benchmark(n, rng):
    c0, offset, width, hidden, fixture_attempts = fixture_window(n, rng)
    started = perf_counter()
    blocks = progression_blocks(n, c0, offset, width)
    preparation_seconds = perf_counter() - started
    assert sum(block.divisor == n for block in blocks) == 1

    started = perf_counter()
    numerators = [progression_moment(blocks, j) for j in range(n)]
    progression_seconds = perf_counter() - started
    direct_digits = sorted({0, 1, n // 2, n - 1})
    started = perf_counter()
    for j in direct_digits:
        assert numerators[j] == direct_moment(n, c0, j, offset, width)
    direct_seconds = perf_counter() - started

    scale = n**8
    bits = [gap_bit(value, scale, j) for j, value in enumerate(numerators)]
    started = perf_counter()
    peak, frequency, dft_frequencies = unit_peak(bits)
    dft_seconds = perf_counter() - started
    expected_frequency = 2 * hidden % n
    assert frequency in (expected_frequency, n - expected_frequency)
    assert peak > 0.5
    divisor = gcd(frequency - 2, n)
    assert 1 < divisor < n and n % divisor == 0

    stationary, = [block for block in blocks if block.divisor == n]
    residues = [(stationary.upper_digit + stationary.velocity * j) % n for j in range(n)]
    remainder = [value - weight(n, z) for value, z in zip(numerators, residues)]
    assert 64 * max(map(abs, remainder)) < scale
    clean_bits = [1 if cos(2 * pi * z / n) > 0 else -1 for z in residues]
    return {
        "n": n,
        "c0": c0,
        "offset": offset,
        "width": width,
        "fixture_attempts": fixture_attempts,
        "unit_blocks": len(blocks),
        "block_divisors": sorted({block.divisor for block in blocks}),
        "direct_digits": len(direct_digits),
        "direct_lift_weights": len(direct_digits) * len(blocks) * n,
        "progression_digits": n,
        "progression_block_weights": n * len(blocks),
        "max_numerator_bits": max(abs(value).bit_length() for value in numerators),
        "denominator_bits": scale.bit_length(),
        "maximum_normalized_background": max(map(abs, remainder)) / scale,
        "gap_digits": sum(64 * abs(value) <= scale for value in numerators),
        "bit_disagreements": sum(a != b for a, b in zip(bits, clean_bits)),
        "dft_frequencies": dft_frequencies,
        "dft_float_significand_bits": float_info.mant_dig,
        "peak": round(peak, 8),
        "frequency": frequency,
        "divisor": divisor,
        "seconds": {
            "setup": round(preparation_seconds, 6),
            "progression": round(progression_seconds, 6),
            "direct": round(direct_seconds, 6),
            "dft": round(dft_seconds, 6),
        },
    }


def check_progressions():
    # Verify each progression, including all residues for partial prime-power
    # divisors. Longer windows also check sums and public-block removal.
    for n in (15, 45, 225, 315, 1155):
        for d in range(1, n + 1):
            if n % d == 0:
                for residue in range(d):
                    assert d * sum(weight(n, z) for z in range(residue, n, d)) == (
                        n * weight(d, residue)
                    )
    for n, c0, j, offset, width in ((45, 2, 24, 17, 45), (225, 7, 112, 3, 67)):
        blocks = progression_blocks(n, c0, offset, width)
        assert progression_moment(blocks, j) == direct_moment(n, c0, j, offset, width)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--seed", type=int, default=261127)
    args = cli.parse_args()
    rng = Random(args.seed)
    check_progressions()
    for modulus in (15, 45, 225, 315, 1155, 1575):
        print(json.dumps(benchmark(modulus, rng), sort_keys=True))
