"""Finite checks for the exact block counts and leaf masses in the AND fiber."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from math import isclose, log2


@cache
def block_arities(width, slope):
    modulus = 1 << width
    return Counter(
        ((s - slope * r) % modulus, r & s)
        for r in range(modulus)
        for s in range(modulus)
    )


def arity(width, slope, r, s):
    modulus = 1 << width
    return block_arities(width, slope % modulus)[(s - slope * r) % modulus, r & s]


def check_slope_sum():
    tested = 0
    for width in range(1, 7):
        modulus = 1 << width
        for r in range(modulus):
            for s in range(modulus):
                total = sum(arity(width, h, r, s) for h in range(1, modulus, 2))
                # Multiply equation (1) by the number of odd slopes.
                predicted = modulus // 2
                for d in range(width):
                    if (r ^ s) & (1 << d):
                        zeros_above = width - d - 1 - ((r & s) >> (d + 1)).bit_count()
                        predicted += (1 << d) * 3**zeros_above
                assert total == predicted, (width, r, s, total, predicted)
                tested += 1
    return tested


def check_product_children():
    tested = 0
    for width in range(1, 5):
        modulus = 1 << width
        for j in (width, width + 2):
            low_modulus = 1 << j
            next_modulus = low_modulus * modulus
            low_pairs = {(1, 1), (1, low_modulus - 1), (low_modulus - 1, 1)}
            low_pairs.add((low_modulus // 2 | 1, low_modulus - 1))
            for x, y in low_pairs:
                slope = -y * pow(x, -1, modulus) % modulus
                for r in range(modulus):
                    for s in range(modulus):
                        p, q = x + low_modulus * r, y + low_modulus * s
                        n, shared = p * q % next_modulus, p & q
                        children = 0
                        for alternative in range(x, next_modulus, low_modulus):
                            partner = n * pow(alternative, -1, next_modulus) % next_modulus
                            assert partner % low_modulus == y
                            if alternative & partner == shared:
                                children += 1
                        assert children == arity(width, slope, r, s)
                        tested += 1
    return tested


def check_tree_mass():
    observations = pairs = 0
    for k, width in ((5, 2), (7, 2), (8, 3), (9, 3)):
        modulus = 1 << k
        block_modulus = 1 << width
        first_arities = Counter(
            (p * q % block_modulus, p & q)
            for p in range(1, block_modulus, 2)
            for q in range(1, block_modulus, 2)
        )
        masses = defaultdict(Fraction)
        for p in range(1, modulus, 2):
            for q in range(1, modulus, 2):
                n, shared = p * q % modulus, p & q
                denominator = first_arities[n % block_modulus, shared % block_modulus]
                slope = -q * pow(p, -1, block_modulus) % block_modulus
                for j in range(width, k, width):
                    actual_width = min(width, k - j)
                    mask = (1 << actual_width) - 1
                    r, s = (p >> j) & mask, (q >> j) & mask
                    denominator *= arity(actual_width, slope, r, s)
                masses[n, shared] += Fraction(1, denominator)
                pairs += 1
        assert all(0 < mass <= 1 for mass in masses.values())
        observations += len(masses)
    return observations, pairs


if __name__ == "__main__":
    gamma = 0.75 * log2(3) - 1
    binary_entropy = -0.25 * log2(0.25) - 0.75 * log2(0.75)
    assert isclose(gamma, 1 - binary_entropy)
    slope_words = check_slope_sum()
    product_blocks = check_product_children()
    fibers, pairs = check_tree_mass()
    print(f"Checked slope identity for {slope_words} word pairs and {product_blocks} product blocks.")
    print(f"Verified exact leaf mass inequality for {fibers} fibers of {pairs} ordered pairs.")
    print(f"Typical ordered-fiber exponent {gamma:.11f}")
