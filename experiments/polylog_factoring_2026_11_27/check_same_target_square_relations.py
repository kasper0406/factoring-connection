"""Check whether combining short rational representations helps on one square target.

Each usable row u/v represents the same unit t modulo N. After removing
factor-base primes, let w be the positive cofactor of |uv|. A subset
removes all rough primes from its product up to a rational square exactly
when its w values multiply to an integer square. An odd subset leaves
an odd power of t, from which a relation for t follows by dividing out
a known square. The sign of uv is retained in the factor base.

Run directly or as a module from the repository root. The default cohort
uses the matched targets and representations from check_rational_height_shell.
"""

from itertools import combinations
from math import gcd, isqrt, prod
from pathlib import Path
from random import Random
import sys

import gmpy2

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.check_rational_height_shell import representations
from src.common.instances import balanced_semiprime
from src.smoothbuild.core import primes_upto


def square(value):
    return isqrt(value) ** 2 == value


def odd_square_subset(values):
    """Return an odd subset mask with square product, or zero if none exists."""
    pending = [(value, 1 << i) for i, value in enumerate(values)]
    coprime = {}
    while pending:
        value, incidence = pending.pop()
        if incidence == 0:
            continue
        if square(value):
            continue
        for atom, other_incidence in coprime.items():
            common = gcd(value, atom)
            if common != 1:
                break
        else:
            coprime[value] = incidence
            continue
        del coprime[atom]
        # Common factors cancel wherever both row bits are set.
        pending.extend(((value // common, incidence),
                        (atom // common, other_incidence),
                        (common, incidence ^ other_incidence)))

    # Distinct surviving atoms are coprime nonsquares. A product is a square
    # exactly when every atom occurs an even number of times.
    rows = [0] * len(values)
    for column, incidence in enumerate(coprime.values()):
        while incidence:
            bit = incidence & -incidence
            rows[bit.bit_length() - 1] |= 1 << column
            incidence ^= bit

    pivots = {}
    for index, row in enumerate(rows):
        subset = 1 << index
        while row:
            lead = row.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = (row, subset)
                break
            pivot, provenance = pivots[lead]
            row ^= pivot
            subset ^= provenance
        if row == 0 and subset.bit_count() % 2:
            return subset
    return 0


def check_solver():
    rng = Random(34718)
    cases = [[6, 10, 15], [2, 2], [8, 18], [27, 35, 105], [4, 7], [8, 3, 6]]
    cases.extend([rng.randrange(1, 260) for _ in range(rng.randrange(1, 10))]
                 for _ in range(10_000))
    for values in cases:
        mask = odd_square_subset(values)
        expected = any(square(prod(subset)) for size in range(1, len(values) + 1, 2)
                       for subset in combinations(values, size))
        assert bool(mask) == expected
        if mask:
            assert mask.bit_count() % 2
            assert square(prod(v for i, v in enumerate(values) if (mask >> i) & 1))
    return len(cases)


def rough_part(value, base_product):
    common = int(gmpy2.gcd(value, base_product))
    while common > 1:
        value //= common
        common = gcd(value, common)
    return value


def square_relation(n, target, pairs, rough):
    """Return (number of rows, unit multiplier, signed smooth part), or None.

    Pairs are unit fractions for target; rough maps each pair to its
    absolute product with all factor-base prime powers removed.
    """
    pairs = sorted(pairs)
    values = [rough[pair] for pair in pairs]
    mask = next((1 << i for i, value in enumerate(values) if square(value)), 0)
    if not mask:
        mask = odd_square_subset(values)
    if not mask:
        return None

    selected = [pair for i, pair in enumerate(pairs) if (mask >> i) & 1]
    rough_product = prod(rough[pair] for pair in selected)
    root = isqrt(rough_product)
    assert root * root == rough_product and len(selected) % 2 == 1
    smooth = prod(u * v for u, v in selected) // rough_product
    denominators = prod(v for _, v in selected)
    multiplier = (pow(target, len(selected) // 2, n)
                  * denominators * pow(root, -1, n)) % n
    assert (target * multiplier * multiplier - smooth) % n == 0
    return len(selected), multiplier, smooth


def evaluate_target(n, target, base_limit):
    """Find a witnessed relation at height N for a unit target, if one exists."""
    target %= n
    if gcd(target, n) != 1:
        raise ValueError("target must be a unit modulo n")
    base_product = prod(primes_upto(base_limit))
    pairs, _ = representations(n, target)
    rough = {pair: rough_part(abs(pair[0] * pair[1]), base_product)
             for pair in pairs if gcd(pair[1], n) == 1}
    return square_relation(n, target, rough.keys(), rough)


def check_combined_relation():
    n, target = 247, 23
    pairs, _ = representations(n, target)
    rough = {pair: rough_part(abs(pair[0] * pair[1]), 6)
             for pair in pairs if gcd(pair[1], n) == 1}
    assert all(not square(value) for value in rough.values())
    length, _, smooth = evaluate_target(n, target, 3)
    assert length >= 3 and rough_part(abs(smooth), 6) == 1


def benchmark(bits):
    base_product = gmpy2.mpz(prod(primes_upto(bits * bits)))
    convergent_single = convergent_any = shell_single = shell_any = 0
    for seed in range(4):
        n = balanced_semiprime(bits, seed).N
        rng = Random(f'relation-baseline:{bits}:{seed}')
        for _ in range(1024):
            root = rng.randrange(1, n)
            while gcd(root, n) != 1:
                root = rng.randrange(1, n)
            target = root * root % n
            candidates, convergents = representations(n, target)
            rough = {pair: rough_part(abs(pair[0] * pair[1]), base_product)
                     for pair in candidates if gcd(pair[1], n) == 1}
            rel = square_relation(n, target, convergents.intersection(rough), rough)
            convergent_single += rel is not None and rel[0] == 1
            convergent_any += rel is not None
            rel = square_relation(n, target, rough.keys(), rough)
            shell_single += rel is not None and rel[0] == 1
            shell_any += rel is not None
    return convergent_single, convergent_any, shell_single, shell_any


if __name__ == "__main__":
    print("Small exhaustive subset comparisons:", check_solver())
    check_combined_relation()
    print("bits convergent_single convergent_any_odd shell_single shell_any_odd")
    for bits in (32, 48, 64, 80, 96):
        print(bits, *benchmark(bits), flush=True)
