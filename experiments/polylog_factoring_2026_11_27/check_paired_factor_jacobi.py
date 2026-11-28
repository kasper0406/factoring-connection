"""Verify independent local signs, the paired sieve and a block-histogram collision."""

from itertools import combinations
from math import gcd, isqrt, lcm
from random import Random

from sympy import primerange

if __package__:
    from .paired_factor_jacobi import (composite_single_coordinate_witness, jacobi, paired_moduli,
                                       paired_scan, quotient_candidates)
else:
    from paired_factor_jacobi import (composite_single_coordinate_witness, jacobi, paired_moduli,
                                      paired_scan, quotient_candidates)


def check_periods():
    checked = shared = 0
    for p, q in combinations(primerange(5, 73), 2):
        if q >= 2*p:
            continue
        n = p*q
        assert (p, q) in quotient_candidates(n)
        for x, y in quotient_candidates(n):
            if x == p:
                continue
            period = lcm(x, y)
            assert gcd(period, n) == 1 and isqrt(x*y)**2 != x*y
            counts = {(s, t): 0 for s in (-1, 1) for t in (-1, 1)}
            for j in range(period):
                d = 1+4*j*n
                signs = jacobi(d, x), jacobi(d, y)
                if 0 not in signs:
                    counts[signs] += 1
            assert len(set(counts.values())) == 1 and counts[(1, 1)] > 0
            checked += 1
            shared += gcd(x, y) > 1
    return checked, shared


def block_example(n):
    moduli = [1+4*t*n for t in (3, 6, 18)]
    histograms = []
    for start in (42697, 39308):
        pairs = list(quotient_candidates(n, start, start+5))
        rows = [[(jacobi(d, x), jacobi(d, y)) for x, y in pairs] for d in moduli]
        histograms.append([tuple(sorted(row)) for row in rows])
        survivors = [pair for pair, signs in zip(pairs, zip(*rows)) if all(z == (1, 1) for z in signs)]
        assert survivors == ([(42697, 70019)] if start == 42697 else [])
    assert histograms[0] == histograms[1] == [((1, -1), (1, 1))]*3
    return histograms[0]


if __name__ == "__main__":
    print("False-pair periods and noncoprime coordinate pairs:", check_periods())
    n = 2989601243
    moduli = paired_moduli(n, 20, rng=Random(20261127))
    survivors, counts = paired_scan(n, moduli)
    assert survivors == [(42697, 70019)]
    one_side = list(quotient_candidates(n))
    single_counts = [len(one_side)]
    for d in moduli:
        one_side = [(x, y) for x, y in one_side if jacobi(d, x) == 1]
        single_counts.append(len(one_side))
    print("Seeded paired survivor counts:", counts)
    print("Matched single-coordinate counts:", single_counts)
    print("Both six-integer blocks have per-row histograms:", block_example(n))
    lower, upper = isqrt(2*n)//2+1, isqrt(n)
    witness = composite_single_coordinate_witness(lower, upper, moduli[:3])
    assert witness is not None
    x, u, v = witness
    assert lower <= x <= upper and x == u*v and 1 < u < v and isqrt(x)**2 != x
    assert all(jacobi(d, x) == 1 for d in moduli[:3])
    print("Three-row one-coordinate false survivor, with factors:", witness)
