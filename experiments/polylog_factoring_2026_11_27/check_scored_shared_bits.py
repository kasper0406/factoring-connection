"""Check mass ordering, exhaustive residues, and noisy-list completion against true factors."""

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb
from random import Random

from sympy import primerange

if __package__:
    from .check_cofactor_zero_queries import direct_zero
    from .hh_residue_lists import factor_from_shared_bits
    from .scored_shared_bits import factor_from_noisy_shared_bits, ranked_shared_residues
else:
    from check_cofactor_zero_queries import direct_zero
    from hh_residue_lists import factor_from_shared_bits
    from scored_shared_bits import factor_from_noisy_shared_bits, ranked_shared_residues


def check_ranking():
    samples = total = 0
    for k in range(2, 7):
        m = 1 << k
        for n in range(1, m, 2):
            for report in range(1, m, 4):
                actual = list(ranked_shared_residues(n, k, report, block_bits=2))
                assert sorted(x for x, _, _ in actual) == list(range(1, m, 2))
                keys = [key for _, key, _ in actual]
                assert keys == sorted(keys)
                assert sum(Fraction(1, key) for key in keys) <= 1
                for x, key, errors in actual:
                    exact = (x & (n*pow(x, -1, m) % m)) ^ report
                    assert exact.bit_count() == errors
                    arities = 1
                    for j in range(1, k, 2):
                        end = min(j+2, k)
                        mod = 1 << end
                        counts = Counter((v & (n*pow(v, -1, mod) % mod)) >> j
                                         for v in range(x % (1 << j), mod, 1 << j))
                        subword = (x & (n*pow(x, -1, mod) % mod) & (mod-1)) >> j
                        arities *= counts[subword]
                    assert key == arities*k*comb(k-1, errors)
                total += len(actual)
                samples += 1
    print("Exhaustive ordered rankings and integer keys:", samples, total)


def check_factoring():
    tests = 0
    rng = Random(2187147)
    for p, q in combinations(primerange(11, 92), 2):
        if q >= 2*p:
            continue
        n, k = p*q, 4
        observed = 1+sum(1<<j for j in range(1, k) if direct_zero(n, 1<<j))
        assert observed == ((p&q) % (1 << k))
        assert factor_from_shared_bits(n, k, observed, complete_list=True) in (p, q)
        for noisy in (observed ^ (1 << rng.randrange(1,k)), rng.randrange(1,1<<k,2)):
            assert factor_from_noisy_shared_bits(n,k,noisy) in (p,q)
            tests += 1
    print("Factoring with altered or random small reports:", tests)


def check_larger():
    p, q = 17179869209, 25769803799
    n, k = p*q, 17
    mask = (p&q) % (1 << k)
    full = {}
    assert factor_from_shared_bits(n, k, mask, complete_list=True, stats=full) in (p,q)
    print("Full exact tree, 69-bit factor-backed report:", full)
    for errors in (0, (1 << 6), (1 << 6) | (1 << 14)):
        work = {}
        assert factor_from_noisy_shared_bits(n, k, mask ^ errors, stats=work) in (p,q)
        print("69-bit fixture, injected errors, leaves, batches, status:",
              errors.bit_count(), work["yielded"], work["batches"], work["completion"]["status"])


if __name__ == "__main__":
    check_ranking()
    check_factoring()
    check_larger()
