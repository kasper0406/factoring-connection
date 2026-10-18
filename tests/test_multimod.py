"""Tests for src/multimod (DCL several-moduli search)."""
from math import gcd, isqrt
import random

import pytest

from src.common.instances import balanced_semiprime
from src.multimod.dcl import (best_design, candidate_range, design, f_range, primes_upto,
                              sqrt_mod_prime, target_window, uncovered_candidates,
                              verify_windows)
from src.multimod.operational import all_but_one_powers, dcl_search
from src.common.metrics import Metrics


def test_sqrt_mod_prime():
    for M in primes_upto(200)[1:]:
        for a in range(M):
            roots = sqrt_mod_prime(a, M)
            assert sorted(roots) == sorted({x for x in range(M) if x * x % M == a})


def test_class_lock_and_window():
    """M | Bp - Aq  =>  Bp + Aq == c (mod M^2) and g in window (true p,q)."""
    hits = 0
    for seed in range(6):
        inst = balanced_semiprime(40, seed)
        N, p = inst.N, inst.labels["p"]
        q = N // p
        for (A, B) in [(1, 1), (1, 2), (2, 3), (3, 1)]:
            F_lo, F_hi = f_range(N, A, B)
            assert F_lo <= B * p + A * q <= F_hi
            for M in primes_upto(400)[1:]:
                if (A * B * N) % M == 0 or (B * p - A * q) % M:
                    continue
                r = p % M
                assert (B * r * r - A * N) % M == 0
                c, g_lo, W = target_window(N, A, B, M, r, F_lo, F_hi)
                assert (B * p + A * q - c) % (M * M) == 0
                g = (B * p + A * q - c) // (M * M)
                assert g_lo <= g < g_lo + W
                hits += 1
    assert hits > 20


def test_all_but_one_powers():
    N = 1000003 * 999983
    exps = [5, 49, 121, 13 * 13, 289]
    met = Metrics()
    got = all_but_one_powers(7, exps, N, met)
    tot = 1
    for e in exps:
        tot *= e
    assert got == [pow(7, tot // e, N) for e in exps]


def test_design_windows_exact():
    N = 2 ** 40 + 2 ** 20 + 1
    d = best_design(N, 3.0)
    assert d["lam"] >= 3.0
    assert verify_windows(N, d["targets"]) == 0
    assert d["cost"] == d["m"] + sum(t["k"] for t in d["targets"])


def test_coverage_sieve_matches_true_p():
    """Harness: the sieve says p covered iff p is in some target class."""
    for seed in range(4):
        inst = balanced_semiprime(36, seed)
        N, p = inst.N, inst.labels["p"]
        d = best_design(N, 3.0)
        lo, hi, unc = uncovered_candidates(N, d["targets"])
        cov = any(p % t["M"] == t["r"] for t in d["targets"])
        assert lo <= p <= hi
        assert cov == (not unc[p - lo])


@pytest.mark.parametrize("bits", [32, 40])
def test_dcl_search_factors(bits):
    ok = 0
    for seed in range(5):
        inst = balanced_semiprime(bits, seed)
        res = dcl_search(inst.N, 3.0, seed, early_exit=False)
        ok += res["factor"] in (inst.labels["p"], inst.N // inst.labels["p"])
    assert ok >= 4


def test_lower_bound_exponents():
    """Arithmetic of Theorems A/B (notes/multimod_results.md section 3)."""
    from fractions import Fraction as F
    # e = 1, per-target coverage (PW)^{1/2}, sum W <= X^2: X >= P^{1/3}
    assert F(1, 3) * F(1, 2) == F(1, 6)
    # e >= 2 branch: e/(K+1) >= 2/5 for all K <= 2e
    for e in range(2, 8):
        for K in range(1, 2 * e + 1):
            assert F(e, K + 1) >= F(2, 5)


def test_design_reproducible_from_m():
    """dcl_search(m=...) rebuilds exactly the best_design target list."""
    from src.multimod.dcl import design
    import math
    for n in (40, 56):
        N = 2 ** n + 2 ** (n // 2) + 1
        for lam in (3.0, math.log(math.isqrt(N)) + 2):
            d = best_design(N, lam)
            d2 = design(N, lam, d["m"])
            assert d2["targets"] == d["targets"]
