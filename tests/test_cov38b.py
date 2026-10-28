"""Tests for src/cov38b (harness code)."""
from fractions import Fraction
from math import isqrt, gcd, lcm

import numpy as np
import pytest
import sympy

from src.cov38b import design as ds
from src.cov38b import fourier as fo
from src.dchworst import coverage as cv


def test_y_of_exact():
    for q in (100003, 10 ** 9 + 7, 1300000003):
        Y = ds.y_of(q)
        assert Y ** 25 <= q ** 8 < (Y + 1) ** 25


def test_mu_table_thue_and_diagonal():
    for l in (43, 103, 179, 787):
        tab = ds.mu_table(l)
        assert tab[1] == 1 and tab[0] == ds.INF
        assert tab[1:].max() <= isqrt(l)
        # brute force for a few t
        for t in (2, 5, l - 1, l // 3):
            best = min(max(abs(A), B) for B in range(1, l) for A in range(-l, l + 1)
                       if A != 0 and gcd(A, B) == 1 and (A - t * B) % l == 0 and max(abs(A), B) <= isqrt(l) + 1)
            assert tab[t] == best


def test_heights_match_bruteforce_coverage():
    q = 100003
    Y = ds.y_of(q)
    L = ds.band(Y)
    tabs = {l: ds.mu_table(l) for l in L}
    ps = [p for p in sympy.primerange(q // 2 + 1, q)][:400]
    h = ds.heights(np.array(ps, dtype=np.int64), q, tabs)
    W = cv.weights(isqrt(2 * Y) + 1, True) + [(1, 1)]
    for p, hp in zip(ps, h):
        hb = min(max(abs(A), B) for (A, B) in W for l in L if (B * p - A * q) % l == 0)
        assert hb == hp


def test_primes_in_segmented():
    got = np.concatenate(list(ds.primes_in(1000, 50000, chunk=4096)))
    assert got.tolist() == list(sympy.primerange(1001, 50000))


def test_reciprocal_test_implies_coverage():
    # coverage38 sec 3.1: B <= H/4, |B p lbar mods q| <= qH/(4Y)  =>  covered at height H
    q = 1000003
    Y = ds.y_of(q)
    L = ds.band(Y)
    H = 8
    k0 = q * H // (4 * Y)
    tabs = {l: ds.mu_table(l) for l in L}
    hits = 0
    for p in list(sympy.primerange(q // 2 + 1, q))[:3000]:
        ok = False
        for l in L:
            x = p * pow(l, -1, q) % q
            for B in range(1, H // 4 + 1):
                k = B * x % q
                if min(k, q - k) <= k0:
                    ok = True
        if ok:
            hits += 1
            assert ds.heights(np.array([p], dtype=np.int64), q, tabs)[0] <= H
    assert hits > 0


def test_freiman_energy_identity():
    # Lemma F0: E(S) mod q = E_Q(R) when 4 H'^3 k0 < q
    for (Hp, k0, q) in ((2, 5, 1009), (3, 6, 10007), (4, 7, 10007), (5, 3, 10007)):
        assert 4 * Hp ** 3 * k0 < q
        assert fo.energy_exact(fo.farey_rationals(Hp, k0)) == fo.energy_mod_q(q, Hp, k0)
        # and |S| = |R|
        S = {(pow(B, -1, q) * k) % q for B in range(1, Hp + 1) for k in range(-k0, k0 + 1)}
        assert len(S) == len(fo.farey_rationals(Hp, k0))


def test_lemma_F1_energy_bound():
    # E(R) <= 8 (2k0+1)^3 * sum_{B} min(B)/lcm(B)   whenever 2k0+1 >= H'^3
    for Hp, k0 in ((2, 4), (3, 13), (4, 32)):
        assert 2 * k0 + 1 >= Hp ** 3
        E = fo.energy_exact(fo.farey_rationals(Hp, k0))
        assert E <= 8 * (2 * k0 + 1) ** 3 * fo.lcm_sum(Hp)


def test_lcm_sum_divisor_bound():
    # sum min/lcm <= H * sum_{M <= H^3} tau(M)^4 / M   (proof of Lemma F1)
    for H in (2, 3, 5, 7):
        rhs = H * sum(Fraction(sympy.divisor_count(M) ** 4, M) for M in range(1, H ** 3 + 1))
        assert fo.lcm_sum(H) <= rhs
        assert fo.lcm_sum(H) >= H  # diagonal B1=..=B4


def test_farey_builders_agree():
    q = 10007
    for Hp, k0 in ((1, 50), (3, 200), (6, 30)):
        for sh in ("tri", "sharp"):
            assert np.abs(fo.farey_tri(q, Hp, k0, sh) - fo.farey_direct(q, Hp, k0, sh)).max() == 0


def test_certificate_is_valid_lower_bound():
    q = 10007
    Y = ds.y_of(q)
    L = ds.band(3 * Y)  # a few more primes
    J = len(L)
    K = fo.kloosterman_all(q, L)
    kappa = K[1:].max()
    g = fo.farey(q, 3, 400, "tri")
    st = fo.stats(g)
    lower = fo.certificate(J, kappa, st, q)
    linv = fo.inv_band(q, L)
    for a in range(1, q, 97):
        cnt = g[(a * linv) % q].sum()
        assert cnt >= lower - 1e-6


def test_holder_lower_bound_le_wiener():
    q = 10007
    for Hp, k0 in ((1, 300), (4, 40)):
        st = fo.stats(fo.farey(q, Hp, k0, "tri"))
        assert st["LB"] <= st["W"] + 1 + 1e-9


def test_mobius_decomposition_of_farey_set():
    # Lemma F2: 1_S = sum_{d<=h} sum_{e|d} mu(e) 1[(d/e)^{-1} [-K/e, K/e]]  (2hK < q)
    q, h, K = 10007, 6, 300
    assert 2 * h * K < q
    S = fo.farey(q, h, K, "sharp")
    acc = np.zeros(q)
    for d in range(1, h + 1):
        for e in sympy.divisors(d):
            mu = sympy.mobius(e)
            if mu == 0:
                continue
            m = d // e
            j = np.arange(-(K // e), K // e + 1, dtype=np.int64)
            np.add.at(acc, (pow(m, -1, q) * j) % q, mu)
    assert np.array_equal(acc, S)


def test_inverse_prime_energy_minimal():
    # Prop B': E({lbar}) = 2J^2 - J when 16 Y^3 < q
    q = 10000019
    Y = 20
    L = [l for l in sympy.primerange(Y + 1, 2 * Y + 1)]
    assert 16 * Y ** 3 < q
    inv = [pow(l, -1, q) for l in L]
    from collections import Counter
    c = Counter((a + b) % q for a in inv for b in inv)
    J = len(L)
    assert sum(v * v for v in c.values()) == 2 * J * J - J
