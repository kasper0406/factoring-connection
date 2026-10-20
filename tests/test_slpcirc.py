"""Tests for src/slpcirc (Lemma S1 / Theorem SP checks, enumeration, classification, E4 fibres)."""
import itertools
import random
from math import gcd, lcm

import flint
import pytest
from sympy import primerange, factorint

from src.slpcirc.enum_slp import enumerate_tau, pmul, padd
from src.slpcirc.sparse import (s1_coset_bound, primitive_root, Lambda, dichotomy_violations,
                                count_zeros_Fp_star, ell0)
from src.slpcirc.e1_e3 import roots_count_sparse, roots_count_sparse_np
from src.slpcirc.classify import classify_factor
from src.slpcirc.quadcomp import fibre, sqrt_table, all_fibres_depth


def test_enumeration_small():
    tau = enumerate_tau(3)
    assert tau[(1,)] == 0 and tau[(0, 1)] == 0
    assert tau[(2,)] == 1
    assert tau[(0, 0, 1)] == 1          # x^2
    assert tau[(0, 0, 0, 0, 1)] == 2    # x^4
    assert tau[(-1, 0, 1)] == 2         # x^2 - 1
    assert tau[tuple([0] * 8 + [1])] == 3  # x^8
    assert len(tau) == 185


def test_s1_univariate_bruteforce():
    rng = random.Random(3)
    for p in primerange(5, 90):
        n = p - 1
        for _ in range(40):
            t = rng.choice([2, 3, 4])
            if n <= t:
                continue
            exps = [0] + rng.sample(range(1, n), t - 1)
            ce = [(1, 0)] + [(rng.randrange(1, p), e) for e in exps[1:]]
            z = count_zeros_Fp_star(ce, p)
            for ell in factorint(n):
                if len({e % ell for e in exps}) == t:
                    assert z <= (n // ell) * s1_coset_bound(ell, t)


def test_s1_multivariate_bruteforce():
    """Lines argument: k = 2, zeros on (F_p^*)^2 bounded when ell > C(t,2) and exponents distinct mod ell."""
    rng = random.Random(5)
    for p in (31, 37, 43, 61, 67):
        n = p - 1
        g = primitive_root(p)
        pw = [pow(g, j, p) for j in range(n)]
        for _ in range(25):
            t = rng.choice([2, 3])
            vecs = [(rng.randrange(n), rng.randrange(n)) for _ in range(t)]
            if len(set(vecs)) < t:
                continue
            cs = [rng.randrange(1, p) for _ in range(t)]
            z = 0
            for i in range(n):
                for j in range(n):
                    if sum(c * pw[(a * i + b * j) % n] for c, (a, b) in zip(cs, vecs)) % p == 0:
                        z += 1
            for ell in factorint(n):
                if ell > t * (t - 1) // 2 and len({(a % ell, b % ell) for a, b in vecs}) == t:
                    assert z <= (n * n // ell) * s1_coset_bound(ell, t)


def test_dichotomy_pollard():
    M = lcm(*range(1, 60))
    for p in primerange(2 ** 12, 2 ** 12 + 3000):
        ce = [(1, M), (-1, 0)]
        r = roots_count_sparse(ce, p)
        assert r == gcd(M, p - 1) + 0  # x = 0 is not a root of x^M - 1
        assert dichotomy_violations(ce, p, 1, 20, r) == []


def test_roots_count_np_agrees():
    rng = random.Random(11)
    for _ in range(100):
        p = rng.choice([101, 103, 1009, 4099])
        ce = [(rng.choice([-3, -1, 1, 2]), rng.randrange(10 ** 9)) for _ in range(3)] + [(1, 0)]
        assert roots_count_sparse(ce, p) == roots_count_sparse_np(ce, p)
        br = sum(1 for x in range(p) if sum(c * pow(x, e, p) for c, e in ce) % p == 0)
        assert br == roots_count_sparse(ce, p)


def test_dickson_root_formula():
    """D_m(x,1) - 2 has (gcd(m,p-1)+gcd(m,p+1))/2 roots in F_p (u <-> u^{-1} on the two tori)."""
    for m in (3, 4, 5, 6, 8, 12):
        T = flint.fmpz_poly.chebyshev_t(m)  # T_m(x/2)*2 = D_m(x)
        coeffs = [int(c) for c in T.coeffs()]
        # D_m(x) = 2*T_m(x/2)
        for p in primerange(7, 120):
            if p <= m:
                continue
            inv2 = pow(2, -1, p)
            cnt = 0
            for x in range(p):
                y = x * inv2 % p
                v = sum(c * pow(y, k, p) for k, c in enumerate(coeffs)) % p
                if (2 * v - 2) % p == 0:
                    cnt += 1
            assert cnt == (gcd(m, p - 1) + gcd(m, p + 1)) // 2


def test_classify():
    P = flint.fmpz_poly
    assert classify_factor(P([1, 1, 1])) == "CYC"
    assert classify_factor(P([1, 0, 4])) == "CYC"          # Phi_4(2x)
    assert classify_factor(P([-2, 0, 0, 1])) == "KUM"
    assert classify_factor(P.cos_minpoly(7)) == "CHEB"
    assert classify_factor(P([-1, 1])) == "LIN"
    assert classify_factor(P([-1, -1, 0, 0, 0, 0, 1])) in ("OTHER", "LOW")


def test_quadcomp_fibre_matches_polynomial():
    rng = random.Random(2)
    for p in (101, 211, 307):
        sq = sqrt_table(p)
        for _ in range(20):
            seq = tuple(rng.choice((-2, -1, 0, 1, 2)) for _ in range(rng.randrange(1, 6)))
            f = flint.nmod_poly([0, 1], p)
            for c in seq:
                f = f * f + c
            assert len(f.roots()) == fibre(seq, p, sq)
        tab = all_fibres_depth(3, p, sq)
        assert len(tab) == 125 and all(tab[s] == fibre(s, p, sq) for s in tab)


def test_ell0_monotone():
    assert ell0(2, 1, 20) <= ell0(3, 1, 20) <= ell0(4, 1, 20)
    assert s1_coset_bound(ell0(3, 1, 20), 3) * 20 < ell0(3, 1, 20)
