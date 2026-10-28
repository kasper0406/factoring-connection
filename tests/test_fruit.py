"""Tests for src/fruit (phase 30).  Factors appear only in test-side checks (flip test)."""
import itertools
import random
from math import gcd

import gmpy2
import pytest

from src.fruit.linalg import hybrid_rank, kernel_basis, sample_forced, sample_uniform
from src.fruit import nfs_hybrid as nh
from src.fruit.harness_e1 import QSDegenerate, beta_zero, instances
from src.fruit.harness_e2 import rational_square_pairs, two_term, sqrt_branch_check
from src.fruit.qs_hybrid import fruitful
from src.fruit.theory import (exponent_hybrid, exponent_pure, count_rational_square_m,
                              rational_square_bound)
from src.nfsfruit.instances import composite


def _span(basis):
    out = set()
    for bits in itertools.product((0, 1), repeat=len(basis)):
        S = 0
        for b, t in zip(basis, bits):
            if t:
                S ^= b
        out.add(S)
    return out


def test_lemma_HY_exact_fraction():
    rng = random.Random(1)
    for _ in range(60):
        C = rng.randrange(3, 9)
        M = rng.randrange(1, 9)
        h = rng.randrange(1, 4)
        A = [rng.getrandbits(C) for _ in range(M)]
        D = [rng.getrandbits(C) for _ in range(h)]
        hr = hybrid_rank(A, D)
        K = _span(kernel_basis(A + D))
        frac = sum(1 for S in K if S >> M) / len(K)
        assert frac == 1 - 2 ** (-hr["r"])
        assert max(0, h - hr["deficit"]) <= hr["r"] <= h


def test_forced_sampler_uses_desc():
    rng = random.Random(2)
    A = [0b011, 0b110, 0b101]
    D = [0b011]
    kern = kernel_basis(A + D)
    for _ in range(20):
        S = sample_forced(kern, 3, rng)
        assert S >> 3 == 1


def test_qs_degenerate_base_is_fruitless_and_descent_flip():
    seed, N, fs = instances("blum1", 40, 1)[0]
    s = QSDegenerate(N, 300, 200000, 2000)
    s.add_qs_beta0(10 ** 7, 80)
    rows = list(s.rows)
    kern = kernel_basis(rows)
    assert kern
    for S in kern:
        X, Y = s.congruence(S)
        assert (X * X - Y * Y) % N == 0
        assert not fruitful(X, Y, N)      # nfsfruit Cor J1: Phi == 0 on beta = 0 rows
    rng = random.Random(3)
    assert s.add_descent_exact(rng, 100000)
    M = len(rows)
    kern2 = kernel_basis(s.rows)
    S = sample_forced(kern2, M, rng)
    if S is None:
        pytest.skip("descent row not in rowspace for this toy")
    X, Y = s.congruence(S)
    # flip check (test-only use of the factors): r -> eps r flips the fruitfulness class
    p, q = fs
    eps = (q * pow(q, -1, p) * (p - 1) + p * pow(p, -1, q)) % N   # -1 mod p, +1 mod q
    i = M
    s.X[i] = s.X[i] * eps % N
    X2, Y2 = s.congruence(S)
    assert Y2 == Y and fruitful(X, Y, N) != fruitful(X2, Y2, N)


def test_beta_zero_formula_matches_jacobi():
    N = 7 * 11  # Blum, N = 5 mod 8
    assert N % 8 == 5
    for x in range(1, 40):
        if gcd(x, N) == 1:
            b = (gmpy2.jacobi(x, N) == -1) ^ bool(x & 1)
            assert beta_zero(x, N) == (not b)


def test_nfs_hybrid_congruence_and_flip():
    N, fs = composite("blum", 40, 3)
    s = nh.build_base(N, 3, 1500, 1000, 400, 40, 60, 12000, 48, 2)
    rng = random.Random(4)
    for _ in range(2):
        assert s.add_descent(rng, 20000)
    s.finalize_hybrid()
    M = len(s.ab)
    kern = kernel_basis(s.rows + s.drows)
    done = 0
    for _ in range(10):
        S = sample_forced(kern, M, rng)
        if S is None:
            break
        bi = [i for i in range(M) if S >> i & 1]
        di = [j for j in range(len(s.drows)) if S >> (M + j) & 1]
        out = s.hybrid_congruence(bi, di)
        if out is None:
            continue
        X, Y = out
        assert (X * X - Y * Y) % N == 0
        ok, nroots = sqrt_branch_check(s, bi, di)
        assert ok
        # flip one used hidden root (test-only eps)
        p, q = fs
        eps = (q * pow(q, -1, p) * (p - 1) + p * pow(p, -1, q)) % N
        j = di[0]
        r0 = s.desc[j]["r"]
        s.desc[j]["r"] = r0 * eps % N
        X2, Y2 = s.hybrid_congruence(bi, di)
        s.desc[j]["r"] = r0
        assert Y2 == Y and (X2 * X2 - Y2 * Y2) % N == 0
        assert fruitful(X, Y, N) != fruitful(X2, Y2, N)
        done += 1
    assert done >= 1


def test_two_term_counters():
    rows = [5, 3, 5, 9, 3, 5]
    pairs, members = two_term(None, None, rows)
    assert pairs == 3 + 1 and members == 5
    rp, rm = rational_square_pairs([12, 3, -3, 27], [{2: 2, 3: 1}, {3: 1}, {3: 1}, {3: 3}])
    assert rp == 3 and rm == 3


def test_exponent_formulas():
    from fractions import Fraction as Fr
    assert [exponent_pure(k) for k in (1, 2, 3)] == [Fr(2, 5), Fr(3, 8), Fr(4, 11)]
    assert [exponent_hybrid(k) for k in (1, 2, 3, 4)] == [Fr(1, 2), Fr(4, 9), Fr(5, 12), Fr(2, 5)]
    for k in range(1, 30):
        assert exponent_hybrid(k) > exponent_pure(k) > Fr(1, 3)


def test_rational_square_bound_small():
    # exact count of m in [L, 2L) for which two fixed box pairs give a rational square product
    for (a, b, a2, b2) in [(1, 5, 2, 7), (0, 3, 1, 4), (2, 9, 5, 8)]:
        L = 4000
        cnt = count_rational_square_m(a, b, a2, b2, L, 2 * L)
        assert cnt <= rational_square_bound(a, b, a2, b2, L, 2 * L)


def test_chain_feasibility_matches_formulas():
    from fractions import Fraction as Fr
    from src.fruit.theory import chain_feasible
    for k in range(1, 8):
        for mode, f in (("pure", exponent_pure), ("hybrid", exponent_hybrid)):
            a = f(k)
            assert chain_feasible(a, k, mode)
            assert not chain_feasible(a - Fr(1, 10 ** 6), k, mode)


def test_degree_t_ladders():
    from fractions import Fraction as Fr
    from src.fruit.theory import chain_feasible_t, exponent_hybrid_t, exponent_pure_t
    assert [exponent_hybrid_t(k) for k in (1, 2, 3)] == [Fr(4, 9), Fr(8, 21), Fr(16, 45)]
    assert [exponent_pure_t(k) for k in (1, 2, 3)] == [Fr(2, 5), Fr(4, 11), Fr(8, 23)]
    for k in range(1, 12):
        for mode, f in (("hybrid", exponent_hybrid_t), ("pure", exponent_pure_t)):
            a = f(k)
            assert chain_feasible_t(a, k, mode)
            assert not chain_feasible_t(a - Fr(1, 10 ** 9), k, mode)
        assert exponent_hybrid_t(k) > exponent_pure_t(k) > Fr(1, 3)
        assert exponent_hybrid_t(k) <= exponent_hybrid(k) and exponent_pure_t(k) <= exponent_pure(k)
