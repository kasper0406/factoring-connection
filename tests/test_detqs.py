"""Tests for src/detqs (sieve exactness, elimination, QS driver, Lemma P, local model)."""
import math
import random

import gmpy2
import pytest
import sympy

from src.detqs.linalg import IncrementalF2
from src.detqs.localmodel import _vdist, model_LN, rho, u_blocks, vdist_square_minus_N
from src.detqs.qs import qs_run
from src.detqs.sieve import factor_base, smooth_mask, smooth_relations


def _brute_smooth(v, y):
    for l in range(2, y + 1):
        while v % l == 0:
            v //= l
    return v == 1


@pytest.mark.parametrize("N", [1000003 * 1000033, 2 ** 40 + 15, 123456789011, 3 ** 5 * 7 * 1000003])
def test_sieve_exact(N):
    y, T = 150, 2500
    mask, m, Q, fb = smooth_mask(N, y, T)
    assert m * m >= N > (m - 1) ** 2
    for t in range(1, T + 1):
        assert int(Q[t - 1]) == (m + t) ** 2 - N
        assert bool(mask[t - 1]) == _brute_smooth((m + t) ** 2 - N, y)


def test_relations_reproduce_Q():
    N = 1000003 * 1000033
    m, fb, ts, exps, Q = smooth_relations(N, 300, 20000)
    assert ts
    for t, e in zip(ts, exps):
        assert math.prod(l ** k for l, k in e.items()) == (m + t) ** 2 - N


def test_incremental_f2_kernel():
    rng = random.Random(1)
    rows = [rng.getrandbits(12) for _ in range(20)]
    el = IncrementalF2()
    for i, r in enumerate(rows):
        c = el.add(r)
        if c is not None:
            acc = 0
            for j in range(i + 1):
                if (c >> j) & 1:
                    acc ^= rows[j]
            assert acc == 0 and (c >> i) & 1
    assert el.rank <= 12


@pytest.mark.parametrize("bits", [32, 40, 48])
def test_qs_factors_and_jacobi(bits):
    rng = random.Random(bits)
    for _ in range(4):
        while True:
            p = int(gmpy2.next_prime(rng.getrandbits(bits // 2) | (1 << (bits // 2 - 1))))
            q = int(gmpy2.next_prime(p + rng.randrange(1, p)))
            if p % 4 == 3 and q % 4 == 3:
                break
        N = p * q
        r = qs_run(N, 1 << (bits // 5), 1 << 17)
        assert r.factor in (p, q)
        for j, o in zip(r.dep_jacobi, r.dep_outcomes):   # Blum rule
            assert (j == -1) == bool(o)


def _sqfree_kernel(v):
    k = 1
    for l, e in sympy.factorint(v).items():
        if e % 2:
            k *= l
    return k


@pytest.mark.parametrize("N", [10007 * 10009, 1009 * 2003, 3 * 5 * 7 * 11 * 13 * 17 * 101 * 103, 1009 * 1013 * 1019])
def test_lemma_P_pair_dependencies_always_fruitful(N):
    """Lemma P: t != t', Q_t Q_t' a square, Q_t != Q_t' => 1 < gcd(X - Y, N) < N."""
    m = math.isqrt(N - 1) + 1
    seen = {}
    pairs = 0
    for t in range(1, 6000):
        v = (m + t) ** 2 - N
        k = _sqfree_kernel(v)
        if k in seen:
            t2 = seen[k]
            v2 = (m + t2) ** 2 - N
            X = (m + t) * (m + t2)
            Y = math.isqrt(v * v2)
            assert Y * Y == v * v2
            g = math.gcd(X - Y, N)
            assert 1 < g < N, (N, t, t2)
            pairs += 1
        else:
            seen[k] = t
    assert pairs >= 1


def test_local_model_basic():
    assert abs(rho(2) - (1 - math.log(2))) < 1e-6
    for l in (3, 5, 7):
        for N in (1, 2, 3, 4, 10, 11):
            d = _vdist(l, N, 6)
            assert abs(sum(d) - 1) < 1e-12
            if N % l:
                e = vdist_square_minus_N(l, N, 5)
                # closed form agrees with enumeration on e = 0..3
                for k in range(4):
                    assert abs(d[k] - e[k]) < 1e-12


def test_lmodel_fast_matches_reference():
    from src.detqs.harness import LModel
    n, y, T = 40, 256, 1 << 16
    lm = LModel(n, y, T)
    us, ws = u_blocks(n, y, T, 8)
    for N in (2 ** 39 + 7, 2 ** 39 + 1234567, 3 * 2 ** 37 + 11):
        assert abs(lm.L(N) - model_LN(N, y, us, ws)) < 1e-9 * max(1.0, lm.L(N))
