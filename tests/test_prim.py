"""Tests for src/prim (phase 25)."""
import random
from math import gcd

import numpy as np
import pytest
import sympy

from src.prim.prog import (C, G, FG, LucasNode, NestModel, PPow, control_nodes, nested_nodes, summed_nodes, vpow)
from src.prim.zerosets import (coset_union, ks_bound, ks_threshold, max_vanishing_coset, primes_b)


def test_vpow_matches_pow():
    p = 10007
    rng = np.random.default_rng(1)
    b = rng.integers(0, p, 500)
    e = rng.integers(0, p - 1, 500)
    out = vpow(b, e, p)
    for x, y, z in zip(b, e, out):
        exp = 0 if x == 0 else pow(int(x), int(y), p)
        assert z == exp


def test_uniform_r_identity_theorem_U():
    """At the true partner q the gate value mod p equals the uniform-r value at r = q mod (p-1)."""
    ps = primes_b(11)
    nodes = [nd for nd, _ in nested_nodes().values() if nd.affine()] + [nd for nd, _ in summed_nodes().values()]
    nodes += list(control_nodes().values())
    rng = random.Random(3)
    for _ in range(40):
        p, q = rng.sample(ps, 2)
        N = p * q
        for nd in nodes:
            assert nd.ev(N, p) == int(nd.vec(p, np.array([q % (p - 1)]))[0])


def test_residue_hit_at_most_twice():
    """each residue mod p-1 is hit by <= 2 primes of the range (Theorem U counting step)."""
    ps = primes_b(12)
    for p in ps[::7]:
        c = {}
        for q in ps:
            if q != p:
                c[q % (p - 1)] = c.get(q % (p - 1), 0) + 1
        assert max(c.values()) <= 2


def test_diagonal_lemma_nested():
    ps = primes_b(10)
    nd, _ = nested_nodes()["D2a_1"]
    for p in ps:
        for q in ps:
            if q <= p:
                continue
            N = p * q
            v = nd.ev(N, N)
            zp, zq = nd.ev(N, p) == 0, nd.ev(N, q) == 0
            g = gcd(v, N)
            assert (1 < g < N) == (zp != zq)
            if zp == zq:
                assert (v == 0) == zp


def test_model_exact_for_uniform_u():
    """NestModel probability equals the exact average over u in F_p^* of [A u^e = B]."""
    p = 211
    for e in (0, 1, 2, 6, 7, 30, 105, 210, 211):
        for A, B in ((1, 1), (3, 5), (2, 4), (7, 0)):
            exact = sum(1 for u in range(1, p) if A * pow(u, e, p) % p == B % p) / (p - 1)
            m = NestModel(C(A), C(1), PPow(C(1), e if e else p - 1), C(B))
            assert abs(m.prob_scalar(1, p) - exact) < 1e-12


def test_coset_union_and_max():
    n = 12
    Z = np.zeros(n, dtype=bool)
    Z[[1, 4, 7, 10]] = True          # coset 1 + 3Z/12 (size 4)
    Z[0] = True
    V = coset_union(Z, 4)
    assert set(np.flatnonzero(V)) == {1, 4, 7, 10}
    assert max_vanishing_coset(Z) == 4


def _sparse_zero_set(p, coeffs, exps):
    r = np.arange(p - 1)
    g = sympy.primitive_root(p)
    gr = vpow(np.full(p - 1, g), r, p)
    val = np.zeros(p - 1, dtype=np.int64)
    for c, e in zip(coeffs, exps):
        val = (val + c * vpow(gr, np.full(p - 1, e % (p - 1)), p)) % p
    return val == 0


def test_lemma_KS_bruteforce():
    """Lemma KS on random sparse polynomials (and on coset-rich ones built from small-order structure)."""
    rng = random.Random(7)
    for p in (101, 241, 421, 601, 1201, 2311):
        n = p - 1
        for _ in range(60):
            t = rng.randint(2, 6)
            if rng.random() < 0.5:
                m = rng.choice(sympy.divisors(n))
                exps = [rng.randrange(n // m) * m + rng.randrange(2) for _ in range(t)]
            else:
                exps = [rng.randrange(n) for _ in range(t)]
            coeffs = [rng.randrange(1, p) for _ in range(t)]
            Z = _sparse_zero_set(p, coeffs, exps)
            V = coset_union(Z, ks_threshold(n))
            assert (Z & ~V).sum() <= ks_bound(n, t)


def test_lucas_node_consistent():
    ps = primes_b(10)
    for p, q in [(ps[0], ps[5]), (ps[3], ps[9])]:
        N = p * q
        vN = LucasNode.ev(N, N)
        assert (vN[0] % p, vN[1] % p) == LucasNode.ev(N, p)
        # V_N = alpha^N + alpha^{-N}: check against the Lucas recurrence for small N
    for n_ in range(1, 30):
        V = [2, 3]
        for _ in range(n_):
            V.append(3 * V[-1] - V[-2])
        assert LucasNode.lucas_V(n_, 10 ** 9 + 7) == V[n_] % (10 ** 9 + 7)


def test_floor_gate():
    N = 3 * 7 * 1
    assert FG(C(2), 1, 1).E(21) == 10
    assert FG(C(2), -1, 1).E(21) == 11


# ----------------------------------------------------------------------------- theorem companions
from src.prim.theory import collision_check, d2_bound, theorem_r_check
from src.prim.zerosets import zero_set


def test_theorem_D2_bound_holds():
    """Theorem D2 zero-set bound for (2^N + 1)^N - w and 2^N (2^N + 3)^N - 5 at every 12-bit prime."""
    nest = nested_nodes()
    for p in primes_b(12):
        for name, (nd, _) in nest.items():
            if name not in ("D2a_1", "D2a_2", "D2a_3", "D2d"):
                continue
            Z = zero_set(nd, p)
            bnd, _ = d2_bound(p, 2, 1, 1)
            assert Z.sum() <= bnd


def test_collision_bound():
    rng = random.Random(5)
    for p in (101, 211, 401, 1009):
        for _ in range(30):
            d = rng.randint(1, 4)
            exps = [rng.randrange(1, p - 1) for _ in range(d)]
            shifts = [rng.randrange(1, p) for _ in range(d - 1)]
            pr, bound = collision_check(p, exps, shifts, rng.randrange(1, p))
            assert pr <= bound + 1e-12


def test_theorem_R_exact():
    primes = list(sympy.primerange(100, 200))
    probes = [([lambda N: N - 1], [], 1),
              ([lambda N: N - 1, lambda N: N + 1], [1], 2),
              ([lambda N: (N - 1) // 2, lambda N: N, lambda N: N - 1], [3, 5], 1)]
    for target in (lambda p, q: int(p % 4 == 1), lambda p, q: int(p % 3 == 1), lambda p, q: int(p < q)):
        r = theorem_r_check(primes, probes, target)
        assert r["holds"], r
