"""Tests for src/below16 (agent below16)."""
import math
import random

import gmpy2
import numpy as np
import pytest

from src.below16 import avgbound as AB
from src.below16 import multibase as MB
from src.below16 import qs_pilot as QS
from src.below16 import experiments as EX


def brute_population(n):
    out = set()
    for p in range(2, 1 << (n // 2 + 1)):
        if not gmpy2.is_prime(p):
            continue
        for q in range(p + 1, 2 * p + 1):
            if gmpy2.is_prime(q) and (1 << (n - 1)) <= p * q < (1 << n):
                out.add((p, q))
    return out


@pytest.mark.parametrize("n", [12, 13, 14])
def test_population(n):
    p, q = AB.population(n)
    assert set(zip(p.tolist(), q.tolist())) == brute_population(n)
    xmin, xmax, _ = AB.class_var_range(n, "p")
    assert p.min() >= xmin and p.max() <= xmax
    ymin, ymax, _ = AB.class_var_range(n, "q")
    assert q.min() >= ymin and q.max() <= ymax


def test_sqrt_mod_prime():
    for l in [3, 5, 7, 11, 13, 17, 97, 101, 1009]:
        for a in range(1, l):
            r = AB.sqrt_mod_prime(a, l)
            if pow(a, (l - 1) // 2, l) == 1:
                assert r * r % l == a
            else:
                assert r is None


@pytest.mark.parametrize("n", [18, 20])
def test_t1_bounds_random_and_dcl(n):
    p, q = AB.population(n)
    rng = random.Random(7 + n)
    gens = [EX.random_t1(rng, n, p, q)[0] for _ in range(25)]
    gens += [EX.oracle_t1(rng, n, p, q, s) for s in (1, 4, 16)]
    gens += AB.dcl_generators(n, 2, math.ceil(2 ** (n / 6) / 2 ** (1 / 3)))[:6]
    for g in gens:
        c = AB.count_T1(g, p, q)
        core = AB.bound_T1_core(g, n)
        assert c <= core
        if AB.t1_conditions(g, n):
            assert core <= AB.bound_T1_closed(g, n) + 1e-6


def test_count_t1_matches_bruteforce():
    n = 16
    p, q = AB.population(n)
    rng = random.Random(3)
    for _ in range(10):
        g, _ = EX.random_t1(rng, n, p, q)
        c = 0
        for a, b in zip(p.tolist(), q.tolist()):
            x, y = (a, b) if g.side == "p" else (b, a)
            k = a * b % g.K
            if k in g.R and k in g.S and x % g.M in g.R[k] and (g.B * x + g.A * y) // g.D in g.S[k]:
                c += 1
        assert c == AB.count_T1(g, p, q)


def test_t2_bounds():
    n = 18
    p, q = AB.population(n)
    rng = random.Random(11)
    for W in (1, 8):
        g = EX.oracle_t2(rng, n, p, q, W)
        c = AB.count_T2(g, n, p, q)
        assert c <= AB.bound_T2_core(g, n) <= AB.bound_T2_closed(g, n) + 1e-6
        for key in list(g.R)[:20]:
            assert AB.cplus_T2(g, n, key) <= AB.cplus_closed(g, n, key) + 1e-9


def test_multibase_superset():
    n = 16
    p, q = AB.population(n)
    rng = random.Random(5)
    for bases in ([2, 3, 6], [2, 4, 8], [6, 10, 15]):
        primes, E = MB.prime_exponents(bases)
        for _ in range(5):
            forms = [(rng.randint(-2, 2), rng.randint(-2, 2)) for _ in bases]
            W = [rng.randint(1, 4) for _ in bases]
            i0 = rng.randrange(len(p))
            c = [b * int(p[i0]) + a * int(q[i0]) for (b, a) in forms]
            cov = MB.covered_multibase(p, q, forms, c, W, E)
            assert cov[i0]
            for j in range(len(primes)):
                cj, sj = MB.covered_coordinate(p, q, forms, c, W, E, j)
                assert np.all(cj | ~cov)
                assert sj <= len(MB.image_of_box(E, W)) <= int(np.prod(W))


def _rp(rng, b, blum=False):
    while True:
        x = rng.getrandbits(b) | (3 << (b - 2)) | 1
        if blum and x % 4 != 3:
            continue
        if gmpy2.is_prime(x):
            return x


def test_qs_deterministic_and_jacobi_rule():
    rng = random.Random(9)
    for n in (40, 48):
        for _ in range(4):
            p, q = _rp(rng, n // 2, True), _rp(rng, n // 2, True)
            if p == q:
                continue
            r1 = QS.qs_factor(p * q, QS.default_config(n))
            r2 = QS.qs_factor(p * q, QS.default_config(n))
            assert r1 == r2
            assert r1["factor"] in (p, q)
            for d in r1["deps"]:
                assert d["nontriv"] == (d["jacobi"] == -1)


def test_qs_kernel():
    rows = [0b011, 0b110, 0b101, 0b111]
    ker = QS._kernel(rows, 4)
    for comb in ker:
        v = 0
        for i in range(4):
            if comb >> i & 1:
                v ^= rows[i]
        assert v == 0 and comb
    assert len(ker) == 4 - 3  # rank 3
