"""Tests for src/dchworst/wcalg.py (DCL-wc, Theorem R).  Tests may use p, q (harness)."""
import random
from math import gcd, isqrt

import gmpy2
import numpy as np
import sympy

from src.common.metrics import Metrics
from src.dchworst import coverage as cv
from src.dchworst import wcalg as wc
from src.dchworst.experiments import blocks
from src.dchworst.experiments_wc import covered_mask


def _orders(b, p):
    return int(sympy.n_order(b, p))


def test_first_r_matches_orders():
    rng = random.Random(5)
    for _ in range(60):
        p = int(gmpy2.next_prime(rng.randint(2000, 40000)))
        q = int(gmpy2.next_prime(rng.randint(p, 2 * p)))
        N = p * q
        for b in (2, 3, 5):
            op, oq = _orders(b, p), _orders(b, q)
            T = rng.choice([50, 500, 5000])
            try:
                kind, r = wc.first_r(b, T, N, Metrics())
            except wc.Found as e:
                assert e.d in (p, q) and op != oq
                continue
            if kind == "equal":
                assert op == oq == r
            else:
                assert min(op, oq) > T


def test_f_bounds_signed():
    rng = random.Random(6)
    for _ in range(200):
        p = int(gmpy2.next_prime(rng.randint(10 ** 5, 10 ** 6)))
        q = int(gmpy2.next_prime(rng.randint(p, 2 * p)))
        N = p * q
        a0, b0 = wc.balanced_range(N)
        assert a0 <= p <= b0
        A = rng.choice([-1, 1]) * rng.randint(1, 9)
        B = rng.randint(1, 9)
        lo, hi = wc.f_bounds(N, A, B, a0, b0)
        assert lo <= B * p + A * q <= hi


def test_square_and_even():
    assert wc.dcl_wc(101 ** 2, [[11]], 2, 64, 5)["factor"] == 101
    assert wc.dcl_wc(2 * 101, [[11]], 2, 64, 5)["factor"] == 2


def _exhaustive(n, design, mode):
    P, Q = cv.all_pairs(n)
    H = 3
    Y = max(3, (1 << (n // 6)) // 2)
    bl = blocks(Y)
    cov = covered_mask(P, Q, bl, H, design)
    T = 1 << (-(-n // 3))
    fails = 0
    for i in range(len(P)):
        p, q = int(P[i]), int(Q[i])
        r = wc.dcl_wc(p * q, bl, H, T, 30, signed=(design == "sgn"), diag=(design == "sgn"),
                      force_base=(2 if mode == "force2" else None))
        if r["factor"] is not None:
            assert (p * q) % r["factor"] == 0 and 1 < r["factor"] < p * q
        elif cov[i]:
            fails += 1
    return fails, int(cov.sum())


def test_theorem_R_exhaustive_n20():
    for design in ("pos", "sgn"):
        for mode in ("full", "force2"):
            fails, ncov = _exhaustive(20, design, mode)
            assert ncov > 0 and fails == 0, (design, mode, fails)


def test_equal_order_branch():
    """p = 2D+1, q = 3D+1 with ord_p(2) = ord_q(2) = D = l0*o: forced base 2 must still
    factor through the equal-order split (or earlier)."""
    found = 0
    Y = 8
    bl = blocks(Y)
    band = [l for b in bl for l in b]
    for l0 in band:
        for o in range(3, 4000):
            D = l0 * o
            p, q = 2 * D + 1, 3 * D + 1
            if not (gmpy2.is_prime(p) and gmpy2.is_prime(q)):
                continue
            if _orders(2, p) != D or _orders(2, q) != D:
                continue
            r = wc.dcl_wc(p * q, bl, 3, 64, 30, signed=True, diag=True, force_base=2)
            assert r["factor"] in (p, q), (p, q, r)
            found += 1
            if found >= 5:
                return
    assert found > 0


def test_equal_order_split_triggered():
    """p - 1 = j l0 o, q - 1 = k l0 o with ord_p(2) = ord_q(2) = l0 o and o < m:
    gamma = 2^L has equal small order o; the split must output p via D-enumeration."""
    bl = blocks(8)
    hits = 0
    for b in bl:
        for l0 in b:
            for o in range(1, 60):
                D = l0 * o
                for j in range(2, 400, 2):
                    p = j * D + 1
                    if not gmpy2.is_prime(p) or _orders(2, p) != D:
                        continue
                    for k in range(j + 2, 2 * j, 2):
                        q = k * D + 1
                        if gmpy2.is_prime(q) and _orders(2, q) == D:
                            r = wc.dcl_wc(p * q, [b], 3, 64, 30, force_base=2)
                            assert r["factor"] in (p, q)
                            hits += bool(r["log"].get("equal_order"))
                            break
    assert hits >= 3
