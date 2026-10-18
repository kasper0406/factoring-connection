"""Tests for src/dchworst (coverage geometry; harness code)."""
import random
from math import gcd, isqrt

import numpy as np

from src.dchworst import coverage as cv
from src.dchworst import theory as th
from src.dcltheorem.params import primes_upto


def test_inverse_table():
    for l in (3, 7, 101, 8191):
        inv = cv.inverse_table(l)
        x = np.arange(1, l)
        assert np.all(x * inv[1:] % l == 1)


def test_hmin_table_bruteforce():
    for l in (101, 103, 211):
        for signed in (False, True):
            Hcap = 12
            tab = cv.hmin_table(l, Hcap, signed)
            for t in range(1, l):
                best = Hcap + 1
                for (A, B) in cv.weights(Hcap, signed):
                    if (A - t * B) % l == 0:
                        best = min(best, max(abs(A), B))
                assert tab[t] == best, (l, signed, t)


def test_minkowski_signed():
    """Signed weights: every nonzero class has mu_l(t) <= floor(sqrt(l)) + 1 (Minkowski;
    the box [-s,s] x [0,s] with s = floor(sqrt l) has (s+1)^2 > l points). t = 1 is the
    excluded weight (1,1); it is reached by e.g. (1 + l k)... so skip t = 1."""
    for l in primes_upto(3000)[10:]:
        s = isqrt(l) + 1
        tab = cv.hmin_table(l, s, True)
        bad = [t for t in range(2, l) if tab[t] > s]
        assert not bad, (l, bad[:5])


def test_positive_weights_no_minkowski():
    """Positive weights alone: t = -1 needs A + B = l, height >= (l+1)/2."""
    l = 211
    tab = cv.hmin_table(l, 60, False)
    assert tab[l - 1] == 61


def test_sweep_matches_bruteforce():
    n = 22
    P, Q = cv.all_pairs(n)
    rng = random.Random(1)
    primes = cv.band(64)[:6]
    out = cv.sweep(P, Q, primes, 12, True, record_every=[len(primes)])
    for _ in range(200):
        i = rng.randrange(len(P))
        p, q = int(P[i]), int(Q[i])
        ev = cv.covering_events(p, q, primes, 12, True)
        h = min([max(abs(A), B) for (A, B, l) in ev], default=13)
        assert out["h"][i] == h


def test_all_pairs_population():
    n = 20
    P, Q = cv.all_pairs(n)
    ps = set(primes_upto(1 << 11))
    brute = [(p, q) for p in ps for q in ps
             if p < q <= 2 * p and (1 << (n - 1)) <= p * q < (1 << n)]
    assert len(brute) == len(P)
    assert set(zip(P.tolist(), Q.tolist())) == set(brute)


def test_cost_exponent_formula():
    # c = 0 -> 1/5 ; c = 1/2 -> 1/6 ; threshold 3/8 -> 1/5
    assert abs(th.worst_case_exponent(0.0) - 0.2) < 1e-12
    assert abs(th.worst_case_exponent(0.5) - 1 / 6) < 1e-12
    assert abs(th.worst_case_exponent(0.375) - 0.2) < 1e-12
    assert th.worst_case_exponent(0.4) < 0.2
    assert abs(th.worst_case_exponent(0.25) - 0.2) < 1e-12   # c=1/4 gives no gain
    # brute-force optimisation over (theta) agrees
    for c in (0.0, 0.1, 0.3, 0.4, 0.45, 0.5):
        assert abs(th.worst_case_exponent(c) - th.worst_case_exponent_numeric(c)) < 1e-3


def test_det_form():
    """Lambda_l = {v : l | det(v, w)} = Z w + l Z^2 for w = (p, q) primitive."""
    p, q, l = 1009, 1871, 103
    for A in range(-20, 21):
        for B in range(0, 21):
            in1 = (B * p - A * q) % l == 0
            # membership in Z w + l Z^2: exists a with (A,B) == a (p,q) mod l
            in2 = any(((A - a * p) % l == 0 and (B - a * q) % l == 0) for a in range(l))
            assert in1 == in2
