import math
import random

import numpy as np
from sympy import factorint

from src.dsr import core


def _brute_cov(N, S, ws):
    G = [a for a in range(1, N) if math.gcd(a, N) == 1]
    unc = set(G)
    for m, w in enumerate(ws, 1):
        unc = {a for a in unc if not S[(a * w) % N]}
        if not unc:
            return m
    return None


def test_smooth_mask_exact():
    P = [2, 3, 5, 7]
    m = core.smooth_mask(500, P)
    for n in range(1, 500):
        assert m[n] == (max(factorint(n)) <= 7 if n > 1 else True)


def test_s_bal_unique_short_split():
    N = 211 * 223
    B, P = 13, [2, 3, 5, 7, 11, 13]
    S = core.s_bal(N, P)
    sm = core.smooth_mask(N, P)
    rng = random.Random(1)
    for _ in range(300):
        x = rng.randrange(1, N)
        if math.gcd(x, N) != 1:
            continue
        found = False
        for v in range(1, math.isqrt(N) + 1):
            u = (x * v) % N
            for uu in (u, u - N):
                if uu != 0 and uu * uu + v * v < N and math.gcd(abs(uu), v) == 1:
                    found = found or (bool(sm[abs(uu)]) and bool(sm[v]))
        assert S[x] == found


def test_coverage_matches_bruteforce():
    N = 101 * 103
    P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    S = core.s_can(N, P) & core.units_mask(N)
    G = core.units_mask(N)
    gen = core.d_box(P)
    ws = [next(gen) for _ in range(3000)]
    cov = core.coverage(N, S, G, iter(ws), cap=3000)
    assert cov["m_star"] is not None and cov["m_star"] == _brute_cov(N, S, ws)


def test_d_slide_deterministic_and_smooth():
    P = [2, 3, 5, 7, 11, 13]
    a = [w for _, w in zip(range(50), core.d_slide(P, 4))]
    b = [w for _, w in zip(range(50), core.d_slide(P, 4))]
    assert a == b
    assert all(max(factorint(w)) <= 13 and sum(factorint(w).values()) == 4 for w in a)
    assert core._beatty(10) == 16 and core._beatty(1) == 1


def test_orbit_hit_gives_dominant_row():
    # the orbit reformulation: q^(H+j) = u (mod N), u smooth < N  =>  row (H+j) e_q - e(u) is q-dominant
    N = 1009 * 1013
    B, P = core.factor_base(N, 0.8)
    S = core.s_can(N, P)
    H = 8 * math.ceil(math.log2(N))
    hits = core.orbit_first_hits(N, S, P, H, 5000)
    for q, j in hits.items():
        assert j is not None
        u = pow(q, H + j, N)
        f = factorint(u)
        assert max(f) <= B
        off = sum(e for p, e in f.items() if p != q)
        assert (H + j) - f.get(q, 0) > off


def test_dsr_output_from_hit():
    N = 1009 * 1013
    B, P = core.factor_base(N, 0.8)
    S = core.s_can(N, P)
    a = pow(3, 8 * 20, N)
    for s, w in enumerate(core.d_box(P), 1):
        if S[(a * w) % N]:
            u, v, ww = core.dsr_output_from_hit(N, a, w, P)
            assert (u - a * ww * v) % N == 0 and abs(u) <= N ** 3
            break


def test_box_obstruction_small_nonsmooth_target():
    # Prop BOX: if the least residue r of a is not B-smooth and w*r < N, then a*w mod N = w*r is not B-smooth.
    N = 1009 * 1013
    B, P = core.factor_base(N, 0.8)
    S = core.s_can(N, P)
    ell = next(x for x in range(B + 1, 10 * B) if all(x % d for d in range(2, math.isqrt(x) + 1)))
    hits = core.first_hits(N, S, [ell], core.d_box(P), 10**6)
    psi = int(core.smooth_mask((N - 1) // ell + 1, P).sum())
    assert hits[ell] is not None and hits[ell] >= psi + 1
    gen = core.d_box(P)
    for _ in range(psi):
        w = next(gen)
        assert w * ell < N and not S[(ell * w) % N]
