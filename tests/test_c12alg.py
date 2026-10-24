"""Tests for c12alg (notes/c12alg_results.md).  Exact arithmetic."""
from __future__ import annotations

import random
from fractions import Fraction as Fr
from math import gcd, prod

import pytest

from src.c12alg import dw, params
from src.c12alg import analyze
from src.thmA.numfield import hom_eval
from src.thmA import numfield as nf


def _prim_pair(rng, V, allow_b0=False):
    while True:
        P1 = (rng.randint(-V, V), rng.randint(0 if allow_b0 else 1, V))
        P2 = (rng.randint(-V, V), rng.randint(0 if allow_b0 else 1, V))
        if gcd(*P1) == 1 and gcd(*P2) == 1 and dw.det2(P1, P2) != 0:
            return P1, P2


def test_short_preimages_and_representatives():
    rng = random.Random(1)
    for _ in range(120):
        d = rng.choice([3, 4, 5])
        P1, P2 = _prim_pair(rng, 4)
        R = dw.representatives(P1, P2, d)
        v1, v2 = dw.v_of(P1, d), dw.v_of(P2, d)
        det = v1[d - 2] * v2[d - 1] - v1[d - 1] * v2[d - 2]
        assert len(R) * abs(dw.det2(P1, P2)) == abs(det)          # [Gamma : Gamma_jk] = |det M| / Delta
        V = max(map(abs, P1 + P2))
        for z in R[:30]:
            x = dw.short_preimage(P1, P2, d, z)
            assert max(map(abs, x)) <= 1 + 4 * V * V


def test_tv_lemma_exact_small():
    rng = random.Random(2)
    from src.c12alg.e1 import _tv_int
    for _ in range(6):
        P1, P2 = _prim_pair(rng, 2)
        if P1[1] * P2[1] < 2:
            continue
        d, H = 3, 30
        k1, c1 = dw.law_phi(P1, P2, d, H)
        R = dw.representatives(P1, P2, d)
        k2, c2 = dw.shift_mix(k1, c1, R)
        n1 = (2 * H + 1) ** d
        tv = _tv_int(k1, c1, n1, k2, c2, n1 * len(R))
        assert tv <= dw.tv_lemma_bound(P1, P2, d, H)


def test_unit_residue_lemma_CU():
    """Fhat_N(P) is a unit mod t = |G(P)| whenever gcd(t, N) = 1, including gcd(t, g1) > 1."""
    rng = random.Random(3)
    seen_g1_common = 0
    for _ in range(4000):
        N = rng.randrange(10 ** 5, 10 ** 8) | 1
        g0, g1 = rng.randint(6, 40), rng.randint(6, 40)
        if gcd(g0, g1) != 1 or gcd(g1, N) != 1:
            continue
        P = (rng.randint(-60, 60), rng.randint(1, 60))
        r = dw.fhat_unit_mod_t(N, rng.choice([3, 5, 7]), (g0, g1), P)
        if r is None:
            continue
        assert r
        if gcd(abs(dw.G_of((g0, g1), P)), g1) > 1:
            seen_g1_common += 1
    assert seen_g1_common > 20


def test_dw_lattice_DWL():
    """(F(P1), F(P2)) over (w, c) satisfies F1 = lam^d F2 mod Delta and hits every class of L / Delta Z^2."""
    rng = random.Random(4)
    for _ in range(25):
        d = 3
        N = rng.randrange(10 ** 3, 10 ** 4) | 1
        g = (rng.randint(5, 12), rng.randint(5, 12))
        if gcd(*g) != 1 or gcd(g[1], N) != 1:
            continue
        P1, P2 = _prim_pair(rng, 3)
        G1, G2 = dw.G_of(g, P1), dw.G_of(g, P2)
        if G1 == 0 or G2 == 0 or gcd(G1 * G2, N) != 1:
            continue
        lam, D = dw.lam_mod(P1, P2)
        seen = set()
        for w in range(1, abs(G1 * G2) * D + 1):
            c = [rng.randint(-50, 50) for _ in range(d)]
            F1 = dw.dw_value(N, w, d, g, c, P1)
            F2 = dw.dw_value(N, w, d, g, c, P2)
            assert (F1 - pow(lam, d, D) * F2) % D == 0
            seen.add((F1 % D, F2 % D))
        # L mod D = {(lam^d x, x)}: D classes, all hit (w ranges over a full period, c random)
        assert len(seen) == D or D == 1


def test_content_lemma_CT_and_IR_DW():
    rng = random.Random(5)
    n_cont = 0
    for _ in range(3000):
        N = rng.randrange(10 ** 5, 10 ** 6) | 1
        d = 3
        g0, g1 = rng.randint(20, 40), rng.randint(20, 40)
        if gcd(g0, g1) != 1 or gcd(g1, N) != 1:
            continue
        w = rng.randint(2, 2000)
        c = [rng.randint(-3, 3) for _ in range(d)]
        f = dw.dw_coeffs(N, w, d, g0, g1, c)
        e = 0
        for x in f:
            e = gcd(e, x)
        e1 = e // gcd(e, N)                      # part of the content coprime to N... (N-part is a factor exit)
        if gcd(e, N) == 1 and e > 1:
            n_cont += 1
            assert w % e == 0
            assert all(ci % e == 0 for ci in c) or any(all(ci % p == 0 for ci in c) for p in range(2, e + 1)
                                                       if e % p == 0)
        # IR-DW: r | w, r not dividing g1 => f mod r vanishes at m_r = g0/g1
        for r in (5, 7, 11, 13):
            if w % r == 0 and g1 % r:
                m_r = g0 * pow(g1, -1, r) % r
                assert sum(fi * pow(m_r, i, r) for i, fi in enumerate(f)) % r == 0
    assert n_cont > 10


def test_root_event_uniform_in_g():
    """Pr_{w,c}[q | F(r, 1) | g] = 1/q exactly over a full period (both cases q | G(r,1) and not)."""
    rng = random.Random(6)
    for _ in range(20):
        q = rng.choice([7, 11, 13])
        N = rng.randrange(10 ** 3, 10 ** 4)
        if N % q == 0:
            continue
        g = (rng.randint(5, 30), rng.randint(5, 30))
        if gcd(*g) != 1 or gcd(g[1], N) != 1:
            continue
        r = rng.randrange(q)
        if rng.random() < 0.3:                      # force q | G(r, 1): r = g0 / g1 mod q
            if g[1] % q == 0:
                continue
            r = g[0] * pow(g[1], -1, q) % q
        hits = tot = 0
        for w in range(q):
            for cl in range(q):                     # c_{d-1} residue (the others fixed)
                c = [3, -2, cl]
                Fv = dw.dw_value(N, w, 3, g, c, (r, 1))
                hits += (Fv % q == 0)
                tot += 1
        assert hits * q == tot


def test_max_minor_lemma_MM():
    rng = random.Random(7)
    cases = []
    for _ in range(3000):
        cases.append(_prim_pair(rng, 30, allow_b0=True))
    cases += [((25, 1), (1, 29)), ((25, 1), (-1, 29)), ((7, 1), (-29, 4)), ((1, 0), (0, 1)), ((5, 3), (-5, 3))]
    for d in (3, 4, 5, 7, 9):
        for P1, P2 in cases:
            if dw.det2(P1, P2) == 0:
                continue
            v1, v2 = dw.v_of(P1, d), dw.v_of(P2, d)
            mx = max(abs(v1[j] * v2[k] - v1[k] * v2[j]) for j in range(d) for k in range(j + 1, d))
            A = max(abs(P1[0]), abs(P1[1])) * max(abs(P2[0]), abs(P2[1]))
            assert 4 * mx >= abs(dw.det2(P1, P2)) * A ** (d - 2)


def test_bezoutian_bound():
    rng = random.Random(8)
    for _ in range(500):
        n = rng.choice([3, 4, 5, 6])
        A = [rng.randint(-50, 50) for _ in range(n + 1)]
        K = [rng.randint(-50, 50) for _ in range(n + 1)]
        P1, P2 = _prim_pair(rng, 9, allow_b0=True)
        ev = lambda f, P: sum(f[i] * P[0] ** (n - i) * P[1] ** i for i in range(n + 1))
        E = ev(A, P1) * ev(K, P2) - ev(A, P2) * ev(K, P1)
        Dl = dw.det2(P1, P2)
        assert E % Dl == 0
        V1 = max(map(abs, P1)); V2 = max(map(abs, P2))
        assert abs(E) <= 2 * n ** 3 * max(map(abs, A)) * max(map(abs, K)) * abs(Dl) * (V1 * V2) ** (n - 1)


def test_pair_sum_identity_PS():
    """C = M + sum_a E1(kappa a) E2(a): cross terms vanish exactly (gcd-class averages, kappa a unit)."""
    rng = random.Random(9)
    for _ in range(40):
        Dl = rng.randint(2, 40)
        kap = rng.randint(1, Dl - 1)
        if gcd(kap, Dl) != 1:
            continue
        A1 = [n for n in range(1, 400) if rng.random() < 0.3]
        A2 = [n for n in range(1, 400) if rng.random() < 0.3]
        C = sum(1 for n1 in A1 for n2 in A2 if (n1 - kap * n2) % Dl == 0)
        cnt1 = [0] * Dl; cnt2 = [0] * Dl
        for n in A1: cnt1[n % Dl] += 1
        for n in A2: cnt2[n % Dl] += 1
        divs = [e for e in range(1, Dl + 1) if Dl % e == 0]
        cls = {e: [a for a in range(Dl) if gcd(a, Dl) == e] for e in divs}
        m1 = {e: Fr(sum(cnt1[a] for a in cls[e]), len(cls[e])) for e in divs}
        m2 = {e: Fr(sum(cnt2[a] for a in cls[e]), len(cls[e])) for e in divs}
        M = sum(len(cls[e]) * m1[e] * m2[e] for e in divs)
        EE = sum((cnt1[kap * a % Dl] - m1[gcd(a, Dl)]) * (cnt2[a] - m2[gcd(a, Dl)]) for a in range(Dl))
        assert C == M + EE


def test_local_factor_product_formula():
    beta = Fr(4, 5)
    for Dl in (6, 30, 210, 15, 77):
        ps = [p for p in range(2, Dl + 1) if Dl % p == 0 and all(p % r for r in range(2, p))]
        # multiplicative model with tau_p = 1 + eps_p (rational stand-in for p^(1-beta))
        eps = {p: Fr(1, p + 3) for p in ps}
        tot = Fr(0)
        for mask in range(1 << len(ps)):
            pi1 = Fr(1); pi0 = Fr(1)
            for i, p in enumerate(ps):
                tp = 1 + eps[p]
                if mask >> i & 1:
                    pi1 *= tp / p; pi0 *= Fr(1, p)
                else:
                    pi1 *= (1 - tp / p); pi0 *= (1 - Fr(1, p))
            tot += pi1 * pi1 / pi0
        assert tot == prod(1 + eps[p] ** 2 / (p - 1) for p in ps)


def test_content_probability_formula():
    p, _ = analyze.p_content(1000, 4, 3)
    # brute force over w in (1000, 2000], c in [-4, 4]^3
    import itertools
    cnt = tot = 0
    for w in range(1001, 2001, 7):
        for c in itertools.product(range(-4, 5), repeat=3):
            g = w
            for ci in c:
                g = gcd(g, ci)
            cnt += g > 1
            tot += 1
    assert abs(cnt / tot - p) < 0.01


def test_params_window():
    for al in [Fr(3, 8) + Fr(k, 400) for k in range(1, 10)]:
        assert params.all_ok(al), (al, params.margins(al))
    m = params.margins(Fr(3, 8))
    assert m["first-order child error D^-alpha << eta = D^-(3-7alpha)"] != "ok"
    assert not params.all_ok(Fr(2, 5))


def test_dw_pipeline_smoke():
    import json
    from src.c12alg import dw_pipeline as dp
    cfg = json.load(open("experiments/configs/c12alg_e2.json"))["pipeline"]
    N = 1000003 * 999983
    r = dp.factor_dw(N, cfg, 12345)
    assert r["factors"] == [999983, 1000003]
