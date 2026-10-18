"""Tests for src/coverage38 (Theorem U cost model, Lehman dichotomy, oracle barrier)."""
import random
from math import gcd, isqrt, log, sqrt

import numpy as np
import pytest
import sympy

from src.coverage38 import theory as th
from src.coverage38 import lehman as lm
from src.coverage38 import oracle as orc


# --------------------------------------------------------------------- theory
def test_theta_exponent():
    assert th.exponent_from_theta(0.0) == pytest.approx(1 / 6)
    assert th.exponent_from_theta(0.25) == pytest.approx(0.2)
    assert th.exponent_from_theta(0.5) == pytest.approx(0.2)        # capped (3/13 > 1/5)
    from src.dchworst.theory import worst_case_exponent
    for c in [0.0, 0.3, 0.375, 0.4, 0.45, 0.5]:
        assert th.exponent_from_theta(th.theta_from_c(c)) == pytest.approx(worst_case_exponent(c))


@pytest.mark.parametrize("J,lam", [(1e3, 3.0), (1e4, 10.0), (3e4, 2.0)])
def test_single_band_optima_match_universal_formula(J, lam):
    P = 2.0 ** 70
    ref = P ** 0.4 * lam ** 0.8 * J ** -0.2
    _, cb = th.optimise(lambda Y: th.box_cost(P, Y, J, lam), J, 1e18)
    _, cl = th.optimise(lambda H: th.lehman_cost(P, J, H, lam), sqrt(lam), 1e12)
    assert cb / ref == pytest.approx(1.0, rel=1e-6)
    assert cl / ref == pytest.approx(1.0, rel=1e-6)


def test_theorem_U_lower_bound_random_designs():
    """Arbitrary heterogeneous designs: moduli at many scales, arbitrary weight sets,
    full or Lehman windows; all satisfy (WIN) with c1 = 1 and R = 1."""
    rng = random.Random(1)
    P = 1e12
    for _ in range(300):
        targets = []
        J = rng.randint(1, 30)
        for _ in range(J):
            M = rng.randint(50, 10 ** 6)
            t = rng.randint(1, 40)
            hs = th.heights_sorted_upto(8)
            rng.shuffle(hs)
            for h in hs[:t]:
                if rng.random() < 0.5:
                    pi = 1.0 / M
                else:
                    pi = rng.uniform(1e-6, 1.0) / M
                targets.append((M, h, pi, pi * pi * h * P))
        s = th.design_stats(targets, P)
        assert s["cost"] >= th.universal_lower_bound(s["lam"], s["J"], P) * (1 - 1e-12)


def test_kth_height_bound():
    hs = sorted(th.heights_sorted_upto(30))
    for k, h in enumerate(hs[:1500], start=1):
        assert h >= sqrt(k / 3.0)


# --------------------------------------------------------------------- Lehman
def _pairs(n, count, seed):
    rnd = random.Random(seed)
    out = []
    while len(out) < count:
        p = sympy.randprime(2 ** ((n - 1) // 2), 2 ** (n // 2))
        q = sympy.nextprime(p + rnd.randint(1, p - 1))
        if q <= 2 * p and 2 ** (n - 1) <= p * q < 2 ** n:
            out.append((p, q))
    return out


def test_lehman_linear_families_easy():
    # (x, x+1)-type: delta = q - p tiny at weight (1,1); (2x+1, 3x+1)-type: 3p - 2q tiny
    for p in [1000003, 1000033, 1000037]:
        q = sympy.nextprime(p)
        assert lm.lehman_factor(p * q, 1, 4) in (p, q)
    p, q = 1000099, 1500151        # 3p - 2q = -5
    assert 3 * p - 2 * q == -5
    assert lm.lehman_factor(p * q, 3, 2) in (p, q)


def test_lehman_score_implies_success():
    for (p, q) in _pairs(40, 40, 3):
        N = p * q
        X0 = 200
        if lm.lehman_score_exact_le(p, q, X0, isqrt(X0)):
            assert lm.lehman_factor(N, isqrt(X0), X0) in (p, q)
        s = lm.lehman_score_vec(np.array([p]), np.array([q]), 40)[0]
        assert (s <= X0) == lm.lehman_score_exact_le(p, q, X0, 40) or abs(s - X0) < 1e-6 * X0


def test_dirichlet_bound_all_pairs_small():
    """max L(p,q) over all of S_24 is O(P^{2/5}) (Dirichlet), and Lehman at that budget
    factors every pair (a correctness check of the dichotomy's easy side)."""
    from src.dchworst.coverage import all_pairs
    Pn, Qn = all_pairs(24)
    sc = lm.lehman_score_vec(Pn, Qn, 64)
    P = 2 ** 12
    assert sc.max() <= 4 * P ** 0.4
    X0 = int(np.ceil(sc.max()))
    for i in range(0, len(Pn), max(1, len(Pn) // 200)):
        p, q = int(Pn[i]), int(Qn[i])
        assert lm.lehman_factor(p * q, isqrt(X0), X0) in (p, q)


# --------------------------------------------------------------------- oracle
def test_farey_set_and_points():
    q = 10007
    H, k0 = 5, 30
    S = orc.farey_arc_set(q, H, k0)
    for x in range(q):
        want = any(min((B * x) % q, q - (B * x) % q) <= k0 for B in range(1, H + 1))
        assert S[x] == want


def test_reciprocal_hit_is_coverage():
    """x_l in S  =>  l | Bp - Aq with B <= H and |A| <= H*P/q + 2Y k0/q + 1."""
    p, q = 7919 * 0 + 104729, 130043
    band = list(sympy.primerange(200, 400))
    H, k0 = 6, 900
    S = orc.farey_arc_set(q, H, k0)
    xs = orc.points(p, q, band)
    for l, x in zip(band, xs):
        if S[x]:
            ok = False
            for B in range(1, H + 1):
                k = (B * x) % q
                if k > q // 2:
                    k -= q
                if abs(k) <= k0:
                    A, r = divmod(B * p - l * k, q)
                    assert r == 0 and (B * p - A * q) % l == 0
                    assert abs(A) <= H * p / q + 2 * 400 * k0 / q + 1
                    ok = True
            assert ok


def test_fourier_identity_and_prop_L():
    q = 40009
    S = orc.farey_arc_set(q, 4, 60)
    band = list(sympy.primerange(1000, 2000))
    J = len(band)
    xs = orc.points(30011, q, band)
    ind = np.zeros(q)
    ind[xs] = 1
    K = np.fft.fft(ind)
    Sh = np.fft.fft(S.astype(float))
    hits = int(S[xs].sum())
    ident = (np.conj(Sh) * K).sum().real / q
    assert abs(ident - hits) < 1e-6
    # Prop L: a random configuration in the complement has all |K(h)| <= 2 lam + 2 sqrt(J log 8q)
    rng = np.random.default_rng(0)
    lam = J * S.sum() / q
    good = 0
    for _ in range(10):
        ys = orc.random_complement(q, S, J, rng)
        assert not S[ys].any()
        ind = np.zeros(q)
        ind[ys] = 1
        mx = np.abs(np.fft.fft(ind))[1:].max()
        good += mx <= 2 * lam + orc.hoeffding_level(J, q - 1)
    assert good >= 5
