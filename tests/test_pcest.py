"""Tests for src/pcest (exact lattice facts of notes/pcest_results.md)."""
import random
from fractions import Fraction
from math import gcd

import pytest

from src.pcest import geometry as gm
from src.pcest import params as pr


def _rand_coprime(rng, lo, hi):
    while True:
        g0, g1 = rng.randrange(lo, hi), rng.randrange(lo, hi)
        if gcd(g0, g1) == 1:
            return g0, g1


@pytest.mark.parametrize("d", [3, 5, 7])
def test_gbe_sizes_and_identity(d):
    rng = random.Random(d)
    for _ in range(30):
        N = rng.randrange(2 ** 80, 2 ** 81) | 1
        M = int(round(N ** (1.0 / (d + 1))))
        g0, g1 = _rand_coprime(rng, M // 2, M)
        a = gm.gbe(N, d, g0, g1)
        assert gm.hom_eval(a, g0, g1) == N
        assert all(abs(x) <= g1 // 2 + 1 for x in a[1:])
        # Lemma GB: |a0| <= N/g1^d + sum_i (g1/2)(g0/g1)^i <= (1 + d 2^(d-1)) M-ish
        assert abs(a[0]) <= 2 * N // g1 ** d + d * g1 * 2 ** d


@pytest.mark.parametrize("d", [3, 4, 5])
def test_field_congruence(d):
    """F(P) g1^d = N b^d (mod G(P)) for the randomized field."""
    rng = random.Random(100 + d)
    for _ in range(40):
        N = rng.randrange(2 ** 60, 2 ** 61) | 1
        M = int(round(N ** (1.0 / (d + 1))))
        g0, g1 = _rand_coprime(rng, M // 2, M)
        R = [rng.randrange(-M, M) for _ in range(d)]
        f = gm.field(N, d, g0, g1, R)
        for _ in range(10):
            a, b = rng.randrange(-500, 500), rng.randrange(1, 500)
            G = gm.G_value(g0, g1, a, b)
            if G == 0:
                continue
            assert (gm.hom_eval(f, a, b) * g1 ** d - N * b ** d) % G == 0


def _pairs_in_special_q(rng, l, g0, g1, V, count):
    """Random primitive pairs P, P' in L_l(g) (a g1 = b g0 mod l), |coords| <= V."""
    t = g0 * pow(g1, -1, l) % l          # a = t b (mod l)
    out = []
    while len(out) < count:
        pts = []
        while len(pts) < 2:
            b = rng.randrange(1, V)
            a = (t * b) % l + l * rng.randrange(-2, 2)
            if abs(a) <= V and gcd(a, b) == 1 and b % l:
                pts.append((a, b))
        P, Pp = pts
        if P[0] * Pp[1] - Pp[0] * P[1] != 0:
            out.append((P, Pp))
    return out


def test_rational_pair_lattice_graph():
    """Lemma Q: image lattice has index q = |Delta|/l and is a graph n' = lam n (mod q)."""
    rng = random.Random(7)
    for l in [101, 211, 1009]:
        g0, g1 = 777, 1234
        for P, Pp in _pairs_in_special_q(rng, l, g0, g1, 3 * l, 25):
            gens, q, lam = gm.rational_pair_lattice(l, P, Pp)
            Delta = P[0] * Pp[1] - Pp[0] * P[1]
            assert q == abs(Delta) // l
            if q > 1:
                assert gcd(lam, q) == 1
                for n, n2 in gens:
                    assert gm.in_gamma_q(n, n2, q, lam)
                # lam = b'/b = kappa mod q (P' = kappa P mod Delta)
                (a, b), (ap, bp) = P, Pp
                if gcd(b, q) == 1:
                    assert (lam - bp * pow(b, -1, q)) % q == 0


def test_rational_pair_lattice_bruteforce():
    """Direct enumeration: the image of L_l(P) contains exactly the Gamma_q points of a box."""
    l = 13
    P, Pp = (5, 1), (-29, 2)                # Delta = 5*2 + 29 = 39 = 13*3
    Delta = P[0] * Pp[1] - Pp[0] * P[1]
    assert Delta % l == 0
    gens, q, lam = gm.rational_pair_lattice(l, P, Pp)
    img = set()
    for g0 in range(-60, 61):
        for g1 in range(-60, 61):
            n, n2 = gm.G_value(g0, g1, *P), gm.G_value(g0, g1, *Pp)
            if n % l == 0 and n2 % l == 0:
                img.add((n // l, n2 // l))
    box = [(n, n2) for n in range(-5, 6) for n2 in range(-5, 6)]
    for z in box:
        assert (z in img) == gm.in_gamma_q(z[0], z[1], q, lam)


@pytest.mark.parametrize("d", [3, 4, 5, 6])
def test_alg_index_and_kernel(d):
    """Prop G / (I6): index of Gamma is |Delta| for primitive pairs; Lemma Kern: the
    coefficient kernel is Q_PP' * Z[x,y]_{d-3} (saturated, rank d-2, inside the kernel)."""
    rng = random.Random(1000 + d)
    for _ in range(25):
        while True:
            P = (rng.randrange(-60, 60), rng.randrange(1, 60))
            Pp = (rng.randrange(-60, 60), rng.randrange(1, 60))
            Delta = P[0] * Pp[1] - Pp[0] * P[1]
            if gcd(*P) == 1 and gcd(*Pp) == 1 and Delta != 0:
                break
        assert gm.alg_pair_index(P, Pp, d) == abs(Delta)
        K = gm.kernel_basis_Q(P, Pp, d)
        v, w = gm.veronese(P, d - 1), gm.veronese(Pp, d - 1)
        for c in K:
            assert sum(ci * vi for ci, vi in zip(c, v)) == 0
            assert sum(ci * wi for ci, wi in zip(c, w)) == 0
        if d >= 3:
            assert gm.is_saturated(K)


def test_divisor_window():
    """Every y-smooth n >= D has a divisor e in [D, D y] (used in Lemma DS)."""
    from src.nfsfruit.smooth import FactorBase
    y, D = 50, 300
    fb = FactorBase(y)
    for n in range(D, 40000):
        fac = fb.factor(n)
        if fac is None:
            continue
        ps = sorted(p for p, e in fac.items() for _ in range(e))
        e = 1
        for p in ps:
            e *= p
            if e >= D:
                break
        assert D <= e <= D * y and n % e == 0


def test_params():
    # d = 9, k = 1000: h = max(1/2000, 1/10 - 7/4000) = 1/10 - 7/4000; theta = 3/10 + 1/4000
    assert pr.theta(9, 1000, "pair") == Fraction(3, 10) + Fraction(1, 4000)
    assert abs(pr.c_limit(9) - 0.91924) < 1e-4
    assert abs(pr.c_limit(7) - 0.97227) < 1e-4
    assert pr.c_limit(5) > 1
    # floors bind only at small k; the infimum is unchanged by the pair floor
    for d in (7, 9, 11, 19):
        k = int(pr.binding_floor_k(d, "pair")) + 1
        assert pr.theta(d, k, "pair") == pr.theta(d, k, "lv")
