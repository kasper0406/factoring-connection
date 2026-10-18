import random

from src.genus2rh.cantor import Jac, IDENTITY, NonUnit
from src.genus2rh.hecm import (exponent_schedule, monic_model, one_curve, primes_upto,
                               sample_curve, resultant_F_dF)
from src.genus2rh.diagnostic import jac_order, squarefree_mod


def _random_curve_point(p, rng):
    while True:
        f, a, b = sample_curve(p, rng)
        if f[5] % p == 0:
            continue
        F, A, B = monic_model(p, f, a, b)
        if squarefree_mod(p, F):
            return F, ((-A % p, 1), (B,) if B else ())


def test_group_law_axioms_small_p():
    rng = random.Random(1)
    for p in (31, 53, 101):
        F, _ = _random_curve_point(p, rng)
        J = Jac(p, F)
        pts = []
        while len(pts) < 3:
            x = rng.randrange(p)
            y2 = sum(F[i] * pow(x, i, p) for i in range(6)) % p
            ys = [y for y in range(p) if y * y % p == y2]
            if ys:
                pts.append(((-x % p, 1), (ys[0],) if ys[0] else ()))
        P, Q, R = pts
        assert J.compose(P, Q) == J.compose(Q, P)
        assert J.compose(J.compose(P, Q), R) == J.compose(P, J.compose(Q, R))
        assert Jac.is_identity(J.compose(P, J.neg(P)))
        assert J.compose(P, IDENTITY) == J.compose(IDENTITY, P)


def test_jacobian_order_annihilates():
    rng = random.Random(2)
    for p in (31, 61, 113, 151):
        for _ in range(3):
            F, D = _random_curve_point(p, rng)
            n = jac_order(p, F)
            assert (p ** 0.5 - 1) ** 4 <= n <= (p ** 0.5 + 1) ** 4
            assert Jac.is_identity(Jac(p, F).smul(n, D))


def test_monic_model_preserves_point():
    rng = random.Random(3)
    N = 101 * 103
    for _ in range(20):
        f, a, b = sample_curve(N, rng)
        assert (sum(f[i] * pow(a, i, N) for i in range(6)) - b * b) % N == 0
        F, A, B = monic_model(N, f, a, b)
        assert F[5] == 1
        assert (sum(F[i] * pow(A, i, N) for i in range(6)) - B * B) % N == 0


def test_schedule_and_primes():
    assert primes_upto(20) == [2, 3, 5, 7, 11, 13, 17, 19]
    s = exponent_schedule(100, 7)
    bound = 12 ** 4
    for r in (2, 3, 5, 7):
        m = s.count(r)
        assert r ** m <= bound < r ** (m + 1)
    assert s == sorted(s)


def test_resultant_detects_square_factor():
    # F = (X-1)^2 (X^3 + X + 1)
    import sympy
    X = sympy.symbols("X")
    poly = sympy.Poly(sympy.expand((X - 1) ** 2 * (X ** 3 + X + 1)), X)
    F = [int(c) for c in reversed(poly.all_coeffs())]
    assert resultant_F_dF(F) == 0
    G = [1, 0, 0, 0, 0, 1]  # X^5 + 1
    assert resultant_F_dF(G) == 3125


def test_nonunit_is_proper_divisor():
    N = 103 * 107
    rng = random.Random(4)
    found = 0
    for _ in range(300):
        d, stage = one_curve(N, 120, 40, rng)
        if d is not None:
            assert 1 < d < N and N % d == 0
            found += 1
    assert found > 0
