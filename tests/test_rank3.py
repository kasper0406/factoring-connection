"""Tests for agent rank3 (rank-r linked modules: Bean-Thunder volume theorem, box count, rank-r Liouville)."""
import itertools
import math
import random
from fractions import Fraction as Fr

import pytest
import sympy
from sympy.ntheory.modular import crt

from src.rank3 import ledger
from src.rank3.forms import (bean_thunder_C, box_abs_values, count_bound_terms, delta_F, discriminant, kappa,
                             liouville_check, log_C_upper_elementary, predicted_valuation, resultant, vp)

POW3 = [[1], [0, 1], [0, 0, 1]]


def _linked_pair(seed, d=3, bits=(10, 11)):
    rnd = random.Random(seed)
    p = int(sympy.nextprime(rnd.randrange(1 << (bits[0] - 1), 1 << bits[0])))
    q = int(sympy.nextprime(rnd.randrange(1 << (bits[1] - 1), 1 << bits[1])))
    N = p * q
    m = rnd.randrange(N)

    def mk():
        c = [0] * d + [1]
        for k in range(1, d):
            c[k] = rnd.randrange(-5, 6)
        v = sum(c[k] * pow(m, k, N) for k in range(1, d + 1)) % N
        c[0] = (-v) % N - (N if (-v) % N > N // 2 else 0)
        return c
    return N, (p, q), m, mk(), mk()


def test_bean_constant_r2():
    C2 = bean_thunder_C(2)
    assert abs(float(C2.mid()) - 15.8997487525690) < 1e-9  # 3 B(1/3,1/3) (Bean 1994)


def test_C_r_elementary_bound():
    for r in range(2, 9):
        assert math.log(float(bean_thunder_C(r).upper())) <= log_C_upper_elementary(r)


def test_kappa_identity():
    for r in range(2, 7):
        for D in range(r + 1, 14):
            P = predicted_valuation(r, D)
            assert Fr(P * math.factorial(D - r), math.factorial(D)) == kappa(r, D) == Fr(r * (r - 1), D * (D - 1))
            # P = r! C(D-2, r-2)
            assert P == math.factorial(r) * math.comb(D - 2, r - 2)


def test_inclusion_exclusion_monotone():
    """S(t) = sum over r-subsets I of (|I cap T| - 1)_+ satisfies S(t) >= (t-1) C(D-2, r-2) (ramified case)."""
    for D in range(4, 10):
        for r in range(2, D):
            for t in range(2, D + 1):
                T = set(range(t))
                S = sum(max(len(T & set(I)) - 1, 0) for I in itertools.combinations(range(D), r))
                assert S >= (t - 1) * math.comb(D - 2, r - 2)
                assert S == t * math.comb(D - 1, r - 1) - math.comb(D, r) + math.comb(D - t, r)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_T1_divisibility_power_basis(seed):
    N, (p, q), m, f, g = _linked_pair(seed)
    info = delta_F(f, POW3, g, POW3, 3)
    assert info["Q2"] != 0
    for pr in (p, q):
        assert vp(info["Q2"], pr) * info["exp"] >= predicted_valuation(3, 6)


def test_T1_non_power_module():
    N, (p, q), m, f, g = _linked_pair(7)
    mm = m - N if m > N // 2 else m
    P = [[1], [1, 1], [0, 2, 1]]
    Q = [[1], [1 - 2 * mm, 3], [0, 2, 1]]  # Q_j(m) = P_j(m) mod N
    ev = lambda h: sum(c * pow(m, k, N) for k, c in enumerate(h)) % N
    assert all(ev(a) == ev(b) for a, b in zip(P, Q))
    info = delta_F(f, P, g, Q, 3)
    assert info["Q2"] != 0
    for pr in (p, q):
        assert vp(info["Q2"], pr) * info["exp"] >= 24


def test_T1_ramified_prime():
    """Eisenstein f at p (totally ramified linking prime): v_p(Delta_F) >= P(r, D) still holds."""
    rnd = random.Random(11)
    p, q = 509, 1033
    N = p * q
    m = int(crt([p, q], [0, 77])[0])
    b1 = 2 * p
    b0 = (-(pow(m, 3, N) + b1 * m)) % N
    b0 -= N if b0 > N // 2 else 0
    f = [b0, b1, 0, 1]
    assert b0 % p == 0 and (b0 // p) % p != 0 and discriminant(f) % p == 0
    g0 = (-(pow(m, 3, N) + 2 * m * m + 3 * m)) % N
    g = [g0 - (N if g0 > N // 2 else 0), 3, 2, 1]
    info = delta_F(f, POW3, g, POW3, 3)
    assert vp(info["Q2"], p) * info["exp"] >= predicted_valuation(3, 6)


def test_T3_count_bound_and_T4():
    N, (p, q), m, f, g = _linked_pair(3)
    info = delta_F(f, POW3, g, POW3, 3)
    vals = box_abs_values(f, POW3, g, POW3, 3, 3)
    for k in range(0, vals[-1].bit_length() + 1, 3):
        X = 1 << k
        c = sum(1 for v in vals if v <= X)
        main, err = count_bound_terms(3, 6, 3, X, info)
        assert c <= main.upper() + err
    assert liouville_check(f, POW3, g, POW3, 3, info)["holds"]


def test_gp_fails_for_shared_side():
    """Delta_F = 0 when the two sides coincide (codex rank3_1 counterexample type: char-0 overlap)."""
    f = [-2, 0, 0, 1]
    info = delta_F(f, POW3, f, POW3, 3)
    assert info["Q2"] == 0


def test_power_basis_delta_is_disc_power():
    """Power basis: |Delta_F| = |disc(fg)|^{r(r-1)(D-2)!/(2 (D-r)!)} (Vandermonde)."""
    f, g = [-2, 0, 0, 1], [-3, 1, 0, 1]
    fg = [int(c) for c in (sympy.Poly(list(reversed(f)), sympy.Symbol('x')) *
                           sympy.Poly(list(reversed(g)), sympy.Symbol('x'))).all_coeffs()][::-1]
    info = delta_F(f, POW3, g, POW3, 3)
    e = 3 * 2 * math.factorial(4) // (2 * math.factorial(3))
    assert abs(info["Q2"]) ** info["exp"] == abs(discriminant(fg)) ** e


def test_ledger_constants():
    assert abs(float(ledger.unscoped_constant(2)) - 1.5262856567) < 1e-9
    assert abs(float(ledger.unscoped_constant(3)) - (128 / 27) ** (1 / 3)) < 1e-12
    assert abs(float(ledger.scoped_constant(3)) - 2.4228274571) < 1e-9
    for r in (2, 3, 4):
        assert abs(ledger.numeric_floor(r, False, 2000) - float(ledger.unscoped_constant(r))) < 5e-3
    # rank 2 optimal: constants increase with r
    vals = [float(ledger.unscoped_constant(r)) for r in range(2, 12)]
    assert vals == sorted(vals)
