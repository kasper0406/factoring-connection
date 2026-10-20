"""Tests for src/nonsmooth (phase 16, route 3)."""
import random
from fractions import Fraction

import gmpy2
import pytest

from src.nonsmooth import census, cheng, equiv, forms

SMALL = [(5, 7), (7, 11), (11, 13), (13, 29), (17, 19), (23, 41), (29, 31)]
MID = [(101, 103), (1009, 2003), (10007, 10009), (65537, 99991)]


def test_forms_class_numbers_known():
    # h(-3)=1, h(-4)=1, h(-23)=3, h(-47)=5, h(-140)=6, h(-71)=7
    assert [forms.class_number(D) for D in (-3, -4, -23, -47, -140, -71)] == [1, 1, 3, 5, 6, 7]


def test_forms_composition_group_law():
    D = -4 * 1009 * 1013
    h = forms.class_number(D)
    rng = random.Random(3)
    fs = forms.reduced_forms(D)
    ident = forms.identity(D)
    for _ in range(10):
        f, g = rng.choice(fs), rng.choice(fs)
        assert forms.compose(f, g) == forms.compose(g, f)
        assert forms.power(f, h) == ident
        assert forms.disc(forms.compose(f, g)) == D


def test_hurwitz_small_values():
    assert forms.hurwitz_H(3) == Fraction(1, 3)
    assert forms.hurwitz_H(4) == Fraction(1, 2)
    assert forms.hurwitz_H(7) == 1
    assert forms.hurwitz_H(8) == 1
    assert forms.hurwitz_H(12) == Fraction(4, 3)


@pytest.mark.parametrize("p,q", SMALL)
def test_kronecker_hurwitz_identity_and_reduction(p, q):
    N = p * q
    S = equiv.oracle_hurwitz_sum(p, q)
    assert S == equiv.kronecker_hurwitz_rhs([1, p, q, N], N) == 2 * N + 2 * q
    assert equiv.reduce_hurwitz_sum(N, S) == p


@pytest.mark.parametrize("p,q", SMALL + MID)
def test_arithmetic_reductions(p, q):
    N = p * q
    assert equiv.reduce_phi(N, equiv.oracle_phi(p, q)) == p
    assert equiv.reduce_sigma(N, equiv.oracle_sigma(p, q)) == p
    assert equiv.reduce_psi(N, equiv.oracle_psi(p, q)) == p
    assert equiv.reduce_r4(N, equiv.oracle_r4(p, q)) == p
    assert equiv.reduce_pillai(N, equiv.oracle_pillai(p, q)) == p


def test_pillai_and_r4_brute_force():
    from math import gcd, isqrt
    for p, q in SMALL[:4]:
        N = p * q
        assert sum(gcd(k, N) for k in range(1, N + 1)) == equiv.oracle_pillai(p, q)
        r = isqrt(N)
        cnt = 0
        for a in range(-r, r + 1):
            for b in range(-r, r + 1):
                for c in range(-r, r + 1):
                    m = N - a * a - b * b - c * c
                    if m < 0:
                        continue
                    d = isqrt(m)
                    if d * d == m:
                        cnt += 1 if d == 0 else 2
        assert cnt == equiv.oracle_r4(p, q)


@pytest.mark.parametrize("p,q", [(5, 7), (7, 11), (11, 13), (13, 29), (101, 103), (1009, 2003)])
def test_genus_X0(p, q):
    N = p * q
    g = equiv.oracle_genus_X0(p, q)
    if N == 35:
        assert g == 3       # genus of X_0(35) is 3
    if N == 77:
        assert g == 7       # genus of X_0(77) is 7
    assert equiv.reduce_genus_X0(N, g) in (p, q)


def test_tau_and_bach_charles_pair():
    tau = equiv.tau_table(400)
    assert tau[1:8] == [1, -24, 252, -1472, 4830, -6048, -16744]
    # multiplicativity spot check
    assert tau[6] == tau[2] * tau[3] and tau[35] == tau[5] * tau[7]
    for p, q in [(5, 7), (11, 13), (13, 29), (17, 19), (101, 103), (199, 197)]:
        A, B = equiv.oracle_tau_pair(p, q, tau)
        if p * q <= 400:
            assert A == tau[p * q]
        f = equiv.reduce_tau_pair(p * q, A, B)
        if tau[p] % p == 0 and tau[q] % q == 0:
            assert f is None      # both non-ordinary (e.g. 5 | tau(5), 7 | tau(7)): documented
        else:
            assert f in (p, q)


@pytest.mark.parametrize("p,q", MID)
def test_quadratic_frobenius_reduction(p, q):
    N = p * q
    f, d, calls = equiv.reduce_quadratic_frobenius(
        N, lambda d: equiv.oracle_quadratic_frobenius(p, q, d))
    assert f in (p, q) and calls <= 1


@pytest.mark.parametrize("p,q", MID[:3])
def test_class_number_reduction(p, q):
    N = p * q
    h = equiv.oracle_class_number_4N(p, q)
    f, _ = equiv.reduce_class_number(N, h, random.Random(N))
    assert f in (p, q)


def test_class_number_half_sum_prime():
    # Dirichlet: h(-p) from the half-interval Legendre sum, p = 3 mod 4
    for p in (7, 11, 19, 23, 31, 43, 47, 59, 71, 103, 1019):
        assert equiv.class_number_half_sum(p) == forms.class_number(-p)


def _special_prime(D, bits, rng):
    while True:
        s = rng.getrandbits(bits // 2) | 1 | (1 << (bits // 2 - 1))
        if (D * s * s + 1) % 4:
            continue
        p = (D * s * s + 1) // 4
        if gmpy2.is_prime(p):
            return int(p)


def test_ladder_matches_affine():
    # x-only ladder against the group law, over a prime field
    p, a, b = 1000003, 5, 7
    rng = random.Random(0)

    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2 and (y1 + y2) % p == 0:
            return None
        lam = ((3 * x1 * x1 + a) * pow(2 * y1, -1, p) if P == Q
               else (y2 - y1) * pow(x2 - x1, -1, p)) % p
        x3 = (lam * lam - x1 - x2) % p
        return x3, (lam * (x1 - x3) - y1) % p

    for _ in range(20):
        while True:
            x = rng.randrange(p)
            r = (x ** 3 + a * x + b) % p
            if r and pow(r, (p - 1) // 2, p) == 1:
                break
        y = equiv._sqrt_mod_prime(r, p)
        k = rng.randrange(1, 10 ** 6)
        R, P, kk = None, (x, y), k
        while kk:
            if kk & 1:
                R = add(R, P)
            P = add(P, P)
            kk >>= 1
        X, Z = cheng.ladder_simple(x, k, a, b, p, cheng.Counters())
        if R is None:
            assert Z % p == 0
        else:
            assert Z % p and X * pow(Z, -1, p) % p == R[0]


@pytest.mark.parametrize("D", cheng.D_LIST)
def test_cheng_special_family(D):
    rng = random.Random(f"t|{D}")
    p = _special_prime(D, 96, rng)
    q = int(gmpy2.next_prime(rng.getrandbits(96) | (1 << 95)))
    N = p * q
    f, DD, c = cheng.cheng_factor(N)
    assert f in (p, q)


def test_cheng_negative_control():
    rng = random.Random("neg")
    for _ in range(3):
        p = int(gmpy2.next_prime(rng.getrandbits(64) | (1 << 63)))
        q = int(gmpy2.next_prime(rng.getrandbits(64) | (1 << 63)))
        f, _, _ = cheng.cheng_factor(p * q, tries_per_D=4, tries_D3=8)
        assert f is None


def test_census_D_is_3_mod_8():
    rows = census.sample(24, 50, 1)
    for r in rows:
        assert r["D_mod8"] == 3
        m = (4 * r["p"] - 1) // r["D"]
        assert gmpy2.is_square(m)


@pytest.mark.parametrize("ell", [5, 7, 11, 13])
def test_gl2_sign_character_is_legendre_det(ell):
    assert equiv.sign_character_check(ell)


def test_division_poly_disc_legendre():
    import sympy as sp
    P = equiv.division_poly_x(5, 2, 3)
    D = int(sp.discriminant(P))
    rng = random.Random(7)
    for _ in range(60):
        p = int(gmpy2.next_prime(rng.randrange(10 ** 5, 10 ** 6)))
        if D % p:
            assert equiv.jacobi(D % p, p) == equiv.jacobi(p, 5)


def test_frob_fixed_space_on_schoof_algebra():
    """Prop FROB(b) on Schoof's algebra (Z/N)[x]/psi_5 with (N/5) = -1."""
    import sympy as sp
    P = equiv.division_poly_x(5, 2, 3)
    lc = int(P.LC())
    rng = random.Random(11)
    done = 0
    while done < 6:
        p = int(gmpy2.next_prime(rng.randrange(10 ** 6, 10 ** 7)))
        q = int(gmpy2.next_prime(rng.randrange(10 ** 6, 10 ** 7)))
        N = p * q
        if p == q or equiv.jacobi(N, 5) != -1 or gcd_(int(sp.discriminant(P)) * lc, N) != 1:
            continue
        inv = pow(lc, -1, N)
        f = [int(c) * inv % N for c in reversed(P.all_coeffs())]   # monic mod N, low first
        M = equiv.oracle_frobenius_matrix(p, q, f)
        g, _ = equiv.reduce_frobenius_fixed_space(N, M)
        assert g in (p, q)
        done += 1


def gcd_(a, b):
    from math import gcd
    return gcd(a, b)
