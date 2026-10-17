"""Tests for Track B/D (agent sketchBD): proved identities checked exactly on small cases."""
import random
from math import comb, gcd

import flint
import pytest

from src.common.instances import balanced_semiprime
from src.common.ring import PolyQuotient, Zmod, binomial_quotient, power_of_linear
from src.trackBD import diagnose as D
from src.trackBD.linalg import bm_field, bm_mod_N, elim_profile, mult_matrix
from src.trackBD.sketch import (binomial_modulus, euclid_readout, lifted_derivative_divided,
                                power_sketch, scalar_F, scalar_G, sketch, sketch_G)


def _Fm_mod(l, f, m):
    fl = flint.nmod_poly([x % l for x in f], l)
    one = flint.nmod_poly([1], l)
    X = flint.nmod_poly([0, 1], l)
    return fl, ((X + one).pow_mod(m, fl) - one - X.pow_mod(m, fl)) % fl


def test_sketch_matches_pure_python():
    N = 3 * 5 * 7 * 11 * 13 + 2 * 0  # 15015 composite
    for r, c in [(1, 2), (3, 5), (7, 1), (12, 3)]:
        R = binomial_quotient(Zmod(N), r, c)
        ref = R.sub(R.sub(power_of_linear(R, N), R.one), R.pow(R.gen(), N))
        assert list(ref) == sketch(N, binomial_modulus(r, c))


@pytest.mark.parametrize("seed", range(12))
def test_theorem_B1_frobenius_identity(seed):
    """Over F_p: F_N == F_m(X)^(p^k) (mod f), N = p^k m.  Includes prime powers and 3 primes."""
    rng = random.Random(seed)
    Ns = [balanced_semiprime(24, seed).N, 3 ** 4 * 7, 5 * 7 * 11 * 13, 49 * 11, 7 ** 5]
    for N in Ns:
        fs = [binomial_modulus(rng.randint(1, 9), rng.randint(1, 30)),
              [int(x) for x in flint.fmpz_poly.cyclotomic(rng.randint(3, 20)).coeffs()],
              [rng.randint(-4, 4) for _ in range(rng.randint(1, 6))] + [1]]
        for f in fs:
            S = sketch(N, f)
            fac = []
            M, d = N, 2
            while d * d <= M:
                while M % d == 0:
                    fac.append(d)
                    M //= d
                d += 1
            if M > 1:
                fac.append(M)
            for p, k, m in D.prime_data(N, fac):
                fl, Fm = _Fm_mod(p, f, m)
                rhs = Fm.pow_mod(p ** k, fl) if not Fm.is_zero() else Fm
                assert flint.nmod_poly([x % p for x in S], p) == rhs % fl


def test_resultant_corollary():
    """Res(f, S) == Res(f, F_q) (mod p) for N = p q."""
    for seed in range(10):
        inst = balanced_semiprime(28, seed)
        N = inst.N
        p, q = inst.labels["factors"]
        f = binomial_modulus(5, 7 + seed)
        S = sketch(N, f)
        eu = euclid_readout(f, S, N)
        if eu["res"] is None:
            assert 1 < eu["factor"] < N
            continue
        for a, b in ((p, q), (q, p)):
            fl, Fb = _Fm_mod(a, f, b)
            assert eu["res"] % a == int(fl.resultant(Fb))


def test_scalar_G_root_count_and_public_root():
    """G_1 has exactly gcd(p-1,q-1)-1 roots mod p; c = -1/2 is always one (public)."""
    for seed in range(6):
        inst = balanced_semiprime(20, seed)
        N = inst.N
        p, q = inst.labels["factors"]
        roots = [c for c in range(p) if scalar_G(N, c, 1) % p == 0]
        assert len(roots) == gcd(p - 1, q - 1) - 1
        assert (-pow(2, -1, p)) % p in roots


def test_anharmonic_symmetry_of_roots():
    """Roots of (1+x)^e = 1 + x^e (e odd) are stable under x -> 1/x and x -> -1-x."""
    for p, e in [(1009, 1481), (199, 241), (4093, 3001), (811, 1201)]:
        Z = {x for x in range(p) if (pow(1 + x, e, p) - 1 - pow(x, e, p)) % p == 0}
        for x in Z:
            if x:
                assert pow(x, -1, p) in Z
            assert (-1 - x) % p in Z


def test_lift_derivative_equals_G():
    for seed in range(5):
        inst = balanced_semiprime(24, seed)
        N = inst.N
        for r, c in [(3, 2), (8, 5), (1, 3)]:
            out = lifted_derivative_divided(N, r, c)
            assert out["G_from_lift"] == sketch_G(N, binomial_modulus(r, c), 1)


def test_quotient_identity_dX_S_equals_minus_r_X_dc_S():
    """d/dX S = - r X^(r-1) d/dc S  (mod N, X^r - c): the quotient fold is determined by S."""
    inst = balanced_semiprime(24, 3)
    N = inst.N
    r, c = 5, 7
    # d/dc S via dual numbers: base ring R_N[eps]/(eps^2), modulus X^r - (c + eps)
    base = PolyQuotient(Zmod(N), [0, 0, 1], name="e")  # eps^2 = 0
    ce = base.from_coeffs([c, 1])
    R = binomial_quotient(base, r, ce)
    S = R.sub(R.sub(power_of_linear(R, N), R.one), R.pow(R.gen(), N))
    S0 = [a[0] for a in S]
    dS_dc = [a[1] for a in S]
    assert S0 == sketch(N, binomial_modulus(r, c))
    dX = [(j * S0[j]) % N for j in range(1, r)] + [0]
    # r X^(r-1) * dS_dc mod (X^r - c)
    prod = [0] * (2 * r)
    for j, v in enumerate(dS_dc):
        prod[j + r - 1] += r * v
    red = [0] * r
    for i in range(2 * r - 1, -1, -1):
        if i >= r:
            prod[i - r] += c * prod[i]
        else:
            red[i] = prod[i] % N
    assert [(a + b) % N for a, b in zip(dX, red)] == [0] * r


def test_hasse_divided_is_unit_times_G():
    """D^[k] F_N / N == binom(N-1,k-1)/k * G_k (mod N) for gcd(k, N) = 1 (integer identity)."""
    N = 35 * 1  # small, explicit
    for k in (1, 2, 3, 4):
        # integer polynomial check on coefficients
        Hk = [comb(j, k) * (comb(N, j) - (1 if j == N else 0) - (1 if j == 0 else 0)) for j in range(N + 1)]
        Hk = [Hk[j + k] for j in range(N + 1 - k)]  # coefficients of D^[k] F_N
        assert all(x % N == 0 for x in Hk)
        coef = comb(N - 1, k - 1) * pow(k, -1, N) % N
        G = [(comb(N - k, j) - (1 if j == N - k else 0)) for j in range(N - k + 1)]
        assert [(x // N) % N for x in Hk] == [(coef * g) % N for g in G]


def test_euclid_and_elim_find_rank_split():
    """N = 199*241, f = X^2 - c with an F_199 root of F_241: rank 1 mod 199, 2 mod 241."""
    p, q = 199, 241
    N = p * q
    roots = [x for x in range(2, p - 1) if (pow(1 + x, q, p) - 1 - pow(x, q, p)) % p == 0
             and (x * x + x + 1) % p]
    assert len(roots) == 6
    x = roots[0]
    # c with c == x^2 mod p and c random-ish mod q (CRT)
    cp, cq = pow(x, 2, p), 5
    c = (cp * q * pow(q, -1, p) + cq * p * pow(p, -1, q)) % N
    f = binomial_modulus(2, c)
    S = sketch(N, f)
    eu = euclid_readout(f, S, N)
    el = elim_profile(mult_matrix(N, f, S), N)
    assert eu["factor"] in (p, q)
    assert el["factor"] in (p, q)


def test_bm_linear_complexity_over_c():
    """c -> F_N(c) mod p is a polynomial of degree e_p - 1, e_p = q mod (p-1): LC_p = e_p."""
    for seed in range(4):
        inst = balanced_semiprime(16, seed)
        N = inst.N
        p, q = inst.labels["factors"]
        s = [scalar_F(N, c) for c in range(3 * p)]
        assert bm_field(s, p)["L"] == q % (p - 1)
        assert bm_field(s, q)["L"] == p
        assert bm_mod_N(s, N)["factor"] in (p, q)


def test_stirling_window():
    """Delta^k F_N(0) = (k+1)! S(N, k+1) == 0 mod p for e_p <= k <= p-2."""
    inst = balanced_semiprime(16, 2)
    N = inst.N
    p, q = inst.labels["factors"]
    e = q % (p - 1)
    F = [scalar_F(N, c) for c in range(p)]
    for k in range(e, p - 1):
        d = sum((-1) ** (k - i) * comb(k, i) * F[i] for i in range(k + 1)) % N
        assert d % p == 0
    k = e - 1
    d = sum((-1) ** (k - i) * comb(k, i) * F[i] for i in range(k + 1)) % N
    assert d % p != 0


def test_prime_power_always_vanishes_mod_p():
    for N, p in [(3 ** 7, 3), (11 ** 3, 11), (101 ** 2, 101)]:
        S = sketch(N, binomial_modulus(4, 3))
        assert all(x % p == 0 for x in S)
