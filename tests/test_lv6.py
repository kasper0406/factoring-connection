"""Unit tests for src/lv6 (exact)."""
import random
from math import gcd

import sympy

from src.lv6 import numfield as nf
from src.lv6.chars import chi, sample_first_degree_prime
from src.lv6.pipeline import atom_columns, left_kernel, smooth_factor, toy_field
from src.lv6.sqrt import choose_prime, hensel_sqrt
from src.pcest.geometry import hom_eval


def _rand_f(rng, d=3, h=50):
    while True:
        f = [rng.randint(-h, h) for _ in range(d)] + [rng.randint(2, h)]
        x = sympy.Symbol("x")
        if sympy.Poly(list(reversed(f)), x).is_irreducible:
            return f


def test_monic_g_identity():
    rng = random.Random(1)
    for _ in range(20):
        f = _rand_f(rng, d=rng.choice([3, 5]))
        d = len(f) - 1
        g = nf.monic_g(f)
        for xv in range(-5, 6):
            # g(f_d x) = f_d^(d-1) f(x)
            lhs = sum(g[i] * (f[d] * xv) ** i for i in range(d + 1))
            rhs = f[d] ** (d - 1) * sum(f[i] * xv ** i for i in range(d + 1))
            assert lhs == rhs


def test_hensel_sqrt_recovers_and_rejects():
    rng = random.Random(2)
    for _ in range(15):
        f = _rand_f(rng, d=rng.choice([3, 5]))
        d = len(f) - 1
        g = nf.monic_g(f)
        lin = [nf.linear(rng.randint(-30, 30), rng.randint(1, 30), f[d], d) for _ in range(rng.randint(2, 12))]
        delta = nf.product(lin, g)
        Gamma = nf.mul(delta, delta, g)
        p = choose_prime(g, 1000, f[d] * nf.disc_poly(g))
        # p must not divide N(Gamma): move to the next admissible prime until Gamma is a unit at p
        while _norm_divisible(Gamma, g, p):
            p = choose_prime(g, p, f[d] * nf.disc_poly(g))
        A = max(abs(c) for c in delta)
        eta = hensel_sqrt(Gamma, g, A, p, rng)
        assert eta is not None and (eta == delta or eta == [-c for c in delta])
        # non-square: multiply by an element of odd valuation at a prime (a linear element not in delta)
        extra = nf.linear(97, 89, f[d], d)
        G2 = nf.mul(Gamma, extra, g)
        if not _norm_divisible(G2, g, p):
            assert hensel_sqrt(G2, g, 10 * A * max(abs(c) for c in extra), p, rng) is None


def _norm_divisible(x, g, p):
    import flint
    ctx = flint.fmpz_mod_poly_ctx(p)
    return ctx([c % p for c in g]).gcd(ctx([c % p for c in x])).degree() > 0


def test_coeff_bound_is_valid():
    rng = random.Random(3)
    for _ in range(20):
        f = _rand_f(rng, d=rng.choice([3, 5]))
        d = len(f) - 1
        g = nf.monic_g(f)
        lin = [nf.linear(rng.randint(-40, 40), rng.randint(1, 40), f[d], d) for _ in range(2 * rng.randint(1, 6))]
        delta = nf.product(lin, g)                       # sqrt of prod(lin)^2 ... use eta = g' * delta
        eta = nf.mul(nf.gprime_elem(g), delta, g)
        A = nf.coeff_bound_sqrt(g, [nf.norm1(x) for x in lin + lin])
        assert max(abs(c) for c in eta) <= A


def test_psi_is_ring_hom():
    rng = random.Random(4)
    N = 1000003 * 999983
    ev = []
    fld = toy_field(N, 3, 4, rng, ev)
    f = fld["f"]
    g = nf.monic_g(f)
    m = fld["g0"] * pow(fld["g1"], -1, N) % N
    bimg = f[3] * m % N
    assert sum(g[i] * pow(bimg, i, N) for i in range(4)) % N == 0
    for _ in range(20):
        x = [rng.randint(-10 ** 6, 10 ** 6) for _ in range(3)]
        y = [rng.randint(-10 ** 6, 10 ** 6) for _ in range(3)]
        assert nf.psi(nf.mul(x, y, g), bimg, N) == nf.psi(x, bimg, N) * nf.psi(y, bimg, N) % N
    # psi(f_d a - b beta) = f_d * G(a, b) / g1
    for a, b in [(3, 7), (-5, 2), (11, 13)]:
        G = fld["g1"] * a - fld["g0"] * b
        assert nf.psi(nf.linear(a, b, f[3], 3), bimg, N) == f[3] * G * pow(fld["g1"], -1, N) % N


def test_good_prime_columns_lemma_C():
    """For r not dividing D = disc(f) f_d and coprime (a, b): r | F(a, b) implies r does not divide b,
    and the single first-degree prime (r, a/b) carries ord_r F(a, b)."""
    rng = random.Random(5)
    f = _rand_f(rng, d=3, h=200)
    D = nf.disc_poly(f) * f[3]
    for _ in range(400):
        a, b = rng.randint(-300, 300), rng.randint(1, 300)
        if gcd(a, b) != 1:
            continue
        F = hom_eval(f, a, b)
        if F == 0:
            continue
        for r in sympy.factorint(abs(F)):
            if D % r:
                assert b % r
                s = a * pow(b, -1, r) % r
                assert sum(f[i] * s ** i for i in range(4)) % r == 0


def test_left_kernel():
    rng = random.Random(6)
    rows = [rng.getrandbits(30) for _ in range(45)]
    ker = left_kernel(rows, 30)
    assert len(ker) >= 15
    for S in ker:
        acc = 0
        for i in range(len(rows)):
            if (S >> i) & 1:
                acc ^= rows[i]
        assert acc == 0 and S != 0


def test_smooth_factor_and_chars():
    P = 1
    for p in sympy.primerange(2, 50):
        P *= p
    assert smooth_factor(2 ** 5 * 3 * 47, P) == {2: 5, 3: 1, 47: 1}
    assert smooth_factor(2 * 53, P) is None
    rng = random.Random(7)
    f = _rand_f(rng, d=3)
    r, s, _ = sample_first_degree_prime(f, 64, rng)
    assert sympy.isprime(r) and sum(f[i] * s ** i for i in range(4)) % r == 0
    assert chi(r, s, 3, 5) in (0, 1)
