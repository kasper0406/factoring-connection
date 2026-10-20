"""Tests for src/resbarrier (agent resbarrier; notes/resbarrier_results.md)."""
from fractions import Fraction as Fr
from math import comb, gcd
import random

import sympy as sp

from src.heurnfs.models import height_lp
from src.resbarrier import forms as F
from src.resbarrier import ledger as Lg


def _gnfs(N, d):
    from src.resbarrier.e1 import gnfs_instance
    return gnfs_instance(N, d)


def test_rank2_constants():
    assert abs(float(Lg.rank2_constant(2)) - 1.5262857) < 1e-6
    assert abs(float(Lg.rank2_constant(1)) - 1.2114137) < 1e-6
    assert abs(float(Lg.rank2_constant(4)) - 1.9229994) < 1e-6


def test_alpha_floors():
    assert Lg.alpha_floor(Fr(1, 2), Fr(1, 2)) == Fr(1, 3)
    assert Lg.alpha_floor(Fr(1, 2), Fr(0)) == Fr(1, 4)
    assert Lg.alpha_floor(Fr(1), Fr(0)) == Fr(1, 2)


def test_lps():
    for d in range(2, 8):
        assert Lg.res_only_min(1, d) == Fr(1, d)
        assert Lg.gae_min(1, d) == Fr(2, d + 2) == height_lp(1, d)
        assert Lg.gae_min(d, d) == height_lp(d, d)
        assert Lg.gr_min(2, 1, d) == Fr(2, d + 2)  # G_2 = G
        # KL: rank r with a rational side costs (r-1)/d, increasing in r
        assert Lg.kl_min(2, d) == Fr(1, d)
        assert all(Lg.kl_min(r, d) == Fr(r - 1, d) for r in range(2, d + 2))
        # G_r for two nonlinear sides grows like (r-1)/(d+1)
        for r in range(2, d + 2):
            assert Lg.gr_min(r, d, d) == max(Fr(r - 1, d + 1), Fr(1, d))
    assert Lg.gr_expected_log_count(3, 10, Fr(1, 10)) < 0
    assert all(Lg.gr_expected_log_count_upper(r, d, Fr(1, d)) < 0 for d in range(3, 30) for r in range(3, d + 1))


def test_res_disc_identity_and_floor():
    rng = random.Random(1)
    for _ in range(20):
        p, q = sp.nextprime(rng.randrange(2 ** 15, 2 ** 16)), sp.nextprime(rng.randrange(2 ** 15, 2 ** 16))
        N = p * q
        d = rng.choice([2, 3, 4])
        f, g, m = _gnfs(N, d)
        rc = F.res_check(f, g, N)
        assert rc["divisible"] and rc["hadamard_ok"]
        dp = F.disc_product_floor(f, g, N)
        assert dp["identity_ok"] and dp["ge_N2"]


def test_c2_bound_small_box():
    for e, d in [(30, 2), (30, 3), (36, 4)]:
        N, f, g = F.snfs_pair(e, d)
        Fm = F.polymul(f, g)
        D = len(Fm) - 1
        vals = F.box_values(Fm, 40)
        for k in range(0, vals[-1].bit_length() + 1, 3):
            c = F.count_le(vals, 1 << k)
            assert F.c2_bound_holds(c, D, 40, 1 << k, N * N)


def test_lemma_nd_etale_gives_trivial_idempotent():
    """Lemma ND (iii)-type sanity check only: when the two linear factors take the root m at different primes
    (Res a unit mod N), gcd(f1(m), N) already splits N. The idempotent mechanism is tested below."""
    p, q = 1000003, 1000033
    N = p * q
    m = 12345
    f1 = [-(m + p), 1]              # x - (m + p): f1(m) = -p
    f2 = [-(m + q), 1]              # f2(m) = -q
    R = F.resultant(f1, f2)
    assert gcd(R, N) == 1
    assert (F.hom_eval(f1, m, 1) * F.hom_eval(f2, m, 1)) % N == 0
    assert gcd(F.hom_eval(f1, m, 1), N) == p


def test_lemma_kl_rank3_power_basis():
    """Lemma KL on rank-3 families (1, alpha, alpha^2) vs (1, m, m^2): N <= |W2(a)| <= C(d+2,2) H(W2) |a|^d."""
    a0, a1, a2, x = sp.symbols("a0 a1 a2 x")
    for (pb, qb, d) in [(14, 15, 3), (16, 16, 4), (18, 18, 3)]:
        N = sp.nextprime(3 * 2 ** (pb - 2) + 17) * sp.nextprime(5 * 2 ** (qb - 3) + 101)
        f, g, m = _gnfs(int(N), d)
        fx = sum(c * x ** i for i, c in enumerate(f))
        W2 = sp.Poly(sp.resultant(fx, a0 + a1 * x + a2 * x ** 2, x), a0, a1, a2)
        H2 = max(abs(int(c)) for c in W2.coeffs())
        v = [1, m, m * m]
        out = F.lemma_kl_check(v, lambda a: int(W2.eval({a0: a[0], a1: a[1], a2: a[2]})), int(N), d, H2)
        assert out["divisible"] and out["nonzero"] and out["ge_N"] and out["upper_ok"]


def test_mnfs_res_only_closed_form():
    c = Lg.mnfs_res_only_closed()
    assert abs(float(c) - 1.7114326376) < 1e-9
    t = sp.symbols("t", positive=True)
    expr = 32 * (sp.Rational(1, 2) + t) / (9 * t ** 2 * (2 - t) ** 2)
    assert sp.simplify(sp.diff(expr, t).subs(t, sp.sqrt(sp.Rational(2, 3)))) == 0
    assert sp.simplify(expr.subs(t, sp.sqrt(sp.Rational(2, 3))) - (36 + 16 * sp.sqrt(6)) / 15) == 0


def test_lemma_nd_idempotent_mechanism():
    """Lemma ND (ii): R = Z[x]/(f1 f2) with gcd(Res, N) = 1. The primitive idempotent e1 = u(x) f2(x)/Res
    (from u f2 + v f1 = Res) maps under x -> m to an idempotent of Z/N computed WITHOUT relations; when m is a
    root of f1 mod p and of f2 mod q (the only way both factors can map), it splits N directly."""
    x = sp.symbols("x")
    p, q = 1000003, 1000033
    N = p * q
    m = 777
    # f1 with root m mod p only, f2 with root m mod q only, coprime over Q
    f1 = sp.Poly((x - m) * (x + 3) + p, x)
    f2 = sp.Poly((x - m) * (x - 5) + q, x)
    R = int(sp.resultant(f1.as_expr(), f2.as_expr(), x))
    assert gcd(R, N) == 1
    s_, t_, h = sp.gcdex(f2.as_expr(), f1.as_expr(), x)  # s f2 + t f1 = h (h constant over Q)
    e1 = sp.Poly(sp.expand(s_ * f2.as_expr() / h), x)
    # evaluate e1(m) mod N exactly with denominators inverted mod N
    acc = 0
    for k, c in enumerate(reversed(e1.all_coeffs())):
        a, b = sp.fraction(sp.nsimplify(c))
        assert gcd(int(b), N) == 1
        acc = (acc + int(a) * pow(int(b), -1, N) * pow(m, k, N)) % N
    assert (acc * acc - acc) % N == 0          # idempotent mod N
    assert gcd(acc, N) in (p, q)               # nontrivial: splits N with no relations
