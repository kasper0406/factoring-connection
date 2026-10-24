"""Tests for src/uprime (exact bookkeeping of Hyp U' at 3/8; toy helpers for E1/E2)."""
import random
from fractions import Fraction as Fr
from math import gcd

from src.common.instances import balanced_semiprime
from src.uprime import hybrid_bv, params, toy


def test_margins_hold_at_reference_points():
    for a, d0, kappa, b2 in [(Fr(3, 5), 1, Fr(1, 8), 2), (Fr(1, 2), 2, Fr(1, 4), 4), (Fr(1, 3), 1, Fr(1, 16), 8)]:
        bad = [n for n, ok, _ in params.margins(a, d0, kappa, b2) if not ok]
        assert bad == [], bad


def test_sq_constraint_is_the_only_constant_constraint():
    # a too small relative to d0/(kappa b2): only (SQ) fails
    bad = [n for n, ok, _ in params.margins(Fr(1, 10), 1, Fr(1, 8), 1) if not ok]
    assert bad == ["(SQ) a > 3 d0/(32 kappa b2)"]
    assert params.sq_constraint(Fr(1, 10), 1, Fr(1, 8), 8)


def test_sigma_formulas():
    assert params.sigma1(Fr(3, 5), 1, Fr(1, 8), 2) == Fr(3, 8) * (Fr(5, 3) + Fr(1, 1))
    assert params.sigma2(1, 2, Fr(1, 4), 2) == Fr(3, 8) * (Fr(1, 2) + Fr(1, 2) + 2)
    # single-y Lemma CS would put the algebraic factor at L[1/2] scale: positive coefficient of D^{1/2} lam
    assert params.sigma1_single_y(Fr(3, 5), 1, Fr(1, 8)) > 0


def test_cm_level_constraint():
    assert not params.cm_level_constraint(1, 2) and params.cm_level_constraint(1, 3)
    assert params.cm_level_constraint(2, Fr(3, 4))


def test_two_pow_3_8_bracket():
    lo, hi = params.two_pow_3_8_bracket()
    assert lo ** 8 < 8 < hi ** 8 and hi - lo == Fr(1, 10 ** 6)
    clo, chi = params.cost_constant_2L(Fr(3, 5), 1, Fr(1, 8), 2, 1)
    assert clo <= chi


def test_smooth_sieve_exact():
    v = hybrid_bv.smooth_upto(20000, 11, seg=999)
    assert list(v) == [n for n in range(1, 20001) if hybrid_bv.is_smooth(n, 11)]


def test_modulus_stats_counts():
    v = hybrid_bv.smooth_upto(5000, 7)
    st = hybrid_bv.modulus_stats(v, 11 * 6)
    assert st["phi"] == hybrid_bv.phi(66) == 20
    assert st["psi_t"] == sum(1 for n in v if gcd(int(n), 66) == 1)


def test_lattice_box_brute():
    rng = random.Random(4)
    for _ in range(6):
        q = rng.choice([101, 251, 1009, 4099])
        lam = rng.randrange(q)
        V = rng.randrange(20, 90)
        pts = toy.lattice_box(q, lam, V)
        brute = set()
        for a in range(-V, V + 1):
            for b in range(-V, V + 1):
                if (a, b) != (0, 0) and (a - lam * b) % q == 0 and gcd(a, b) == 1:
                    brute.add((a, b) if (b > 0 or (b == 0 and a > 0)) else (-a, -b))
        assert pts == sorted(brute)


def test_children_are_roots_and_child_values_divisible():
    N = balanced_semiprime(48, 4242).N
    A = toy.sample_field(N, random.Random("t"), d=3, M_bits=20)
    ctx = toy.Ctx(2 ** 11, 2 ** 16)
    seen = 0
    for ell in [p for p in range(2 ** 17 + 1, 2 ** 17 + 4000) if toy.is_prime(p)][:20]:
        r = toy.classify_level1(A, ell, 16, ctx)
        for nd in r.get("nodes", []):
            for q, rr in nd["children"]:
                for (a, b) in toy.lattice_box(q, rr, int(q ** 0.5 * 4))[:20]:
                    assert toy.F_val(A, a, b) % q == 0
                seen += 1
    assert seen > 0
