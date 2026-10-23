"""Tests for src/multidescent (size model, toy children / fiber algebra)."""
import random
from fractions import Fraction as Fr
from math import gcd

from src.common.instances import balanced_semiprime
from src.multidescent import model, toy


def test_alpha_j_closed_form_matches_grid():
    for j in range(1, 7):
        a = model.alpha_j(j)
        assert model.chain_feasible(a, j) is not None
        assert model.chain_feasible(a - Fr(1, 10 ** 6), j) is None
        g = model.min_alpha_grid(j, 1320)
        assert a <= g <= a + Fr(1, 1320)


def test_two_level_parameters():
    p = model.two_level_parameters()
    assert p["alpha"] == Fr(3, 8) and p["delta"] == Fr(1, 4) and p["zeta1"] == Fr(1, 2)
    assert p["level1_alg_value"] == Fr(7, 8) and p["level2_value"] == Fr(3, 4)
    assert p["children"] == Fr(3, 8)          # d^{3/2} children, d = D^{1/4}
    assert [model.alpha_j(j) for j in (1, 2, 3, 4)] == [Fr(2, 5), Fr(3, 8), Fr(4, 11), Fr(5, 14)]


def test_union_bound_margins():
    assert model.union_bound_margin(2) == [Fr(3, 8)]
    assert model.union_bound_margin(3) == [Fr(4, 11), Fr(0)]
    assert model.union_bound_margin(4)[-1] < 0


def _field(d=3, M_bits=20, seed=5):
    N = balanced_semiprime(48, 777).N
    return N, toy.sample_field(N, random.Random(seed), d=d, M_bits=M_bits)


def test_smooth_parts_exact():
    P = toy.primorial(toy.primes_upto(50))
    xs = [2 ** 5 * 3 * 47, 53 * 2, 1, 97 * 89, 2 * 3 * 5 * 7 * 11 * 13]
    sp = toy.smooth_parts(xs, P)
    assert [sp[i] == xs[i] for i in range(len(xs))] == [True, False, True, False, True]
    assert sp[1] == 2


def test_box_points_in_lattice():
    q, lam = 10007, 1234
    for (a, b) in toy.box_points(q, lam, 6):
        assert (a - lam * b) % q == 0 and gcd(a, b) == 1 and (b > 0 or (b == 0 and a > 0))


def test_fiber_preserves_parent_and_labels():
    N, A = _field()
    ctx = toy.Ctx(2 ** 12, 2 ** 16)
    rng = random.Random(3)
    ells = [p for p in range(2 ** 17 + 1, 2 ** 17 + 4000) if toy.is_prime(p)]
    done = 0
    for ell in ells:
        nodes = toy.level1_nodes(A, ell, 16, ctx)
        if not nodes or not nodes[0]["alg"]:
            continue
        nd = nodes[0]
        P0 = (nd["a"], nd["b"])
        assert toy.G_val(A, *P0) % ell == 0
        for (q, r) in nd["alg"]:
            assert toy.F_val(A, *P0) % q == 0 and r in toy.roots_mod(A["f"], q)
        for _ in range(4):
            B = toy.fiber_field(A, P0, rng, 3)
            if B is None:
                continue
            assert toy.F_val(B, *P0) == toy.F_val(A, *P0)
            assert toy.hom_eval(B["f"], A["g0"], A["g1"]) == N * A["w"]
            for (q, r) in nd["alg"]:
                if (B["f"][-1] * B["disc"]) % q:
                    assert r in toy.roots_mod(B["f"], q)
        done += 1
        if done >= 3:
            break
    assert done >= 1


def test_line_keeps_parent_rational_value():
    N, A = _field()
    a0, b0 = 1234, 777
    G0 = toy.G_val(A, a0, b0)
    for j in range(-5, 6):
        g0, g1 = A["g0"] + j * a0, A["g1"] + j * b0
        assert g1 * a0 - g0 * b0 == G0
