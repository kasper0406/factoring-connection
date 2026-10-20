"""Tests for agent smoothenum (target programs, barrier helpers, mechanisms)."""
import random
from math import gcd

import pytest
import sympy as sp

from src.smoothenum import target as T
from src.smoothenum.barrier import lq_gain, rho, wagner_max_k
from src.smoothenum.mech import (Counters, bottom_up_hit, crt_root, lattice_count, lattice_points,
                                 lattice_sample, reduced_basis)
from src.smoothenum.nfsdata import SideTest, base_m_pair, primes_upto


def test_factory_closed_form():
    c = T.factory_closed()
    assert abs(float(c) - 1.6385862) < 1e-6
    assert sp.simplify(c ** 3 - (20 + 8 * sp.sqrt(6)) / 9) == 0
    assert abs(float(T.factory_precomp_closed()) - 2.0068501) < 1e-6


def test_alg_infimum_equals_cuberoot_4_9():
    assert sp.simplify(T.alg_output_infimum() - sp.cbrt(sp.Rational(4, 9))) == 0


def test_costs_at_published_mnfs_optimum():
    b = ((46 + 13 * 13 ** 0.5) / 108) ** (1 / 3)
    x = [b, b, b * (13 ** 0.5 - 1) / 3, b * (4 - 13 ** 0.5) / 3, b * (4 * 13 ** 0.5 - 10) / 3]
    c = T.costs(x, "box")
    assert abs(c["enum"] - T.C_MNFS) < 1e-6 and abs(c["la"] - T.C_MNFS) < 1e-6
    assert abs(c["yield_slack"]) < 1e-6


def test_factory_point_feasible_in_rat_mode():
    beta = ((2 + 6 ** 0.5) / 6) ** (2 / 3)
    delta = beta ** -0.5
    e = beta + 1 / (6 * beta * delta)
    c = T.costs([e, beta, beta, 0.0, delta], "rat", theta=1.0)
    assert abs(max(c["enum"], c["test"], c["la"]) - 1.6385862) < 1e-6
    assert c["yield_slack"] > -1e-9


def test_doubly_oracle_any_h():
    # codex construction: beta = beta_a = h, delta = 3h, r = 1/(9h^2), e = h + 2r gives cost 2h
    for h in (0.05, 0.2, 0.5):
        x = [h + 2 / (9 * h * h), h, h, 0.0, 3 * h]
        c = T.costs(x, "both")
        assert abs(c["la"] - 2 * h) < 1e-12 and c["yield_slack"] > -1e-9


def test_rho_values():
    assert abs(rho(2) - (1 - sp.log(2).evalf())) < 1e-7
    assert abs(rho(10) / 2.77017183772596e-11 - 1) < 1e-3
    assert lq_gain(5.0, 2.0) > 1


def test_wagner_only_k2_asymptotically():
    assert wagner_max_k(10000.0, 100.0) == 2
    assert wagner_max_k(50.0, 20.0) >= 4


def test_side_test():
    t = SideTest(1000, 12, 24)
    assert t.ok(2 * 3 * 997 * 1009 * 4001)
    assert not t.ok(4099 * 4091)
    assert t.ok(1009 * 4001)


def _naive(R, M, A, Bm):
    return sorted((a, b) for b in range(1, Bm + 1) for a in range(-A, A + 1) if (a - R * b) % M == 0)


@pytest.mark.parametrize("R,M,A,Bm", [(5, 7, 30, 10), (1234, 10007, 2000, 40), (3, 97 * 89, 500, 60),
                                      (77, 101 * 103 * 107, 4000, 300)])
def test_lattice_enumeration_exact(R, M, A, Bm):
    cnt = Counters()
    pts = sorted(lattice_points(R, M, A, Bm, cnt))
    assert pts == _naive(R, M, A, Bm)
    assert lattice_count(R, M, A, Bm) == len(pts)
    smp = lattice_sample(R, M, A, Bm, 50, random.Random(1))
    assert all(p in set(pts) for p in smp)


def test_crt_root():
    R, M = crt_root([(3, 7), (5, 11), (2, 13)])
    assert M == 1001 and R % 7 == 3 and R % 11 == 5 and R % 13 == 2


def test_bottom_up_unique_hit():
    P = base_m_pair(2 ** 61 - 1, 3)
    m = -P.Y0
    A, Bm = 200, 50
    box = {(a, b) for a in range(-A, A + 1) for b in range(1, Bm + 1)}
    vals = {}
    for a, b in box:
        vals.setdefault(P.G(a, b), []).append((a, b))
    assert all(len(v) == 1 for v in vals.values())  # G injective on the box since 2A < m
    rng = random.Random(3)
    for _ in range(300):
        a, b = rng.choice(sorted(box))
        assert bottom_up_hit(P, P.G(a, b), A, Bm) == (a, b)
    assert bottom_up_hit(P, P.G(A + 5, 3), A, Bm) is None or abs(P.G(A + 5, 3)) > 0


def test_theta_closed_form_endpoints_and_numeric():
    assert abs(T.c_rat_theta_closed(0.0) - T.C_MNFS) < 1e-6
    assert abs(T.c_rat_theta_closed(1.0) - float(T.factory_closed())) < 1e-6
    num = T.program("rat", theta=0.5, restarts=15, seed=1)["c"]
    assert num >= T.c_rat_theta_closed(0.5) - 1e-4 and num - T.c_rat_theta_closed(0.5) < 3e-3
    assert abs(T.theta_for_gain(0.01) - 0.04696) < 5e-4
