"""Tests for src/thmAB (notes/thmAB_results.md): exact checks of the finite statements of the algebraic cell device."""
from fractions import Fraction as Fr
import random

from src.thmAB import cells, params


def test_margins_hold_on_window_and_endpoint():
    for al in [Fr(3, 8), Fr(19, 50), Fr(7, 18), Fr(39, 100), Fr(399, 1000)]:
        assert params.all_ok(al), (al, params.margins(al))
    assert params.acceptability_ok(Fr(3, 8)) == "needs a > 2 c2 (O4)"
    assert params.acceptability_ok(Fr(7, 18)) == "ok"
    # theta = 2 alpha - 1/2 on [3/8, 2/5); separate parent thresholds would fail if shared for alpha >= 7/18
    for al in [Fr(3, 8), Fr(7, 18), Fr(2, 5) - Fr(1, 1000)]:
        assert params.theta(al) == 2 * al - Fr(1, 2)
    assert params.theta(Fr(7, 18)) >= 3 - 7 * Fr(7, 18)          # shared K0 fails the parent saving (O3)
    assert params.theta_par(Fr(7, 18)) < 3 - 7 * Fr(7, 18)


def test_class_bijection_CL():
    for Dp, kappa, X in [(12, 5, 30), (8, 3, 25), (9, 2, 27), (30, 7, 35), (1, 0, 5), (25, 3, 30)]:
        if Dp == 1:
            kappa = 0
        out = cells.class_bijection(Dp, kappa, X)
        assert out["points"] > 0


def test_cell_factorisation_CF():
    rng = random.Random(7)
    n = 0
    for Dp in [6, 8, 12, 15, 18, 30]:
        units = [k for k in range(1, Dp) if cells.gcd(k, Dp) == 1]
        for delta in [t for t in range(1, Dp + 1) if Dp % t == 0]:
            for _ in range(3):
                kappa = rng.choice(units)
                lo1 = rng.randint(-400, 300)
                lo2 = rng.randint(-400, 300)
                cells.cf_check(Dp, kappa, delta, (lo1, lo1 + rng.randint(60, 200)), (lo2, lo2 + rng.randint(60, 200)), 7)
                n += 1
    assert n > 50


def test_parallelogram_geometry_PG():
    rng = random.Random(11)
    done = 0
    while done < 60:
        d = rng.choice([3, 4, 5])
        P1 = (rng.randint(-6, 6), rng.randint(1, 6))
        P2 = (rng.randint(-6, 6), rng.randint(1, 6))
        if cells.gcd(*P1) != 1 or cells.gcd(*P2) != 1:
            continue
        if P1[0] * P2[1] - P1[1] * P2[0] == 0:
            continue
        H = rng.randint(5, 40)
        G1, G2 = rng.randint(1, 9), rng.randint(1, 9)
        cells.pg_check(P1, P2, d, H, G1, G2)
        done += 1
    # nearly parallel pairs (large kappa)
    for P1, P2 in [((50, 49), (51, 50)), ((7, 1), (8, 1)), ((100, 1), (99, 1))]:
        cells.pg_check(P1, P2, 3, 10, 3, 2)


def test_boundary_area_BL():
    for P1, P2, d in [((1, 2), (3, 1), 3), ((5, 4), (6, 5), 3), ((2, 3), (1, 5), 4)]:
        poly, _ = cells.parallelogram(P1, P2, d, 3, 2, 1)
        for h in [Fr(1, 4), Fr(1, 10), Fr(1, 25)]:
            cells.boundary_area_check(poly, h)


def test_fibre_masses_PG_c():
    for P1, P2, d, H, G1, G2 in [((1, 2), (3, 1), 3, 12, 2, 3), ((2, 1), (1, 3), 3, 15, 1, 2),
                                 ((3, 2), (1, 4), 3, 10, 1, 1)]:
        Det = abs(P1[0] * P2[1] - P1[1] * P2[0])
        lam = next(l for l in range(Det) if (l * P2[0] - P1[0]) % Det == 0 and (l * P2[1] - P1[1]) % Det == 0)
        out = cells.fibre_mass_check(P1, P2, d, H, G1, G2, pow(lam, d, Det), Det)
        assert out["Z"] > 0
