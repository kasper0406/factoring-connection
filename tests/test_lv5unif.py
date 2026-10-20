"""Tests for src/lv5unif (notes/lv5unif_results.md). Exact arithmetic only."""
from fractions import Fraction as Fr
from math import gcd
import random

from src.lv5unif.kern import (fibre_count_lattice, fibre_histogram_brute, kernel_basis, phi,
                              triangular_bound, vvec)
from src.lv5unif import params as P
from src.thmA.params import h_star, theta as theta_thmA


def _pairs(rnd, V, n):
    out = []
    while len(out) < n:
        b, bp = rnd.randint(1, V), rnd.randint(1, V)
        a, ap = rnd.randint(-V, V), rnd.randint(-V, V)
        if gcd(a, b) == 1 and gcd(ap, bp) == 1 and a * bp - ap * b != 0:
            out.append(((a, b), (ap, bp)))
    return out


def test_kernel_basis_is_in_kernel():
    rnd = random.Random(1)
    for d in range(3, 12):
        for Pp in _pairs(rnd, 9, 10):
            v, vp = vvec(d, *Pp[0]), vvec(d, *Pp[1])
            for e in kernel_basis(d, *Pp):
                assert phi(e, v, vp) == (0, 0)


def test_lattice_count_equals_brute_force_and_bound():
    rnd = random.Random(2)
    for d, H, V in [(3, 6, 3), (3, 9, 4), (4, 4, 3), (4, 5, 2)]:
        for Pp in _pairs(rnd, V, 6):
            hist = fibre_histogram_brute(d, Pp[0], Pp[1], H)
            bound = triangular_bound(d, Pp[0][1], Pp[1][1], H)
            assert max(hist.values()) <= bound
            for _ in range(5):
                c0 = [rnd.randint(-H, H) for _ in range(d)]
                z = phi(c0, vvec(d, *Pp[0]), vvec(d, *Pp[1]))
                assert fibre_count_lattice(d, Pp[0], Pp[1], H, c0) == hist[z]


def test_probability_form_of_kern_u():
    # Pr <= e/((2H+1)^2 (bb')^{d-2}) when bb' <= (2H+1)/d ; check (2H+1)^2 (bb')^{d-2} maxfib/(2H+1)^d <= 2.72
    rnd = random.Random(3)
    for d, V in [(3, 2), (4, 2), (3, 3)]:
        H = d * V * V
        for Pp in _pairs(rnd, V, 5):
            b, bp = Pp[0][1], Pp[1][1]
            if b * bp * d > 2 * H + 1:
                continue
            hist = fibre_histogram_brute(d, Pp[0], Pp[1], H)
            stat = Fr((2 * H + 1) ** 2 * (b * bp) ** (d - 2) * max(hist.values()), (2 * H + 1) ** d)
            assert stat <= Fr(272, 100)


def test_two_pow_bracket_and_factor():
    lo, hi = P.two_pow_2_5_bracket()
    assert lo ** 5 < 4 < hi ** 5 and hi - lo == Fr(1, 10 ** 6)
    flo, fhi = P.recursion_factor_bracket()
    assert Fr(412, 100) < flo <= fhi < Fr(414, 100)


def test_theta0_continuous_at_regime_boundary():
    for d0 in [Fr(1), Fr(2), Fr(3, 2), Fr(7, 3)]:
        kb = d0 * d0 / 4
        assert P.theta0(d0, kb) == Fr(3) / d0
        assert Fr(2) / d0 + d0 / (4 * kb) == Fr(3) / d0


def test_all_margins_hold():
    for a, d0, kappa in [(Fr(1), Fr(2), Fr(1, 2)), (Fr(3, 5), Fr(1), Fr(1, 8)), (Fr(2), Fr(1, 3), Fr(5)),
                         (Fr(1, 4), Fr(4), Fr(1, 100))]:
        for name, ok, where in P.margins(a, d0, kappa):
            assert ok, (name, a, d0, kappa)


def test_cost_constant_examples():
    lo, hi = P.cost_constant(Fr(3, 5), Fr(1), Fr(1, 8), Fr(0))
    # s = 1/4, sigma = 2*4/(5*3/5) = 8/3 ; p = a + s + sigma = 211/60 ; rec = (17/20)*4.13 = 3.51
    assert lo == hi == Fr(3, 5) + Fr(1, 4) + Fr(8, 3)
    lo, hi = P.cost_constant(Fr(1), Fr(1), Fr(1, 8), Fr(5))
    # s = 2*(1/8)*6 = 3/2 ; recursion term (1 + 3/2)*4.1298 = 10.32 binds
    assert lo < hi and Fr(1032, 100) < lo and hi < Fr(1033, 100)


def test_finite_N_parameters_match_asymptotics():
    # For a large 'D' (as an exact integer), the exact finite-(d,k) theta (thmA formula at mu_opt, eta0 = 0)
    # times D^{1/5} is within 5% of theta0(d0, kappa).
    for D, d0, kappa in [(10 ** 15, Fr(2), Fr(1, 2)), (10 ** 20, Fr(1), Fr(1, 8)), (10 ** 25, Fr(3), Fr(4))]:
        D15 = Fr(round(D ** 0.2))      # D^{1/5} to integer accuracy (only used to define d, k; not a decision)
        d = int(d0 * D15) | 1
        k = max(int(kappa * D15 * D15) + 1, d + 1)
        nu = Fr(1, 4 * k)
        mus = [(1 - nu) / d, (1 + (d - 1) * nu) / (d + 1)]
        th = min(theta_thmA(d, k, mu) for mu in mus if mu >= nu)
        ratio = th * D15 / P.theta0(d0, kappa)
        assert Fr(95, 100) < ratio < Fr(105, 100), (D, d0, kappa, float(ratio))


def test_tau_table_all_zero():
    assert all(row[-1] == 0 for row in P.TAU_TABLE)
