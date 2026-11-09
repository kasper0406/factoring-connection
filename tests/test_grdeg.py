"""Tests for agent grdeg (degeneracy-driven G_r failure, frame-less type-D)."""
import math
from fractions import Fraction as Fr

from src.grdeg import construct as C
from src.grdeg import ledger as L
from src.grdeg import obstruction as O
from src.rank3 import forms as RF


# --- ledger (deliverable #1) --------------------------------------------------

def test_ledger_aggregate_is_quarter():
    assert L.aggregate_alpha()["alpha"] == Fr(1, 4)


def test_ledger_r4_is_third():
    assert L.r4_alpha()["alpha"] == Fr(1, 3)


def test_ledger_crossover_endpoints():
    # c=0 (floor 1/d, independent of r) -> 1/4 ; c=1 (floor ~r/d, R4) -> 1/3
    assert L.crossover_floor(Fr(0)) == Fr(1, 4)
    assert L.crossover_floor(Fr(1)) == Fr(1, 3)
    # monotone increasing in c
    vals = [L.crossover_floor(Fr(k, 4)) for k in range(5)]
    assert all(vals[i] < vals[i + 1] for i in range(4))


# --- construct (E1) -----------------------------------------------------------

def test_D1_is_frame_free_degenerate():
    mod = C.weil_restriction_q2q3(1001)
    geo = C.geometry(mod)
    assert geo["Q2_is_zero"] is True          # Delta_F = 0
    assert geo["frame_free"] is True          # no r+1 rows in general position
    assert geo["num_frames"] == 0


def test_D2_and_GP_are_framed():
    for mod in (C.cyclotomic_z5(1001), C.gp_control_cubics(1001)):
        geo = C.geometry(mod)
        assert geo["Q2_is_zero"] is False
        assert geo["frame_free"] is False


def test_D1_small_norms_only_on_sublattice():
    # genuine rank-4 candidates (shift active) pay |W'| ~ N^Theta(1); small norms live on a3=a4=0 sublattice.
    E = 4
    sub_mins, gen_over_N = [], []
    for N in (101 * 103, 10007 * 10009):
        mod = C.weil_restriction_q2q3(N)
        sub_min = gen_min = None
        for a in RF.box_points(4, E):
            w, w2 = RF.F_value(mod["f"], mod["P"], mod["g"], mod["Q"], a)
            if w == 0 or w2 == 0:
                continue
            v = abs(w * w2)
            if a[2] == 0 and a[3] == 0:
                sub_min = v if sub_min is None else min(sub_min, v)
            else:
                gen_min = v if gen_min is None else min(gen_min, v)
        sub_mins.append(sub_min)
        gen_over_N.append(gen_min / N)
    # sublattice minimum is N-independent (same for both N); genuine minimum grows faster than N
    assert sub_mins[0] == sub_mins[1]
    assert all(x > 1 for x in gen_over_N)       # genuine |F| >> N
    assert gen_over_N[1] > gen_over_N[0]          # scales up with N


def test_D1_coincidence_sublattice_is_rank2_gp():
    # the small-norm sublattice a3=a4=0 is where the two sides COINCIDE (type-A on that plane):
    # there W'=W, so it is a rank-2 general-position sub-family (a linear exception, rankdegen LP_exc(ii)).
    mod = C.weil_restriction_q2q3(10007 * 10009)
    for a in RF.box_points(4, 3):
        if a[2] == 0 and a[3] == 0 and (a[0] or a[1]):
            w, w2 = RF.F_value(mod["f"], mod["P"], mod["g"], mod["Q"], a)
            assert w == w2            # sides coincide on the sublattice


# --- obstruction (E2) ---------------------------------------------------------

def test_link_requirement_is_sqrt_for_quadratic():
    req = O.link_requirement(7 * 11, [-2, 0, 1])
    assert req["subfield_degree_e"] == 2
    assert "factoring" in req["hardness"].lower()
    assert req["charged_subroutine"] == "root_mod_N(mu_k, N)"


def test_rabin_split_structural():
    # p,q == +-1 mod 8 so 2 is a QR at both; CRT-inequivalent root splits N, diagonal does not.
    p, q = 1000033, 2000039  # both == 1 or 7 mod 8
    assert pow(2, (p - 1) // 2, p) == 1 and pow(2, (q - 1) // 2, q) == 1
    out = O.rabin_split_harness(2, p, q)
    assert out["splittable"]
    byname = {r["pair"]: r for r in out["results"]}
    mixed = byname["mixed root (CRT-inequivalent)"]
    assert mixed["gcd"] in (p, q) and mixed["class"].startswith("STRUCTURAL")
    diag = byname["-s (diagonal)"]
    assert diag["gcd"] == 1 and diag["class"] == "no split"


def test_naive_root_search_fails_for_semiprime():
    # honest operational cost: no root of x^2-2 mod N found below sqrt(smallest prime)
    p, q = 1000033, 2000039
    out = O.naive_link_attempt(p * q, [-2, 0, 1], search_bound=10 ** 4)
    assert out["found"] is False
