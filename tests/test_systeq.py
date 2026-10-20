"""Tests for src/systeq (agent systeq, phase 16)."""
from fractions import Fraction

import pytest

from src.systeq.ffield import GF, split_count_formula
from src.systeq import bgjt
from src.systeq.intanalog import (descartes_count, ev, is_smooth, phi_n, vfamily_candidates)
from src.systeq.barriers import (additive_exponents, substitution_homography,
                                 systematic_degree_lower_bound, wilson_frobenius_check,
                                 zsigmondy_check)


def test_field_axioms_small():
    F = GF(3, 2)
    for x in range(1, F.Q):
        assert F.m(x, int(F.inv[x])) == 1
    # Frobenius x -> x^3 is additive on F_9
    for x in range(F.Q):
        for y in range(F.Q):
            assert F.pw(F.a(x, y), 3) == F.a(F.pw(x, 3), F.pw(y, 3))


def test_systematic_equation_ff():
    # x^q - x = prod_{a in F_q} (x - a) over F_q, q = 5 inside F_25
    F, q = bgjt.make_field(5, 1)
    roots = F.root_multiset(F.trim([0, int(F.neg[1])] + [0] * (q - 2) + [1]))
    assert sorted(roots) == sorted(x for x in range(F.Q) if F.pw(x, q) == x)
    assert len(roots) == q


def test_split_formula_bruteforce():
    F = GF(3, 1)
    s, t = split_count_formula(3, 2)
    # binary quadratic forms over F_3 up to scalar: count split ones directly
    cnt = 0
    for a in range(3):
        for b in range(3):
            for c in range(3):
                f = F.trim([c, b, a])
                if not f:
                    continue
                if f[-1] != 1:
                    continue
                cnt += F.splits(f)
    # monic deg<=2 forms representative count = t; split = s
    assert t == 13 and s == 10
    assert cnt == s


def test_positive_control_q5():
    s = bgjt.control_summary(5, 1)
    assert s["P_V"] == (25 ** 4 - 1) // 24
    assert s["split_in_V"] >= s["orbit_q3_plus_q"]
    assert s["relations"] == s["verified_mod_I"] > s["Q"] + 1
    assert s["G_FF"] > 10


def test_orbit_elements_split():
    F, q = bgjt.make_field(7, 1)
    # m * (x^q - x) = (cx+d)(ax+b)^q - (ax+b)(cx+d)^q splits for any invertible matrix over F_49
    import random
    rng = random.Random(1)
    for _ in range(20):
        a, b, c, d = (rng.randrange(F.Q) for _ in range(4))
        if F.s(F.m(a, d), F.m(b, c)) == 0:
            continue
        ab, cd = [b, a], [d, c]
        abq = [F.pw(b, q)] + [0] * (q - 1) + [F.pw(a, q)]
        cdq = [F.pw(d, q)] + [0] * (q - 1) + [F.pw(c, q)]
        f = F.psub(F.pmul(F.trim(cd), F.trim(abq)), F.pmul(F.trim(ab), F.trim(cdq)))
        assert f and F.splits(f)


def test_smoothness_exact():
    assert is_smooth(2 ** 40 * 3 ** 5, 3)
    assert not is_smooth(2 * 101, 100)
    assert is_smooth(97 * 89 * 2 ** 30, 100)
    assert is_smooth(-1, 5)


def test_vfamily_congruence():
    N, c = vfamily_candidates(7, 130, 8, 500, 0)
    assert N == ev(phi_n(7), 130)
    assert pow(130, 7, N) == 1
    assert all((v - r) % N == 0 for v, r in c)
    # the systematic side is at least as large as N (Dixon regime), the small side tiny
    assert min(abs(v) for v, _ in c) > N


def test_descartes_barrier():
    for n in (7, 8):
        r = descartes_count(n, 2)
        assert r["split_nondegenerate"] == 0
        assert r["split_degenerate"] == 4 * 5  # A in [-2,2]\{0}, B in [-2,2]
    assert descartes_count(2, 2)["split_nondegenerate"] > 0


def test_prop_F_additive_exponents():
    # A = (Z/15)[t]/(t^2+1) = F_9 x F_5 x F_5 ; additive s need (2|s-1),(4|s-1) and s = 3^i mod 8
    s_list = additive_exponents(15, [1, 0, 1], 40)
    assert s_list == [s for s in range(2, 41) if (s - 1) % 4 == 0 and s % 8 in (1, 3)]
    for s in s_list:
        assert (s - 1) % 2 == 0 and (s - 1) % 4 == 0
    # on Z/N itself (degree-1 algebra) additive power maps are s == 1 mod lambda(N)
    s_list = additive_exponents(21, [0, 1], 30)  # lambda(21) = 6
    assert s_list == [7, 13, 19, 25]


def test_lemma_CF_and_wilson():
    assert wilson_frobenius_check(11, 2)
    assert wilson_frobenius_check(13, 5)
    assert not systematic_degree_lower_bound(list(range(10)), 11)
    assert systematic_degree_lower_bound(list(range(11)), 11)


def test_substitution_homography_collapses():
    for q in (3, 5, 7):
        assert substitution_homography(q)


def test_integer_systematic_small_q():
    # the only integer systematic equations: T^2 - T = T(T-1), T^3 - T = (T-1)T(T+1)
    from flint import fmpz_poly
    for q in range(2, 9):
        _, facs = fmpz_poly([0, -1] + [0] * (q - 2) + [1]).factor()
        allsplit = all(g.degree() == 1 for g, _ in facs)
        assert allsplit == (q <= 3)


def test_zsigmondy_antismooth():
    hits, tot = zsigmondy_check(16, 12, 10)
    assert tot == 110 and hits == 0
