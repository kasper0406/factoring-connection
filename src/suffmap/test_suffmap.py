"""Exact sanity checks for notes/suffmap_results.md. Run: .venv/bin/python -m pytest -q src/suffmap/test_suffmap.py"""
from fractions import Fraction as Fr

import sympy as sp

from src.suffmap import listing as L
from src.suffmap import principle as P


def test_fixed_size_optimum_matches_closed_form():
    for al, s in [(Fr(1, 2), 1), (Fr(1, 2), Fr(1, 2)), (Fr(1, 4), 1), (Fr(1, 5), Fr(5, 4)), (Fr(1, 10), 3)]:
        opt = P.optimum_fixed_size(al, s)
        assert sp.simplify(opt["c"] - P.c_sampling(al, s)) == 0
        assert sp.simplify(opt["trials"] - opt["linear_algebra"]) == 0  # balanced at the optimum


def test_dixon_and_qs_constants():
    assert sp.simplify(P.c_sampling(Fr(1, 2), 1) - sp.sqrt(2)) == 0      # Dixon L[1/2, sqrt 2]
    assert P.c_sampling(Fr(1, 2), Fr(1, 2)) == 1                        # QS L[1/2, 1]


def test_sigma_table():
    assert [P.sigma_needed(a) for a in P.TARGETS] == [Fr(2, 3), Fr(1, 2), Fr(2, 5), Fr(1, 5)]
    assert P.s_needed(Fr(1, 4), 1) == 1 and P.s_needed(Fr(1, 5), 1) == sp.Rational(5, 4)


def test_size_entropy_laws_reproduce_known_constants():
    t = P.law_table()
    assert t["Dixon"]["alpha"] == Fr(1, 2) and sp.simplify(t["Dixon"]["c"] - sp.sqrt(2)) == 0
    assert t["QS"]["c"] == 1
    assert t["GNFS"]["alpha"] == Fr(1, 3)
    assert sp.simplify(t["GNFS"]["c"] ** 3 - sp.Rational(64, 9)) == 0
    assert sp.simplify(t["SNFS"]["c"] ** 3 - sp.Rational(32, 9)) == 0
    assert sp.simplify(t["GNFS t=3"]["c"] ** 3 - sp.Rational(128 * 2, 27)) == 0   # heurnfs (128(t-1)/(9t))^{1/3}
    assert t["aggregate rank-r (stylised, A=1)"]["alpha"] == Fr(1, 4)
    for v in t.values():
        assert v["balance_ratio"] == 1


def test_law_exponent_rule():
    assert P.exponent_from_law(Fr(1, 2), Fr(1, 2)) == Fr(1, 3)
    assert P.exponent_from_law(Fr(2, 5), 0) == Fr(1, 5)
    assert P.exponent_from_law(Fr(3, 10), Fr(1, 2)) == Fr(1, 5)
    assert P.exponent_from_law(Fr(1, 5), 0) == Fr(1, 10)
    assert P.log_power_of_law(Fr(1, 2), Fr(1, 2)) == 0


def test_one_sided_free_lister_equals_sampling():
    for al, nu in [(Fr(1, 4), 1), (Fr(1, 5), 2), (Fr(1, 3), Fr(1, 2))]:
        assert sp.simplify(P.one_sided_optimum_free_lister(al, nu) - P.c_sampling(al, nu)) == 0


def test_box_exponents():
    assert P.congruence_box_exponent(Fr(1, 4)) == Fr(3, 4)
    assert P.nfs_doubly_box_exponent(Fr(1, 3)) == Fr(1, 3)
    assert P.nfs_doubly_box_exponent(Fr(1, 10)) == Fr(4, 5)


def test_lemma_sl_call_bound():
    for pick in ("first", "last", "middle"):
        r = L.check_lemma_sl(10 ** 6, 10 ** 6 + 1500, 40, pick)
        assert r["complete"] and r["calls"] <= r["bound"]
    r = L.check_lemma_sl(2 ** 31, 2 ** 31 + 300, 7)   # sparse case, possibly K = 0
    assert r["complete"] and r["calls"] <= r["bound"]


def test_rational_rows_bound():
    r = L.rational_side_rows(m=10 ** 5 + 3, E=150, y=50, bmax=12)
    assert r["calls"] <= r["bound"]


def test_congruence_lower_bound():
    for N in (15, 21, 35, 77):
        assert L.congruence_max_lower_bound(N)


def test_single_interval_program():
    r = P.single_interval_program()
    assert sp.simplify(r["c_cubed"] - sp.Rational(16, 9)) == 0          # (16/9)^{1/3} = 1.2114
    assert sp.simplify(r["kappa"] ** 3 - sp.Rational(2, 9)) == 0
    assert r["intersection_left_of_min"]
    assert min((P.single_interval_exponent(Fr(k, 60)), Fr(k, 60)) for k in range(1, 60))[1] == Fr(1, 3)


def test_promise_constructor_bounds():
    lo, hi, y = 10 ** 5, 10 ** 5 + 4000, 40
    S = [n for n in range(lo, hi + 1) if L.is_smooth(n, y)]
    G = max(b - a for a, b in zip([lo - 1] + S, S + [hi + 1]))      # max gap incl. ends
    h0 = 4 * G
    H = hi - lo + 1
    adv = L.AdversarialPromiseOracle(S, lo, hi, h0)
    # tiling, then sliding windows: the adversary caps distinct outputs at 4H/h0 + 1
    outs = set(L.tiling_outputs(lo, hi, h0, adv))
    for a in range(lo, hi - h0 + 2, max(1, h0 // 7)):
        outs.add(adv(a, a + h0 - 1))
    assert len(outs) <= 4 * H // h0 + 1
    # an honest (leftmost) constructor tiled gives >= floor(H/h0) distinct outputs
    honest = lambda a, b: min(s for s in S if a <= s <= b)
    assert len(set(L.tiling_outputs(lo, hi, h0, honest))) >= H // h0
    revealed = sorted(adv.revealed)
    assert all(b - a >= h0 // 4 for a, b in zip(revealed, revealed[1:]))
