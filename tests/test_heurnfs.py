"""Tests for src/heurnfs (agent heurnfs): closed forms, numeric optima, exact height LP, barriers."""
from fractions import Fraction as Fr

import pytest
import sympy as sp

from src.heurnfs import models as M


def test_closed_forms_known_values():
    assert abs(float(M.gnfs_closed()) - 1.9229994) < 1e-6
    assert sp.simplify(M.gnfs_closed() ** 3 - sp.Rational(64, 9)) == 0
    assert abs(float(M.mnfs_closed()) - 1.9018836) < 1e-6
    assert abs(float(M.gnfs_closed(1)) - 1.5262857) < 1e-6          # SNFS (32/9)^{1/3}
    assert sp.simplify(M.tdim_closed(2) - M.gnfs_closed()) == 0       # t = 2 is classical GNFS
    assert float(M.tdim_closed(3)) > float(M.tdim_closed(2))          # higher-dim sieving is worse
    assert sp.simplify(M.init_split_closed() ** 3 - 3) == 0


def test_init_split_numeric_matches():
    assert abs(M.init_split_numeric() - float(M.init_split_closed())) < 1e-6


@pytest.mark.parametrize("kappa", [1.0, 2.0, 3.0])
def test_gnfs_numeric_matches_closed(kappa):
    r = M.gnfs(kappa=kappa)
    assert r["viol"] < 1e-9
    assert abs(r["c"] - float(M.gnfs_closed(kappa))) < 2e-4


def test_tdim3_numeric_matches_closed():
    r = M.tdim(3)
    assert abs(r["c"] - float(M.tdim_closed(3))) < 2e-4


def test_mnfs_numeric_matches_coppersmith():
    r = M.mnfs()
    assert r["viol"] < 1e-9
    assert abs(r["c"] - float(M.mnfs_closed())) < 3e-4
    assert r["test_cost"] < r["c"] - 0.05          # ECM/batch test cost does not bind at the optimum


def test_factory():
    r = M.factory()
    assert abs(r["c"] - 1.6386) < 3e-4
    assert 1.99 < r["precomp"] < 2.02


def test_descent_rows_never_beat_gnfs():
    assert abs(M.descent_rows()["c"] - float(M.gnfs_closed())) < 2e-4
    assert M.descent_rows(c_init=float(M.init_split_closed()))["c"] > 2.1


def test_symmetric_mnfs_worse():
    assert M.mnfs_symmetric()["c"] > 2.3


def test_height_lp_exact():
    for d in range(2, 12):
        assert M.height_lp(1, d) == Fr(2, d + 2)
        assert M.height_lp(d, d) == Fr(1, d)
    assert M.height_lp(2, 10) == Fr(2, 13)


def test_tower_ledger_values():
    # optimistic tower model: equals the 1-D monic-g LP (eta does not enter by construction)
    assert M.tower_ledger(3, 5)["heights"] == Fr(9, 29)
    # genuine TNFS analogue is strictly worse for eta >= 2
    assert M.tnfs_genuine_coeff(1) == 2 and all(M.tnfs_genuine_coeff(e) > 2 for e in (2, 3, 6))


def test_height_lp_sum_bound():
    # summing the two constraints: (D+1)(h1+h2) >= 2, for every degree pair
    for d1 in range(1, 9):
        for d2 in range(max(d1, 2), 15):
            assert M.height_lp(d1, d2) >= Fr(2, d1 + d2 + 1)


def test_omega_curve_mnfs_le_gnfs():
    cur = M.la_omega_curve(omegas=(1.5, 2.373))
    for w, v in cur.items():
        assert v["mnfs"] <= v["gnfs"] + 1e-9
    assert cur[1.5]["mnfs"] < 1.87 and cur[1.5]["mnfs_no_test_term"] < cur[1.5]["mnfs"] - 0.02


def test_alpha_floors():
    assert M.alpha_floor_generic() == Fr(1, 3)
    assert M.alpha_floor_one_level_descent() == Fr(2, 5)
    assert M.alpha_floor_aggregate() == Fr(1, 4)
    # generic ledger: aggregate with r ~ 1 (r_exp = 0) and best d_exp gives back 1/3
    best = min(M.alpha_floor_aggregate(Fr(0), Fr(k, 60)) for k in range(1, 60))
    assert best == Fr(1, 3)
