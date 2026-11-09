"""Tests for src/caudit/coverage.py (agent caudit). Exact arithmetic; no leakage (no factor-aware code)."""
from fractions import Fraction as Fr

import sympy as sp

from src.caudit import coverage as cov


def test_sub_third_wall():
    # NFS sits exactly at 1/3; the wall slack is 0.
    assert cov.sub_third_wall(Fr(1, 2), Fr(1, 2))["alpha"] == Fr(1, 3)
    assert cov.sub_third_wall(Fr(1, 2), Fr(1, 2))["wall_slack"] == 0
    # To get below 1/3 at b = 1/2 you need a < 1/2 (smaller norm growth).
    assert cov.sub_third_wall(Fr(2, 5), Fr(1, 2))["below_third"] is True
    assert cov.sub_third_wall(Fr(1, 2), Fr(1, 2))["below_third"] is False
    # Dixon/QS at (a,b)=(1,0) give 1/2.
    assert cov.sub_third_wall(Fr(1), Fr(0))["alpha"] == Fr(1, 2)


def test_catalogue_exponents_and_the_formal_sub_third_hole():
    tab = cov.catalogue_exponents()
    assert tab, "catalogue is non-empty"
    for name, d in tab.items():
        assert d["balance_ok"], f"{name} fails the size-entropy balance identity"
    # Every REALIZED construction (Dixon, QS, the NFS family) sits at or above 1/3.
    for name in ("Dixon", "QS", "GNFS", "SNFS", "GNFS t=3"):
        assert tab[name]["at_or_above_third"], f"{name} claims alpha < 1/3"
    assert tab["GNFS"]["alpha"] == Fr(1, 3)
    assert tab["SNFS"]["alpha"] == Fr(1, 3)
    # KEY AUDIT POINT: the law itself does NOT forbid sub-1/3. The stylised aggregate rank-r law (a=1/2, b=0)
    # formally gives alpha = 1/4 and passes the balance identity. Theorem SS is a reduction, not a prohibition:
    # the barrier against it is the ABSENCE of a family with that law (resbarrier G_r), believed-false only.
    assert tab["aggregate rank-r (stylised, A=1)"]["alpha"] == Fr(1, 4)
    assert tab["aggregate rank-r (stylised, A=1)"]["at_or_above_third"] is False


def test_mnfs_constants_exponent_one_third():
    g = float(cov.gnfs_constant())
    m = float(cov.mnfs_constant())
    assert abs(g - 1.9229994) < 1e-6
    assert abs(m - 1.9018836) < 1e-6
    d = cov.outside_S_is_constant_only()
    assert d["exponent"] == Fr(1, 3)
    assert d["mnfs_below_gnfs"] is True
    assert d["factory_c"] < m  # factory improves the constant further, still exponent 1/3


def test_prop_si_constant():
    c = cov.prop_si_constant()
    assert sp.simplify(c - sp.Rational(16, 9) ** sp.Rational(1, 3)) == 0
    chk = cov.prop_si_check()
    assert abs(chk["c_float"] - 1.2114137) < 1e-6
    assert chk["exponent"] == Fr(1, 3)
    assert chk["below_mnfs"] is True
    assert chk["matches_suffmap"], "Prop SI constant disagrees with suffmap's single_interval_program"


def test_leak_floor_threshold():
    # 2048-bit N: Coppersmith needs 512 known top bits of p.
    d = cov.coppersmith_known_bits_threshold(2048)
    assert d["bits_needed"] == 512
    assert d["fraction"] == Fr(1, 4)
    # ceil behaviour on a non-multiple of 4.
    assert cov.coppersmith_known_bits_threshold(2050)["bits_needed"] == 513


def test_cross_family_union_is_unproven():
    d = cov.cross_family_union_status()
    assert d["status"] == "UNPROVEN EXTRAPOLATION"
