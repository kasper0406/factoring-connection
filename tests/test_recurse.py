"""Tests for agent recurse (src/recurse)."""
import json
from fractions import Fraction as Fr

import sympy as sp

from src.recurse import constants as C
from src.recurse import pipeline_rec as pr
from src.thmA import params as tpar


def test_fixpoint_identity_exact():
    for th in (Fr(3, 10), Fr(16, 39), Fr(1, 7), Fr(2, 5)):
        t = sp.Rational(th.numerator, th.denominator)
        c, a = C.cstar(th), C.astar(th)
        # both relation terms balance at a*: a + theta/(2a) = (2+sqrt2) a = c*
        assert sp.simplify(a + t / (2 * a) - c) == 0
        assert sp.simplify((2 + C.SQ2) * a - c) == 0
        # fixpoint of c = a + c/sqrt2 and of g(c) = (theta/c + c)/sqrt2
        assert sp.simplify(a + c / C.SQ2 - c) == 0
        assert sp.simplify((t / c + c) / C.SQ2 - c) == 0
        # linear algebra slack: 2a*/c* = 2 - sqrt2
        assert sp.simplify(2 * a / c - (2 - C.SQ2)) == 0


def test_c_of_a_minimised_at_astar():
    th = Fr(3, 10)
    t = sp.Rational(3, 10)
    a0 = C.astar(th)
    c0 = C.cstar(th)
    for q in (Fr(1, 10), Fr(1, 4), Fr(7, 30), Fr(2, 9), Fr(1, 2)):
        a = sp.Rational(q.numerator, q.denominator)
        val = sp.Max(a + t / (2 * a), (2 + C.SQ2) * a)
        assert sp.N(val - c0, 50) >= 0
    assert abs(float(a0) - 0.2493) < 1e-3


def test_rational_certificates_match_radicals():
    for th in (Fr(16, 39), Fr(23, 54), Fr(3, 10), Fr(5, 12), Fr(2, 5), Fr(41, 99)):
        assert C.rec_lt_1(th) == bool(sp.N(C.cstar(th), 60) < 1)
        cA = (1 + sp.Rational(th.numerator, th.denominator)) / C.SQ2
        assert C.rec_beats_A(th) == bool(sp.N(C.cstar(th) - cA, 60) < 0)
        # c* < 1 has the same threshold as Theorem A's (1+theta)/sqrt2 < 1
        assert C.rec_lt_1(th) == tpar.c_lt_1(th)


def test_fixed_depth_starts_at_LP_and_theorem_A():
    th = Fr(3, 10)
    cs = C.fixed_depth(th, 6)
    assert cs[0] == 1
    assert sp.simplify(cs[1] - (1 + sp.Rational(3, 10)) / C.SQ2) == 0
    vals = [float(c) for c in cs]
    assert all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
    assert vals[-1] > float(C.cstar(th))
    # geometric contraction factor sqrt2 - 1 at the fixed point
    r = (vals[-1] - float(C.cstar(th))) / (vals[-2] - float(C.cstar(th)))
    assert abs(r - (2 ** 0.5 - 1)) < 0.02


def test_table_d9():
    rows = {r["d"]: r for r in C.table(ds=(9,), ks=(50,))}
    r = rows[9]
    assert r["kmin"] == 13 and r["theta_kmin"] == Fr(16, 39)
    assert abs(float(r["c_rec_inf"]) - 0.85104) < 1e-5


def test_oracle_substitution_toy():
    cfg = json.load(open("experiments/configs/thmA_e1.json"))["pipeline"]
    from src.common import instances as I
    inst = I.balanced_semiprime(46, 900)
    a = pr.factor(inst.N, cfg, 900, "pb")
    b = pr.factor(inst.N, cfg, 900, "rec")
    assert a["factors"] == b["factors"] == sorted(inst.factors)
    assert pr.outer_transcript(a) == pr.outer_transcript(b)
