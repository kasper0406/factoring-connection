"""Tests for src/primality (owner: primality).  Run: .venv/bin/python -m pytest -q tests/test_primality.py"""
import random
from math import gcd

import pytest
import sympy

from src.common.instances import balanced_semiprime
from src.primality import dlog, formulas as F, menu, ops, verify, verify_codex


@pytest.mark.parametrize("primes", [(11, 31), (5, 13), (7, 43), (101, 103), (3, 11, 17), (5, 29, 73), (7, 13, 19), (5, 17, 29)])
def test_power_tests_exhaustive(primes):
    r = verify.check_power_tests(primes)
    assert r["ok"] and r["ok_closed"] is not False
    assert r["tally"]["minus1_final"] == 0  # a^(N-1) = -1 mod N never happens


@pytest.mark.parametrize("pq", [(11, 31), (13, 37), (17, 59), (29, 41)])
def test_lucas_and_codex_pairs(pq):
    r = verify.check_lucas(pq)
    assert "skip" in r or (r["ok"] and r["ok_closed_uniform"])
    c = verify_codex.check_pair(pq)
    assert not c["bad"], c["bad"]
    l5 = verify_codex.check_lucas_codex(pq)
    assert l5 is None or l5["ok"]


def test_frobenius_small():
    for pq in [(7, 11), (5, 23), (13, 17)]:
        r = verify.check_frobenius(pq)
        assert "skip" in r or r["ok"]


def test_torus_all_patterns():
    cases = verify.torus_cases(300, per_pattern=1, seed=3)[:60]
    pats = set()
    for c in cases:
        r = verify.check_torus(c)
        assert r["ok"]
        pats.add(tuple(r["pattern"]))
        rc = verify_codex.check_torus_codex(c)
        assert not rc["bad"]
    assert len(pats) == 4


def test_aks_twist_and_period():
    for c in verify.aks_cases(30, seed=5):
        r = verify.check_aks(c)
        assert r["ok_twist"] and r["ok_period"] is not False


def test_carmichael_561():
    ex = F.mr_exact([3, 11, 17])
    assert ex["prefix"] == sympy.Rational(31, 32) or float(ex["prefix"]) == 31 / 32
    assert F.fermat_exact([3, 11, 17]) == 0


def test_kangaroo_factors_from_fermat_residue():
    for bits in (32, 40):
        inst = balanced_semiprime(bits, 7)
        for meth in ("kangaroo", "bsgs"):
            r = dlog.factor_via_fermat_residue(inst.N, 1, meth)
            assert tuple(sorted(r["factors"])) == tuple(inst.labels["factors"])
        r = dlog.factor_via_lucas_residue(inst.N, 1)
        assert tuple(sorted(r["factors"])) == tuple(inst.labels["factors"])
    r = dlog.factor_via_fermat_residue(10403, 1)
    assert r["x"] == 202 and r["factors"] == (101, 103)


def test_sqrt_minus1_gives_sos():
    rng = random.Random(1)
    ps = [p for p in sympy.primerange(5, 400) if p % 4 == 1]
    for _ in range(40):
        p, q = rng.sample(ps, 2)
        N = p * q
        x = int(sympy.sqrt_mod(-1 % N, N))  # harness-side root
        uv = dlog.sos_from_sqrt_minus1(N, x)
        assert uv is not None and uv[0] ** 2 + uv[1] ** 2 == N


def test_order_multiple_split():
    N = 101 * 103
    hits = 0
    for a in range(2, 200):
        if gcd(a, N) == 1:
            o = sympy.n_order(a, N)
            hits += dlog.split_with_order_multiple(N, a, 7 * o) is not None
    assert hits > 50


def test_menu_identities_and_ec():
    assert menu.check_identities(400)["bad"] == 0
    r = menu.check_ec_formula(seed=2)
    assert r["bad"] == 0 and r["checked"] > 100
