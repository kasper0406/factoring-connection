"""Tests for agent `newfamily` (phase 24)."""
import ast
import os
import random
from fractions import Fraction
from math import gcd

import pytest

from src.newfamily.recip import cf, dedekind6, jacobi
from src.newfamily.modsym import modsym_plus, setup
from src.newfamily import diag

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_dedekind_reciprocity_matches_brute():
    for c in range(2, 90):
        for a in range(1, c):
            if gcd(a, c) == 1:
                assert Fraction(dedekind6(a, c), 6 * c) == diag.dedekind_brute(a, c)


def test_dedekind_classical_congruences():
    # 12 c s(a,c) = a + a^{-1} (mod c) for 3 not dividing c; 12 c s(a,c) = c + 1 - 2 (a/c) (mod 8) for c odd
    rng = random.Random(1)
    for _ in range(400):
        c = rng.randrange(5, 10 ** 6) | 1
        if c % 3 == 0:
            continue
        a = rng.randrange(1, c)
        if gcd(a, c) != 1:
            continue
        D2 = 2 * dedekind6(a, c)
        assert (D2 - a - pow(a, -1, c)) % c == 0
        assert (D2 - (c + 1 - 2 * jacobi(a, c))) % 8 == 0


def test_dedekind_unit_decomposition():
    for p, q in [(5, 7), (7, 11), (11, 13), (13, 29)]:
        N = p * q
        for a in range(1, N):
            if gcd(a, N) != 1:
                continue
            sp, sq, u = diag.dedekind_split(a, p, q)
            assert u == diag.unit_part_brute(a, N)
            assert Fraction(dedekind6(a, N), 6 * N) == sp + sq + u


def test_cf():
    assert cf(3, 7) == [2, 3]
    assert cf(1, 5) == [5]


def test_modsym_hecke_relations():
    tbl, k, eigs = setup(11)
    assert k == 2 and eigs == [-2, 3]
    rng = random.Random(2)

    def br(x):
        return modsym_plus(x.numerator, x.denominator)
    for _ in range(150):
        N = rng.randrange(2, 400)
        a = rng.randrange(1, N)
        if gcd(a, N) != 1:
            continue
        x = Fraction(a, N)
        for l, al in ((3, -1), (5, 1), (7, -2), (13, 4)):
            assert sum(br((x + j) / l) for j in range(l)) + br(l * x) == al * br(x)


def test_salie_closed_form_small():
    import cmath
    import math

    def e(t):
        return cmath.exp(2j * math.pi * t)
    for N, (p, q) in [(15, (3, 5)), (21, (3, 7)), (35, (5, 7)), (77, (7, 11)), (91, (7, 13))]:
        for a, b in [(1, 1), (2, 3), (4, 1), (9, 2)]:
            if gcd(a * b, N) != 1:
                continue
            S = sum(jacobi(x, N) * e((a * x + b * pow(x, -1, N)) / N) for x in range(1, N) if gcd(x, N) == 1)
            roots = diag.salie_exact_terms(a, b, N, p, q)
            epsN = 1 if N % 4 == 1 else 1j
            cfm = epsN * math.sqrt(N) * jacobi(b, N) * sum(e(2 * v / N) for v in roots)
            assert abs(S - cfm) < 1e-8


def test_salie_reduction_factors():
    from src.common.instances import balanced_semiprime
    ok = 0
    for seed in range(8):
        inst = balanced_semiprime(48, 24_500_000_000_000 + seed)
        p, q = inst.labels["factors"]
        orc = diag.SalieOracle(p, q, delta=1e-4, seed=seed)
        g = diag.salie_factor(p * q, orc, c=3 + seed, delta=1e-4)
        ok += g in (p, q)
    assert ok >= 7


def test_theorem_RA_cubic_symbol_depends_on_N_mod_norm():
    from sympy import primerange
    rng = random.Random(3)
    for (u, v) in [(1, 3), (2, 3), (4, 3)]:
        nrm = u * u - u * v + v * v
        ps = [x for x in primerange(5, 1500) if nrm % x]
        tab, vals = {}, set()
        for _ in range(1500):
            a, b = rng.sample(ps, 2)
            k = diag.cubic_symbol_N(u, v, a, b)
            key = (a * b) % nrm
            assert tab.setdefault(key, k) == k
            vals.add(k)
        assert len(vals) == 3      # non-trivial symbol, yet a function of N mod Norm(alpha)


def test_battery_menu_and_isolation():
    from src.newfamily.battery import MENU_NAMES, compute_stats
    assert len(MENU_NAMES) == len(set(MENU_NAMES)) == 98
    st, aux = compute_stats(1000000016000000063)   # (10^9+7)(10^9+9)
    assert set(st) == set(MENU_NAMES)
    forbidden = {"factorint", "factors", "labels", "balanced_semiprime", "instances", "diag", "harness", "sqrt_mod"}
    for fn in ("recip.py", "modsym.py", "battery.py"):
        tree = ast.parse(open(os.path.join(ROOT, "src", "newfamily", fn)).read())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert not any(f in (node.module or "") for f in forbidden), (fn, node.module)
                assert not any(a.name in forbidden for a in node.names), fn
            if isinstance(node, ast.Attribute):
                assert node.attr not in forbidden, (fn, node.attr)
