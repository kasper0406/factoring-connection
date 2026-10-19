"""Tests for src/specialfam (agent specialfam, phase 7)."""
from __future__ import annotations

import ast
import pathlib
import random
from fractions import Fraction
from math import gcd, isqrt, prod

import pytest
import sympy

from src.specialfam.polyprod import block_products, norm_mod
from src.specialfam.normbox import (normbox_factor, polybox_search, disc_box_factor,
                                    near_square_factor, cubic_form_factor, ternary_factor,
                                    band_box)
from src.specialfam.smoothcopp import (smooth_part_ps, unit_roots, divisors, smoothcopp_factor,
                                       common_factor_factor, balanced_band)
from src.specialfam.copper import best_params, rigorous_ok
from src.specialfam import harness as H

ROOT = pathlib.Path(__file__).resolve().parents[1]
OPERATIONAL = ["polyprod.py", "normbox.py", "copper.py", "smoothcopp.py"]


def test_operational_modules_do_not_import_harness_or_factoring():
    banned = {"sympy", "src.specialfam.harness", "src.common.instances"}
    for fn in OPERATIONAL:
        tree = ast.parse((ROOT / "src" / "specialfam" / fn).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert node.module not in banned, (fn, node.module)
            if isinstance(node, ast.Import):
                for a in node.names:
                    assert a.name not in banned, (fn, a.name)
        src = (ROOT / "src" / "specialfam" / fn).read_text()
        for bad in ("factorint", "is_prime(", "isprime", ".factor(", "inverse_mod("):
            assert bad not in src, (fn, bad)


@pytest.mark.parametrize("seed", range(6))
def test_block_products_match_direct(seed):
    rng = random.Random(seed)
    N = rng.getrandbits(70) | 1
    S = [rng.randrange(N) for _ in range(rng.randint(1, 9))]
    chis = []
    for _ in range(rng.randint(1, 7)):
        d = rng.randint(1, 4)
        chis.append([rng.randrange(N) for _ in range(d)] + [1])
    got = block_products(N, S, chis)
    for c, g in zip(chis, got):
        want = prod(sum(ci * pow(s, i, N) for i, ci in enumerate(c)) for s in S) % N
        assert g == want


def test_norm_mod_is_resultant():
    rng = random.Random(7)
    y = sympy.symbols("y")
    for _ in range(20):
        d = rng.randint(1, 5)
        chi = [rng.randint(-50, 50) for _ in range(d)] + [1]
        r = [rng.randint(-50, 50) for _ in range(d)]
        N = 10 ** 12 + 39
        want = sympy.resultant(sum(c * y ** i for i, c in enumerate(chi)),
                               sum(c * y ** i for i, c in enumerate(r)) + 0 * y, y)
        assert norm_mod([c % N for c in chi], [c % N for c in r], N) == int(want) % N


def test_gaussian_ring_mechanism_equals_resultant_packaging():
    """r18's (Z/N)[i] evaluation of F(mv+ib) has norm equal to our block product."""
    N = 1000003 * 999983
    m = 7
    F = [1]
    for j in range(m):  # F(X) = prod (X + j), integer coefficients
        F = [((F[i - 1] if i else 0) + j * (F[i] if i < len(F) else 0)) for i in range(len(F) + 1)]
    for v in range(3):
        for b in range(1, 4):
            z = (m * v % N, b % N)
            acc = (0, 0)
            for c in reversed(F):  # Horner in (Z/N)[i]
                acc = ((acc[0] * z[0] - acc[1] * z[1] + c) % N, (acc[0] * z[1] + acc[1] * z[0]) % N)
            nrm = (acc[0] ** 2 + acc[1] ** 2) % N
            bp = block_products(N, list(range(m)), [[(m * v) ** 2 + b * b, 2 * m * v, 1]])[0]
            assert nrm == bp == prod((m * v + j) ** 2 + b * b for j in range(m)) % N


@pytest.mark.parametrize("seed", range(5))
def test_normbox_factor_on_family(seed):
    eta = Fraction(119, 1000)
    inst = H.normbox_instance(40, Fraction(1, 2), eta, seed)
    res = normbox_factor(inst["N"], eta, kmin=8)
    f = res["factor"]
    assert f and inst["N"] % f == 0 and 1 < f < inst["N"]
    assert res["band"] <= inst["kp"]


def test_normbox_gcd_N_scan():
    """Two promised primes in the same block force gcd = N; the scan recovers a factor."""
    eta = Fraction(1, 4)
    # p = 101^2 + 1, q = 105^2 + 4 ... build N = p*q with both in one band/block
    cands = [(a, b) for a in range(100, 140) for b in (1, 2) if sympy.isprime(a * a + b * b)]
    (a1, b1), (a2, b2) = cands[0], [c for c in cands if c[1] == cands[0][1]][1]
    N = (a1 * a1 + b1 * b1) * (a2 * a2 + b2 * b2)
    res = polybox_search(N, lambda t, b: [t * t + b * b, 2 * t, 1], lambda a, b: a * a + b * b,
                         0, 141, [b1], m=141)
    assert res["factor"] in (a1 * a1 + b1 * b1, a2 * a2 + b2 * b2)
    assert res["how"] == "scan"


def test_controls_never_return_wrong_factor():
    eta = Fraction(119, 1000)
    for seed in range(4):
        inst = H.control_instance(36, seed)
        res = normbox_factor(inst["N"], eta, kmin=8)
        f = res["factor"]
        if f:
            assert inst["N"] % f == 0 and 1 < f < inst["N"]


def test_other_packaged_families():
    rng = random.Random(3)
    # HKL-type box near the diagonal
    while True:
        a, b = 1000 + rng.randint(0, 30), 1000 + rng.randint(0, 30)
        if sympy.isprime(a * a + b * b):
            break
    q = sympy.nextprime(3 * 10 ** 6)
    N = (a * a + b * b) * q
    assert disc_box_factor(N, 1000, 1000, 32)["factor"] in (a * a + b * b, q)
    # p = a^2 + k (not a norm form)
    while True:
        a, k = rng.randint(500, 900), rng.randint(1, 40)
        if sympy.isprime(a * a + k):
            break
    N = (a * a + k) * int(sympy.nextprime(10 ** 6))
    assert near_square_factor(N, 900, 40)["factor"] in (a * a + k, N // (a * a + k))
    # x^3 + 2 y^3
    while True:
        x, y = rng.randint(1, 60), rng.randint(1, 60)
        if sympy.isprime(x ** 3 + 2 * y ** 3):
            break
    N = (x ** 3 + 2 * y ** 3) * int(sympy.nextprime(10 ** 7))
    assert cubic_form_factor(N, 60, 60)["factor"] in (x ** 3 + 2 * y ** 3, N // (x ** 3 + 2 * y ** 3))
    # ternary x^2 + y^2 + z^2 with long x side
    while True:
        x, y, z = rng.randint(100, 400), rng.randint(0, 3), rng.randint(1, 3)
        if sympy.isprime(x * x + y * y + z * z):
            break
    pp = x * x + y * y + z * z
    N = pp * int(sympy.nextprime(10 ** 6))
    assert ternary_factor(N, 400, 3, 3)["factor"] in (pp, N // pp)


def test_smooth_part_ps():
    rng = random.Random(11)
    for _ in range(10):
        M = rng.getrandbits(60) + 2
        B = rng.randint(2, 3000)
        want = {ell: e for ell, e in sympy.factorint(M).items() if ell <= B}
        assert smooth_part_ps(M, B) == want


def test_unit_roots_complete():
    for fac in ({2: 1}, {2: 2}, {2: 3}, {2: 5, 3: 1}, {3: 2, 5: 1, 7: 1}, {2: 4, 3: 1, 11: 1}):
        d = prod(l ** e for l, e in fac.items())
        brute = sorted(r for r in range(d) if r * r % d == 1)
        assert unit_roots(fac) == brute


def test_root_count_bound():
    """Sum over d | H of #{r: r^2 = 1 mod d} <= 2 tau(H)^2 (Lemma R)."""
    for fac in ({2: 7, 3: 3, 5: 2, 7: 1, 11: 1}, {2: 2, 3: 1, 13: 1, 17: 1, 19: 1, 23: 1}):
        tot = sum(len(unit_roots(f)) for d, f in divisors(fac))
        tau = prod(e + 1 for e in fac.values())
        assert tot <= 2 * tau * tau


def test_coppersmith_params_rigorous():
    N = (2 ** 40 + 15) * (2 ** 40 + 75)
    D, m, Y = best_params(N, 2 ** 39, 9)
    assert rigorous_ok(N, 2 ** 39, Y, D, m) and not rigorous_ok(N, 2 ** 39, Y + 1, D, m)


@pytest.mark.parametrize("seed", range(3))
def test_smoothcopp_mixed_signs(seed):
    inst = H.smooth_gcd_instance(56, {2: 4, 3: 2, 5: 1, 7: 1}, "+-+-", "-+--", seed)
    res = smoothcopp_factor(inst["N"], 16, inst["g"], balanced_band(inst["N"]))
    assert res["factor"] in (inst["p"], inst["q"])


@pytest.mark.parametrize("seed", range(3))
def test_common_factor_general_g(seed):
    inst = H.common_gcd_instance(56, 10, seed)
    res = common_factor_factor(inst["N"], balanced_band(inst["N"]))
    assert res["factor"] in (inst["p"], inst["q"])


def test_band_box_monotone():
    eta = Fraction(119, 1000)
    for k in range(8, 60):
        A, B = band_box(k, eta)
        assert A * A <= 2 ** k < (A + 1) ** 2
        assert B ** 1000 <= 2 ** (k * 119) < (B + 1) ** 1000
