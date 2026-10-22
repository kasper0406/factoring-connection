"""Tests for src/selfred (phase 22)."""
import ast
import random
from math import gcd, prod
from pathlib import Path

import pytest
import sympy

from src.common.instances import balanced_semiprime
from src.selfred import census as C
from src.selfred import harness as H
from src.selfred import reductions as R
from src.selfred import rigidity as G

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("mod", ["rigidity", "reductions"])
def test_operational_isolation(mod):
    tree = ast.parse((ROOT / f"src/selfred/{mod}.py").read_text())
    names = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            names |= {a.name for a in n.names}
        elif isinstance(n, ast.ImportFrom):
            names.add(n.module or "")
    assert not any(x.startswith(("sympy", "flint")) or "harness" in x or "census" in x for x in names), names


# ----------------------------------------------------------------------------- integer roots
def test_integer_roots_random():
    rng = random.Random(1)
    for _ in range(200):
        roots = [rng.randint(-10 ** 6, 10 ** 6) for _ in range(rng.randint(1, 5))]
        P = [rng.choice([1, -1, 2, 3])]
        for r in roots:
            P = G.pmul(P, [-r, 1])
        for _ in range(rng.randint(0, 2)):  # irreducible quadratic factors, no integer roots
            a = rng.randint(1, 10 ** 5)
            P = G.pmul(P, [a, 0, 1])
        if rng.random() < 0.3:  # repeated root
            P = G.pmul(P, [-roots[0], 1])
        lo, hi = -2 * 10 ** 6, 2 * 10 ** 6
        assert G.integer_roots(P, lo, hi) == sorted(set(roots))


def test_integer_roots_critic_regression():
    """P = 2(x+1)(7x^2+3x-1): P' has two roots in (-1, 0); v1 bracketing missed x = -1."""
    assert G.integer_roots([-2, 4, 20, 14], -40, 40) == [-1]


def test_integer_roots_bruteforce_fuzz():
    rng = random.Random(7)
    for _ in range(3000):
        P = G.trim([rng.randint(-30, 30) for _ in range(rng.randint(2, 6))])
        if len(P) < 2:
            continue
        want = [x for x in range(-60, 61) if G.peval(P, x) == 0]
        assert G.integer_roots(P, -60, 60) == want, P


def test_solve_rigid_critic_case():
    p, q = 2166859425403, 1887537611287
    P = [-2, 4, 20, 14]
    # A(x) = P(x - p - 1) + 2, B = 1, N2 = 2: valid promise N2 = A(p) B(q)
    sh = [-(p + 1), 1]
    A = [0]
    for c in reversed(P):
        A = G.psub(G.pmul(A, sh), [-c]) if A != [0] else [c]
    A = G.trim(A)
    A[0] += 2
    assert G.peval(A, p) == 2
    assert G.solve_rigid(p * q, 2, A, [1]) in (p, q)


def test_codex2_counterexamples():
    assert G.integer_roots([0, 1, 3, 2], -5, 5) == [-1, 0]          # x(x+1)(2x+1)
    assert G.integer_roots(G.pmul([0, 1, 3, 2], [0, 1, 3, 2]), -5, 5) == [-1, 0]
    assert G.solve_rigid(35, 6, [-324, 181, -33, 2], [1]) in (5, 7)
    assert C.unit_fixed_divisor([2 ** 20, 2 ** 20]) == 2 ** 21
    d, _ = R.factor_via_phi_mod_ell(3027, H.phi_mod, balanced=False)
    assert d in (3, 1009)


def test_integer_roots_close_roots_and_bigcoeffs():
    P = G.pmul(G.pmul([-(10 ** 30), 1], [-(10 ** 30 + 1), 1]), [-(10 ** 30 + 3), 7])
    assert G.integer_roots(P, 0, 10 ** 31) == [10 ** 30, 10 ** 30 + 1]


# ----------------------------------------------------------------------------- Theorem SR2
@pytest.mark.parametrize("seed", range(12))
def test_rigidity_factors(seed):
    rng = random.Random(seed)
    inst = balanced_semiprime(64 + 32 * (seed % 3), seed)
    p, q = inst.factors
    N = inst.N
    a, b = rng.randint(1, 4), rng.randint(1, 4)
    A = [rng.randint(-2 ** 40, 2 ** 40) for _ in range(a + 1)]
    B = [rng.randint(-2 ** 40, 2 ** 40) for _ in range(b + 1)]
    A[-1] = A[-1] or 1
    B[-1] = B[-1] or 1
    if seed % 4 == 0:   # torus-type map (p+1)(q+1), and p^2+p+1 Lucas-type
        A, B = [1, 1], [1, 1]
    if seed % 4 == 1:
        A, B = [1, 1, 1], [1, 1, 1]
    if seed % 2:        # attach to either factor
        p, q = q, p
    N2 = H.rigid_map_output(p, q, A, B)
    d = G.solve_rigid(N, N2, A, B)
    assert d is not None and N % d == 0 and 1 < d < N


def test_rigidity_degenerate_multiplier_class():
    inst = balanced_semiprime(64, 7)
    p, q = inst.factors
    for j in (1, 2, 3):
        A, B = [0] * j + [5], [0] * j + [-3]
        assert G.is_degenerate(inst.N, A, B)
        assert G.solve_rigid(inst.N, H.rigid_map_output(p, q, A, B), A, B) is None
    assert not G.is_degenerate(inst.N, [0, 1], [0, 0, 1])   # p q^2: not constant in x
    assert G.solve_rigid(inst.N, p * q * q, [0, 1], [0, 0, 1]) in (p, q)


# ----------------------------------------------------------------------------- reductions
def test_rad_decomposition():
    rng = random.Random(3)
    for _ in range(30):
        fac = {int(sympy.randprime(2, 10 ** 6)): rng.randint(1, 5) for _ in range(rng.randint(1, 4))}
        N = prod(p ** e for p, e in fac.items())
        dec, calls = R.squarefree_decomposition_via_rad(N, H.rad)
        want = {}
        for p, e in fac.items():
            want[e] = want.get(e, 1) * p
        assert dec == want
        assert calls == max(fac.values())


def test_valuation_split():
    inst = balanced_semiprime(64, 11)
    p, q = inst.factors
    N = inst.N
    assert R.valuation_split(p ** 3 * q ** 2 * 17, N)[0] == "split"
    kind, a, Y = R.valuation_split(N ** 4 * 1234567, N)
    assert (kind, a, Y) == ("balanced", 4, 1234567)
    assert R.valuation_split(1234567, N) == ("balanced", 0, 1234567)


@pytest.mark.parametrize("seed", range(6))
def test_phi_mod_ell(seed):
    inst = balanced_semiprime(80, seed)
    d, calls = R.factor_via_phi_mod_ell(inst.N, H.phi_mod)
    assert d in inst.factors
    assert calls <= 40


def test_euclid_discriminant_split():
    """f with a double root mod p only: Euclid on (f, f') over Z/N must hit a zero divisor (Lemma EUC)."""
    inst = balanced_semiprime(64, 5)
    p, q = inst.factors
    N = inst.N
    r = 12345
    # f = (x - r)(x - r - p) has disc p^2: double root mod p, simple mod q
    f = G.pmul([-r, 1], [-(r + p), 1])
    kind, d = R.poly_gcd_mod_N(f, G.pderiv(f), N)
    assert kind == "split" and d in (p, q)


# ----------------------------------------------------------------------------- census sanity
def test_e1_small_bound_holds():
    rows = C.e1_census(10)
    assert all(r["P1_ok"] for r in rows)
    assert rows[-1]["j"] == 20


def test_unit_fixed_divisor_and_sr1_small():
    h = G.pmul([0, 1], G.pmul([-1, 1], G.pmul([-2, 1], [-3, 1])))   # x(x-1)(x-2)(x-3)
    h = [6 * c for c in h]
    U = C.unit_fixed_divisor(h)
    # on units: 2-adic: x odd -> (x-1)(x-3) has v2 >= 3 ; 3-adic: x unit -> one of x-1, x-2 divisible by 3
    assert U == 6 * 8 * 3
    for p in (1009, 1013, 2003):
        r = C.e2_case(h, p)
        assert r["P1_ok"] and r["Gr"] == U


def test_weighted_null_fft_matches_exact():
    for o, d in [(300, 40), (57, 0), (0, 31), (5, 5)]:
        assert abs(C.weighted_null_excess(o, d) - C.weighted_null_excess(o, d, exact_limit=0)) < 1e-9
