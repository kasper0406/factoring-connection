"""Tests for src/sicheck (agent sicheck, phase 20)."""
import inspect
import random
from fractions import Fraction as Fr

import sympy as sp

from src.sicheck import core, model


def test_sieve_matches_brute_force():
    P = core.primes_upto(50)
    for A in (2, 1000, 10 ** 9 + 7, 2 ** 40 + 12345):
        H = 20000
        brute = [n for n in range(A, A + H) if core.is_smooth(n, P)]
        assert core.smooth_in_interval(A, H, P) == brute


def test_jacobi_matches_sympy():
    rng = random.Random(1)
    for _ in range(300):
        n = rng.randrange(3, 10 ** 6) | 1
        a = rng.randrange(0, 10 ** 6)
        assert core.jacobi(a, n) == sp.jacobi_symbol(a, n)


def test_model_exact_constants():
    P = model.si_program()
    assert P["kappa_star_cubed"] == sp.Rational(2, 9)
    assert P["c_cubed"] == sp.Rational(16, 9)
    assert P["kappa_min_cubed"] == sp.Rational(4, 9)
    assert sp.simplify(P["nu_threshold"] - 2 / (3 * model.kap)) == 0
    assert model.si_optimum_is_global()
    S = model.size_principle_consistency()
    assert S["b_opt_minus_kappa_star"] == 0 and S["c_ss_cubed"] == sp.Rational(16, 9)
    grid = [Fr(i, 60) for i in range(1, 30)]
    assert min(grid, key=model.single_interval_exponent) == Fr(1, 3)
    assert model.single_interval_exponent(Fr(1, 3)) == Fr(1, 3)
    E = model.single_interval_exponent_program()
    assert sp.simplify(E["residue_penalty_level"] - (1 - 2 * model.al)) == 0


def _blum_instance():
    # p, q = 3 mod 4 (harness-only knowledge)
    p, q = 1000003, 1000039
    assert p % 4 == 3 and q % 4 == 3
    return p * q


def test_factor_si_end_to_end_small():
    N = 1000003 * 1000033
    res = core.factor_si(N, 150, 1_100_000, constructor="rand", seed=0, excess=15)
    f = res["factor"]
    assert f is not None and 1 < f < N and N % f == 0
    assert res["stats"]["relations"] >= res["stats"]["fb_size"] + 15


def test_beta_formula_and_lemma_H_blum():
    N = _blum_instance()
    res = core.factor_si(N, 150, 1_100_000, constructor="next", seed=1, excess=15, n_random_deps=0)
    rels, deps, primes = res["relations"], res["deps"], res["primes"]
    beta = [core.beta_blum(r.u, r.v, primes, N) for r in rels]
    rng = random.Random(3)
    masks = list(deps.basis)
    for _ in range(30):
        m = 0
        for b in deps.basis:
            if rng.getrandbits(1):
                m ^= b
        if m:
            masks.append(m)
    for m in masks:
        Z, g = core.square_root_of_unity(rels, m, primes, N)
        assert g == 1
        fr = int(core.classify_root(Z, N)[0])
        pred = 0
        for i in range(len(rels)):
            if m >> i & 1:
                pred ^= beta[i]
        assert fr == pred


def test_qs_factor_base_degenerate_for_blum():
    p, q = 1000003, 1000099          # both 3 mod 8, so N = 1 mod 8 and J_N(2) = +1 (2 is in the J+ base)
    N = p * q
    assert core.jacobi(2, N) == 1
    primes = core.primes_upto(330)
    qfb = core.qs_factor_base(primes, N)
    res = core.factor_si(N, 330, 1_100_000, constructor="next", seed=0, excess=10, factor_base=qfb,
                         n_random_deps=20, k_max=6000)
    st = res["stats"]
    assert st["dim_deps"] >= 1
    assert not st["phi_nonzero"] and st["random_fruitful"] == 0


def test_operational_code_has_no_leakage():
    code = [ln for ln in inspect.getsource(core).splitlines()
            if ln.strip().startswith(("import", "from")) or "(" in ln and not ln.strip().startswith(("#", '"'))]
    imports = [ln for ln in code if ln.strip().startswith(("import", "from"))]
    assert all("diag" not in ln and "sympy" not in ln and "instances" not in ln for ln in imports)
    body = "\n".join(code)
    for bad in ("factorint", ".labels", "isprime", "totient", "_diag("):
        assert bad not in body, bad
