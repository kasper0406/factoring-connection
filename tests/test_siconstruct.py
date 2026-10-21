"""Tests for agent `siconstruct` (phase 20)."""
import random

import mpmath as mp
import sympy as sp

from src.common.instances import balanced_semiprime
from src.siconstruct.constructors import aim2d_construct, brute_interval, mitm_construct
from src.siconstruct.core import (factor_over, is_smooth, lpath_params, primes_upto, random_smooth_list,
                                  sieve_window)
from src.siconstruct.theory import (cvp_bounds, linear_oracle_cost, mitm_expected_hits, si_cost, si_cost_slack,
                                    tau_threshold)


def _inst(bits=24, seed=7):
    N = balanced_semiprime(bits, seed).N
    pp = lpath_params(bits)
    return N, pp["y"], pp["h"], primes_upto(pp["y"])


def test_sieve_complete_and_exact():
    N, y, h, P = _inst()
    for X in (N, 3 * N + 17):
        got, cost = sieve_window(X, 6000, P, chunk=1000)
        ref = [z for z in range(X, X + 6000) if is_smooth(z, P)]
        assert got == ref and cost == 6000


def test_brute_outputs_valid():
    N, y, h, P = _inst()
    r = brute_interval(N, 2, h, P, 4000)
    for u in r["outputs"]:
        assert 2 * N <= u <= 2 * N + h and is_smooth(u, P)


def test_mitm_finds_exactly_the_window_pairs():
    N, y, h, P = _inst()
    rng = random.Random(3)
    A = random_smooth_list(1 << 11, 1 << 12, P, 80, rng)
    B = random_smooth_list(3 * N >> 11, 3 * N >> 10, P, 80, rng)
    r = mitm_construct(N, h, A, B)
    ref = sorted({a * b for a in A for b in B if a * b >= N and (a * b) % N <= h})
    assert [u for u, _ in r["outputs"]] == ref
    for u, (k, w) in r["outputs"]:
        assert u == k * N + w and 0 <= w <= h and is_smooth(u, P)


def test_aim2d_outputs_valid_and_cofactor_bound():
    N, y, h, P = _inst()
    rng = random.Random(4)
    S = random_smooth_list(1 << 8, 1 << 9, P, 10, rng)
    K = 4
    r = aim2d_construct(N, h, K, S, P, max_per_sk=50)
    for u, (k, w, s) in r["outputs"]:
        assert u == k * N + w and 0 <= w <= h and u % s == 0 and is_smooth(u, P)
    # Lemma RC: every tested cofactor is >= kN/s >= N / max(S)
    assert 2 ** r["cofactor_bits_min"] > N // max(S)


def test_hits_formula_exact_rational():
    assert mitm_expected_hits(10, 20, 99, 1000) == sp.Rational(10 * 20 * 100, 1000)


def test_si_cost_constants():
    r = si_cost(0)
    assert sp.simplify(r["c"] ** 3 - sp.Rational(16, 9)) == 0
    assert abs(float(tau_threshold(sp.Rational(19019, 10000))) - 0.757186) < 1e-5
    assert sp.simplify(linear_oracle_cost() - 6 ** sp.Rational(1, 3)) == 0
    # with tau below the threshold the constant beats MNFS 1.9019
    assert float(si_cost(sp.Rational(3, 4))["c"]) < 1.9019 < float(si_cost(sp.Rational(76, 100))["c"])
    assert float(si_cost_slack(0)["c"]) == float(si_cost(0)["c"])
    assert float(si_cost_slack(2)["c"]) > 1.89


def test_prop_cv_inequalities_on_witnesses():
    """Prop CV(a),(b): every lattice vector has dist^2 >= ln kN - 1/(4C^2); a squarefree smooth u in I_k has
    dist^2 <= ln kN + h/kN + 1.  Checked at 80 digits on witnesses found by the sieve and on random vectors."""
    N, y, h, P = _inst(28, 1)
    P = primes_upto(400)  # a larger base so that squarefree witnesses exist at toy size
    k = 3
    found, _ = sieve_window(k * N, h + 1, P)
    sqf = [u for u in found if all(e == 1 for e in factor_over(u, P).values())]
    assert sqf
    for u in sqf[:20]:
        f = factor_over(u, P)
        e = [f.get(p, 0) for p in P]
        b = cvp_bounds(e, P, N, k, h)
        assert b["lower"] <= b["dist2"] <= b["witness_upper"]
        assert b["dist2"] / b["lower"] - 1 <= mp.mpf(1.01) / mp.log(k * N)
    rng = random.Random(5)
    for _ in range(200):
        e = [rng.randint(-2, 2) for _ in P]
        b = cvp_bounds(e, P, N, k, h)
        assert b["dist2"] >= b["lower"]


def test_prop_cv_converse_small_v():
    """Prop CV(c): a vector with dist^2 <= ln kN + c' has v <= exp((c' + |Delta|)/2) and |w| <= v kN |Delta| e^|Delta|."""
    N, y, h, P = _inst(28, 2)
    k = 5
    rng = random.Random(9)
    for _ in range(300):
        e = [rng.randint(-1, 1) for _ in P]
        b = cvp_bounds(e, P, N, k, h)
        cprime = b["dist2"] - mp.log(k * N)
        u = 1
        v = 1
        for ei, p in zip(e, P):
            if ei > 0:
                u *= p ** ei
            elif ei < 0:
                v *= p ** (-ei)
        assert mp.log(v) <= (cprime + abs(b["Delta"])) / 2 + mp.mpf(10) ** -30
        w = u - v * k * N
        assert abs(w) <= v * k * N * abs(b["Delta"]) * mp.e ** abs(b["Delta"]) + 1
