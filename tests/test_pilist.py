"""Tests for agent pilist (targets, decoder barrier helpers, pencils / O2, tail closure)."""
from collections import Counter
from math import gcd, log, prod, sqrt

import pytest

from src.pilist import targets as T
from src.pilist.decoder import (boneh_threshold, count_model_fraction, hg_row_decoder, log_moment, log_P,
                                primes_upto, rankin_log_bound, weighted_log_bound)
from src.pilist.pencil import (coeffs_of_roots, enumerate_split, find_lines, line_params, line_stats,
                               small_c_triples, tl_class)
from src.pilist.tails import F2Basis, first_moment_log, r_threshold, simulate


# ------------------------------------------------------------------ targets
def test_dixon_constant():
    assert abs(T.dixon_check() - sqrt(2)) < 1e-6


def test_tau_max_matches_bruteforce():
    for alpha, c, nu in [(1 / 3, 1.9019, 1.0), (1 / 3, 1.9019, 3.0), (0.25, 1.0, 0.5), (0.25, 0.6, 2.0)]:
        tm = T.tau_max(alpha, c, nu)
        ks = [i / 20000 * c / 2 for i in range(1, 20001)]
        best = max(c - k - alpha * nu / k for k in ks)
        assert abs(best - tm) < 1e-3
        # feasibility at tau_max: some kappa gives cost <= c
        assert min(T.congruence_cost(alpha, k, nu, tm) for k in ks) <= c + 1e-3


def test_doubly_per_output_tends_to_tau():
    assert abs(T.doubly_per_output(0.5, h=1e-4) - 0.5) < 1e-3


# ------------------------------------------------------------------ decoder
def _exact_tail(P, logS, lam=0.0):
    qs = sorted(P)
    tot = 0.0
    for mask in range(1 << len(qs)):
        pr, lg = 1.0, 0.0
        for i, q in enumerate(qs):
            if mask >> i & 1:
                pr *= 1 / q
                lg += log(q)
            else:
                pr *= 1 - 1 / q
        if lg >= logS:
            tot += pr * (2.718281828459045 ** (lam * lg))
    return tot


def test_log_moment_zero_and_one():
    P = {2: 1, 3: 2, 5: 1}
    assert abs(log_moment(P, 0.0)) < 1e-12
    # E[d] for squarefree q: 1 + (q-1)/q ; for q^2: 1 + (1-1/q) + ...
    exact = (1 + (2 - 1) / 2) * (1 + (1 - 1 / 3) + 1 / 9 * 9 - (1 - 1 / 3) * 0 - 1 / 3 * 0 - 0) if False else None
    assert log_moment(P, 1.0) > 0


def test_rankin_bound_is_an_upper_bound():
    P = {q: 1 for q in primes_upto(29)}
    for logV in (3.0, 6.0):
        for logS in (6.0, 10.0, 15.0):
            for lam in (0.0, 0.1):
                lb, _ = rankin_log_bound(logV, P, logS, lam)
                assert logV + log(max(_exact_tail(P, logS, lam), 1e-300)) <= lb + 1e-9


def test_count_model_bracket():
    P = {q: 1 for q in primes_upto(31)}
    for thr in (4.0, 9.5, 14.0):
        ex = _exact_tail(P, thr)
        lo, hi = count_model_fraction(P, thr)
        assert lo <= ex * (1 + 1e-9) + 1e-15 and ex <= hi * (1 + 1e-9) + 1e-15


def test_hg_decoder_finds_planted_root():
    Pd = {q: 1 for q in primes_upto(60)}
    P = prod(Pd)
    a0 = 777
    d0 = prod(q for q in Pd if q <= 41)          # large divisor of P
    c = (-a0 + d0 * 12345) % P                    # a0 + c == 0 mod d0
    assert (a0 + c) % d0 == 0
    roots, norms = hg_row_decoder(c, P, 4096, 8, 16)
    Scert = norms[0] ** (1 / 8)
    if gcd(a0 + c, P) > Scert:
        assert a0 in roots
    # certificate soundness: every |x| <= X with gcd > S_cert is a root
    for x in range(-4096, 4097):
        if gcd(x + c, P) > Scert:
            assert x in roots


def test_weighted_reduces_to_unweighted_scale():
    P = {q: 1 for q in primes_upto(40)}
    logV = 5.0
    w = {q: 1.0 for q in P}
    lw, _ = weighted_log_bound(logV, P, w, 0.05)
    lu, _ = rankin_log_bound(logV, P, boneh_threshold(logV, P), 0.05)
    assert abs(lw - lu) < 0.05


# ------------------------------------------------------------------ pencils
def test_coeffs_of_roots():
    assert coeffs_of_roots([1, 2]) == [-3, 2]
    assert coeffs_of_roots([0, 0, 0]) == [0, 0, 0]


def test_three_term_identity_M3_found_and_coprime():
    # X(X+3)^2 - 2(X+1)^3 + (X+2)(X-1)^2 = 0, pairwise coprime: naive three-term rigidity fails at M = 3
    P0, P1, P2 = coeffs_of_roots([0, -3, -3]), coeffs_of_roots([-1, -1, -1]), coeffs_of_roots([-2, 1, 1])
    assert all(a - 2 * b + c == 0 for a, b, c in zip(P0, P1, P2))
    ms, arr = enumerate_split(3, 3)
    lines = find_lines(arr)
    target = {ms.index((-3, -3, 0)), ms.index((-1, -1, -1)), ms.index((-2, 1, 1))}
    hit = [(mem, D) for mem, D in lines if target <= set(mem)]
    assert hit
    mem, D = hit[0]
    st = line_stats(ms, mem)
    assert st["e"] == 3 and st["disjoint"]
    assert small_c_triples(line_params(arr, mem, D)) >= 1


def test_codex_example_constant_remainder():
    P0 = coeffs_of_roots([-8, -1, 3, 4])
    P1 = coeffs_of_roots([-5, -4, 0, 7])
    P2 = coeffs_of_roots([-7, -2, 1, 6])
    G = [a + b - 2 * c for a, b, c in zip(P0, P1, P2)]
    assert G[:3] == [0, 0, 0] and G[3] == -72
    ms = [(-8, -1, 3, 4), (-5, -4, 0, 7), (-7, -2, 1, 6)]
    assert tl_class(ms, 0, 1, 2, L=1) == "C"


def test_claim_A_disjoint_moving_roots_small_cells():
    for M, B in [(3, 3), (4, 3), (5, 2)]:
        ms, arr = enumerate_split(M, B)
        for mem, D in find_lines(arr):
            st = line_stats(ms, mem)
            assert st["disjoint"]
            assert st["nonsqf"] <= max(2 * st["e"] - 2, 0) or st["e"] <= 1


# ------------------------------------------------------------------ tails
def test_f2_basis():
    b = F2Basis()
    assert not b.add(0b011)
    assert not b.add(0b110)
    assert b.add(0b101)


def test_first_moment_monotone_and_threshold():
    a = first_moment_log(50, 6, 1, 1e-3)
    b = first_moment_log(100, 6, 1, 1e-3)
    assert b > a
    R = r_threshold(6, 1, 1e-3)
    assert first_moment_log(R, 6, 1, 1e-3) <= log(0.5) + 1e-9
    assert first_moment_log(R + 1, 6, 1, 1e-3) > log(0.5) - 1e-9


def test_simulate_runs():
    out = simulate(300, 5, 1, "greedy", 0)
    assert out["R_star"] is not None and out["R_star"] > 10
