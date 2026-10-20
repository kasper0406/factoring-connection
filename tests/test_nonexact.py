"""Tests for src/nonexact (phase 17, O1)."""
import random
from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb

import numpy as np
import pytest
from flint import fmpz_poly

from src.nonexact import core as C


def brute_box_points(w, N, B, s, V):
    bnd = C.e_bounds(B, s)
    cnt = 0
    import itertools
    for e in itertools.product(*[range(-b, b + 1) for b in bnd]):
        val = (w[0] + sum(e[i] * w[i + 1] for i in range(s))) % N
        if abs(C.centred(val, N)) <= V:
            cnt += 1
    return cnt


def test_lin_weights_match_evaluation():
    rng = random.Random(1)
    N = 1000003 * 7
    for _ in range(50):
        s = rng.randint(1, 4)
        S = [rng.randint(-9, 9) for _ in range(s)]
        mu, c = rng.randrange(N), rng.randrange(1, N)
        w = C.lin_weights(c, mu, N, s)
        e = [(-1) ** i * C.coeffs(C.from_roots(S), s + 1)[s - i] for i in range(1, s + 1)]
        lhs = (c * int(C.from_roots(S)(mu))) % N
        rhs = (w[0] + sum(e[i] * w[i + 1] for i in range(s))) % N
        assert lhs == rhs


@pytest.mark.parametrize("s,B,V", [(1, 5, 3), (2, 4, 7), (2, 3, 40), (3, 2, 5)])
def test_lin_box_points_vs_brute(s, B, V):
    rng = random.Random(s * 100 + B)
    N = 10007 * 101
    for _ in range(5):
        mu, c = rng.randrange(N), rng.randrange(1, N)
        w = C.lin_weights(c, mu, N, s)
        if np.gcd(w[s], N) != 1:
            continue
        assert C.lin_box_points(w, N, B, s, V) == brute_box_points(w, N, B, s, V)


def test_lin_split_hits_vs_direct():
    rng = random.Random(5)
    N = 65537 * 257
    B, s, V = 6, 3, 2000
    ms = C.s_multisets_e(B, s)
    E = np.array([e for _, e in ms], dtype=np.int64)
    for _ in range(20):
        mu, c = rng.randrange(N), rng.randrange(1, N)
        w = C.lin_weights(c, mu, N, s)
        direct = sum(1 for S, _ in ms if abs(C.centred(c * int(C.from_roots(S)(mu)), N)) <= V)
        assert C.lin_split_hits(w, N, V, E) == direct


@pytest.mark.parametrize("s,B,V", [(1, 6, 50), (2, 4, 30), (3, 2, 20)])
def test_lattice_enumeration_finds_all_box_points(s, B, V):
    rng = random.Random(11 * s)
    N = 4099 * 4111
    for _ in range(4):
        mu, c = rng.randrange(N), rng.randrange(1, N)
        w = C.lin_weights(c, mu, N, s)
        cnt = {}
        pts = C.lin_lattice_enumerate(w, N, B, s, V, cnt)
        assert len(pts) == brute_box_points(w, N, B, s, V)
        for p in pts:
            val = (w[0] + sum(p[i] * w[i + 1] for i in range(s))) % N
            assert C.centred(val - p[s], N) == 0 and abs(p[s]) <= V


def test_multiplier_lattice_index_and_membership():
    rng = random.Random(3)
    N = 1009 * 1013
    mu, mu2 = rng.randrange(N), rng.randrange(N)
    F = C.quadratic_presentation(mu, mu2, N)
    assert int(F(mu)) % N == 0 and int(F(mu2)) % N == 0
    for (j, k) in [(1, 1), (2, 2), (4, 2)]:
        A = [rng.randint(-5, 5) for _ in range(6)]
        P = C.from_roots(A)
        Lb = C.multiplier_lattice(P, F, N, j, k)
        assert C.det_abs(Lb) == N ** 2  # surjective case (unit P, j + 1 + k >= 2)
        for v in Lb:
            Q = fmpz_poly(v[: j + 1])
            r = fmpz_poly(v[j + 1:])
            R = C.coeffs(C.red_mod(Q * P - r, F), 2)
            assert all(x % N == 0 for x in R)


def test_split_detection():
    assert C.split_over_Z_small(C.coeffs(3 * C.from_roots([1, -2, 2])), 2)
    assert not C.split_over_Z_small(C.coeffs(C.from_roots([1, 5])), 4)
    assert not C.split_over_Z_small([1, 0, 1], 10)
    assert C.split_over_Z_small([7], 1)


def test_ih_codex_example_k1():
    # ptebarrier's example: P - P' = 14400, N = 25, r = r' = -1 => one relation, u = 14400 == 0 = v? (delta = 0)
    A = [0, 0, 7, 7, -7, -7]
    A2 = [8, -8, 5, -5, 3, -3]
    rels, E = C.ih_relations(A, A2, [-1], [-1], 25, 1)
    vals = {(abs(u), abs(v)) for _, u, v, nt in rels if nt}
    assert vals == {(14400, 0)}
    for _, u, v, nt in rels:
        assert (u - v) % 25 == 0


def _k2_pairs(M, R, limit=50):
    groups = {}
    for A in combinations_with_replacement(range(-R, R + 1), M):
        key = (sum(A), sum(a * a for a in A)) + tuple(sum(a ** t for a in A) for t in range(3, M - 1))
        groups.setdefault(key, []).append(A)
    out = []
    for g in groups.values():
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                D = C.from_roots(g[i]) - C.from_roots(g[j])
                if D.degree() == 1:
                    out.append((g[i], g[j]))
                    if len(out) >= limit:
                        return out
    return out


def test_ih_k2_genuine_degree_one_difference():
    from sympy import primefactors
    pairs = _k2_pairs(4, 6)
    assert pairs
    tested = 0
    for A, A2 in pairs:
        c = C.coeffs(C.from_roots(A) - C.from_roots(A2), 2)
        for d1 in range(-2, 3):
            for d0 in range(-2, 3):
                g = int(np.gcd(c[1] - d1, c[0] - d0))
                for N in (primefactors(g) if g else []):
                    if N <= 9:
                        continue
                    rels, E = C.ih_relations(A, A2, [d0, d1], [0, 0], N, 2)
                    nontriv = {(abs(u), abs(v)) for _, u, v, nt in rels if nt}
                    for _, u, v, nt in rels:
                        assert (u - v) % N == 0 and nt == (u != v)
                    n_supp = len(set(A) ^ set(A2))
                    assert len(nontriv) >= C.ih_distinct_lower_bound(4, 2, 1, n_supp)
                    tested += 1
    assert tested > 0


def test_runner_smoke(tmp_path, monkeypatch):
    """Tiny non-registered config: checks that E1-E3 run end to end (not an experiment)."""
    from src.nonexact import run_experiment as RE
    monkeypatch.setenv("NONEXACT_RAW", str(tmp_path))
    cfg = {
        "E1": {"cells": [{"n": 20, "V": 63, "B": [4], "s": [1, 2]}], "target_hits": 5, "T_max": 20000,
               "box_subset": 2, "enum_subset": 2, "enum_p_max": 50, "instance_seed": 9, "rng_seed": 1},
        "E2": {"n": [24], "d": [1, 2], "jk_d1": [[2, 1]], "jk_d2": [[1, 2]], "B": 6, "M": 5, "T": 3,
               "baseline_draws": 2, "instance_seed": 9, "rng_seed": 2},
        "E3": {"cells": [[1, 4, 8], [2, 4, 6]], "H": 3, "N_min": 50, "pairs_cap": 3, "r_tries": 30, "rng_seed": 3},
    }
    RE.run_E1(cfg)
    RE.run_E2(cfg)
    RE.run_E3(cfg)
    import json
    rows = [json.loads(l) for l in open(tmp_path / "nonexact_E3.jsonl")]
    assert all(r["violations"] == 0 for r in rows)
