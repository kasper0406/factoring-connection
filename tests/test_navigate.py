"""Tests for src/navigate (agent `navigate`, phase 25)."""
import math
import random
from collections import Counter
from itertools import product

import numpy as np
import pytest

from src.common.instances import balanced_semiprime
from src.navigate import e2_counts as E2
from src.navigate.core import (SmoothTester, factor_over, gauss_split, primes_upto, reduced_key, run_g1, run_j2,
                               run_s0, verify_relation)


def test_smooth_tester_matches_trial_division():
    fb = primes_upto(97)
    T = SmoothTester(fb)
    rng = random.Random(1)
    for z in list(range(1, 3000)) + [rng.randrange(1, 1 << 40) for _ in range(3000)] + [2 ** 40 * 97 ** 3]:
        assert T.is_smooth(z) == (factor_over(z, fb) is not None)


def test_gauss_split_is_relation():
    N = balanced_semiprime(40, 7).N
    rng = random.Random(2)
    for _ in range(200):
        t = rng.randrange(1, N)
        u, v = gauss_split(t, N)
        assert (u - t * v) % N == 0 and u * u < N and v * v <= 4 * N


def test_e1_algorithms_output_valid_relations():
    N = balanced_semiprime(40, 11).N
    fb = primes_upto(1 << 10)
    rng = random.Random(3)
    for _ in range(5):
        r = rng.randrange(1, N)
        a = r * r % N
        for fn in (run_s0, run_g1, run_j2):
            R = fn(N, a, fb, "t", budget=200)
            assert R["tests"] <= 200
            for u, v in R["outputs"]:
                assert verify_relation(u, v, a, N, fb)


def test_verify_rejects_bad():
    N = balanced_semiprime(40, 11).N
    fb = primes_upto(1 << 10)
    assert not verify_relation(2, 3, 5, N, fb)          # wrong congruence
    a = 6 * pow(5, -1, N) % N
    assert verify_relation(6, 5, a, N, fb)
    assert not verify_relation(6 * 1031, 5 * 1031, a, N, fb)  # 1031 > 2^10 not in fb
    assert reduced_key(6, -4) == (-3, 2) and reduced_key(-3, 2) == (-3, 2)


def _brute_counts(N, p, q, fb, B):
    """C(t) by brute force over signed exponent vectors with weight <= B (independent of e2_counts)."""
    lns = [math.log(l) for l in fb]
    ranges = [range(-int(B / x), int(B / x) + 1) for x in lns]
    C = Counter()
    for e in product(*ranges):
        if sum(abs(k) * x for k, x in zip(e, lns)) <= B + 1e-9:
            t = 1
            for l, k in zip(fb, e):
                t = t * pow(l, k, N) % N
            C[t] += 1
    return C


def test_e2_counts_match_brute_force(monkeypatch):
    # tiny instance: 14-bit semiprime, fb = primes <= 7, budget 1.5 ln N
    seed = 5
    inst = balanced_semiprime(14, seed)
    monkeypatch.setattr(E2, "BITS", 14)
    res_cells = E2.run_instance(seed, 7, [1.0, 1.5])
    N, p, q = inst.N, inst.labels["p"], inst.labels["q"]
    fb = primes_upto(7)
    for beta, cell in zip([1.0, 1.5], res_cells["cells"]):
        C = _brute_counts(N, p, q, fb, beta * math.log(N))
        # H = <fb>: compute by closure
        H = {1}
        frontier = [1]
        while frontier:
            x = frontier.pop()
            for l in fb:
                for z in (x * l % N, x * pow(l, -1, N) % N):
                    if z not in H:
                        H.add(z)
                        frontier.append(z)
        assert res_cells["H"] == len(H)
        vals = np.array([C.get(t, 0) for t in H], dtype=float)
        assert abs(cell["mu"] - vals.mean()) < 1e-9
        assert abs(cell["h"] - (vals >= 1).mean()) < 1e-9
        assert cell["max"] == int(vals.max())


def test_lemma_RT_random_target_size_floor():
    # Lemma RT: #{t unit : exists u = t v (mod N), 0 < |u||v| <= X} <= 4 X (1 + ln X)
    for seed in range(3):
        N = balanced_semiprime(16, seed).N
        for X in (N // 64, N // 16, N // 4):
            hit = set()
            for v in range(1, X + 1):
                if math.gcd(v, N) != 1:
                    continue
                vinv = pow(v, -1, N)
                for u in range(1, X // v + 1):
                    for su in (u, -u):
                        hit.add(su * vinv % N)
            assert len(hit) <= 4 * X * (1 + math.log(X))
            assert len(hit) <= 2 * X * (1 + math.log(X))  # signs of v are absorbed: (u,v) ~ (-u,-v)


def test_prop_SQ_unconditional_direction():
    # a uniform unit is a square with probability exactly 2^-omega (omega = 2), for any fixed shift q
    inst = balanced_semiprime(16, 1)
    N, p, q = inst.N, inst.labels["p"], inst.labels["q"]
    units = [t for t in range(1, N) if math.gcd(t, N) == 1]
    sq = {t * t % N for t in units}
    assert len(sq) * 4 == len(units)
    for shift in (1, 2, 3, 5):
        assert sum(1 for t in units if t * shift % N in sq) * 4 == len(units)
