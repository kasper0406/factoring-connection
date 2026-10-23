"""Tests for src/smoothbuild (agent smoothbuild, phase 24)."""
import math
import random
from itertools import combinations

import flint
import pytest

from src.common.instances import balanced_semiprime
from src.smoothbuild.core import (FAMILIES, GROUPS, atom_product_count, burnside_rank, gauss_split,
                                  hidden_root_class_factor, inv_square_sum_fixed, left_kernel_rank, poly_eval,
                                  primes_upto, psi_count, root_list, trial_factor)


def test_burnside_ranks():
    expect = {"S4": (24, 2), "A4": (12, 2), "D4": (8, 3), "C4": (4, 4), "V4": (4, 4), "S6": (720, 2), "C6": (6, 6)}
    for name, (order, rank) in expect.items():
        o, orbits, r = burnside_rank(GROUPS[name], len(GROUPS[name][0]))
        assert (o, orbits, r) == (order, 1, rank)


def test_galois_rank_bound_regular_iff_max():
    # rank - 1 <= d - 1, equality iff |G| = d (regular action)
    for name, gens in GROUPS.items():
        d = len(gens[0])
        o, _, r = burnside_rank(gens, d)
        assert r - 1 <= d - 1
        assert (r == d) == (o == d)


def test_root_sieve_matches_factorisation():
    """c_p from roots mod p equals divisibility read off flint factorisations (pre-factored oracle)."""
    P1, P2 = 1 << 16, 1 << 18
    primes = [p for p in primes_upto(P2) if p > P1]
    for name in ("C4:Phi5", "S4:x4+x+1", "D4:x4-2"):
        coeffs, _ = FAMILIES[name]
        roots = {p: set(root_list(coeffs, p)) for p in primes}
        for a in range(70000, 70060):
            v = poly_eval(coeffs, a)
            fac = {int(q) for q, _e in flint.fmpz(v).factor() if P1 < int(q) <= P2}
            sieve = {p for p in primes if a % p in roots[p]}
            assert fac == sieve


def test_gauss_split_invariants():
    rng = random.Random(5)
    for bits in (40, 64, 100):
        N = balanced_semiprime(bits, 1).N
        for _ in range(200):
            t = rng.randrange(1, N)
            u, v = gauss_split(t, N)
            assert (u - t * v) % N == 0
            assert u * u < N
            assert v != 0 and v * v <= 4 * N


def test_hidden_root_factoring_structural():
    for bits in (32, 40, 48):
        inst = balanced_semiprime(bits, 3)
        res = hidden_root_class_factor(inst.N, 2 ** (bits // 4 + 1), 0)
        assert res["factor"] in inst.labels["distinct"]
        assert res["kind"] == "STRUCTURAL"


def test_lemma_AT_atom_optimality():
    """Any K pairwise coprime atoms > 1 compose at most Psi(X, p_K) integers <= X."""
    rng = random.Random(11)
    P = primes_upto(200)
    for X in (10 ** 4, 10 ** 5):
        for _ in range(40):
            K = rng.randrange(1, 6)
            atoms = []
            while len(atoms) < K:
                a = rng.randrange(2, 120)
                if all(math.gcd(a, b) == 1 for b in atoms):
                    atoms.append(a)
            assert atom_product_count(atoms, X) <= psi_count(X, P[K - 1])
        assert atom_product_count(P[:5], X) == psi_count(X, P[4])


def test_lemma_KI_kernel_inclusion():
    """Extra free-side columns (large atoms) never enlarge the left kernel."""
    rng = random.Random(3)
    for _ in range(200):
        R = rng.randrange(2, 12)
        small = [rng.getrandbits(6) for _ in range(R)]
        large = [rng.getrandbits(5) for _ in range(R)]
        combined = [s | (l << 6) for s, l in zip(small, large)]
        assert left_kernel_rank(combined) <= left_kernel_rank(small)


def test_lemma_DSC_size_floor():
    """u = t v (mod N), (u, v) not proportional to (t, 1) => |u| + |t||v| >= N."""
    rng = random.Random(9)
    N = balanced_semiprime(32, 0).N
    for _ in range(300):
        t = rng.randrange(2, N // 3)
        v = rng.randrange(-1000, 1000)
        k = rng.randrange(-3, 4)
        u = t * v + k * N
        if u != t * v:
            assert abs(u) + abs(t) * abs(v) >= N


def test_inv_square_sum_fixed():
    s = inv_square_sum_fixed([2, 3, 5], 60)
    exact = (1 << 60) * (1 / 4 + 1 / 9 + 1 / 25)
    assert abs(s - exact) <= 3


def test_trial_factor():
    P = primes_upto(50)
    assert trial_factor(2 ** 5 * 3 * 47, P) == {2: 5, 3: 1, 47: 1}
    assert trial_factor(53, P) is None
    assert trial_factor(-12, P) == {2: 2, 3: 1}
