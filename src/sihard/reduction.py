"""Prop XS: X3C -> "S-smooth integer in an interval" with an explicit prime list S (parsimonious).

Reduction (results note, Prop XS). Universe {0, ..., 3q-1}; sets T_1..T_m (3-subsets, deduplicated).
  b = q + 1; w_j = sum_{e in T_j} b^e; W = sum_{e < 3q} b^e; Bmax = b^{3q}.
  R (chosen), R' = R // (4q); Amax = (Bmax + 1) * R; Q >= 2^{q+2} * Amax^2 // R + 2.
  p_j: a prime in [Q + w_j R, Q + w_j R + R').
  x = Q^q + Q^{q-1} R W,  h = Q^{q-1} R // 2.
Claim: S-smooth z in [x, x+h]  <->  exact covers, bijectively (z = product of the cover's primes).
Prime finding is randomized (seeded). The theorem uses AKS (zero-error); the code uses sympy.isprime (BPSW above 2^64),
so the code itself is not a certified zero-error reduction. See build_instance_theory for the provable parameters.
"""
from __future__ import annotations

import itertools
import random
from typing import List, Sequence, Tuple

import sympy


def build_instance(q: int, sets: Sequence[Tuple[int, int, int]], rng: random.Random, R: int | None = 1 << 20):
    """R = None: the provable choice R = K^2, K = 2^{q+2}(b^{3q}+1)^2, so Q ~ K^3 and R' ~ Q^{2/3}/(4q)
    (primes of density ~1/log Q in [Q', Q' + R') by Huxley). R = 2^20 is the registered E4 toy choice (heuristic
    prime supply)."""
    sets = sorted(set(tuple(sorted(s)) for s in sets))
    b = q + 1
    w = [sum(b ** e for e in s) for s in sets]
    W = sum(b ** e for e in range(3 * q))
    Bmax = b ** (3 * q)
    if R is None:
        K = (1 << (q + 2)) * (Bmax + 1) ** 2
        R = K * K
    Rp = R // (4 * q)
    Amax = (Bmax + 1) * R
    Q = (1 << (q + 2)) * Amax * Amax // R + 2
    primes = []
    for wj in w:
        lo = Q + wj * R
        for _ in range(100000):
            c = lo + rng.randrange(Rp)
            if sympy.isprime(c):
                primes.append(c)
                break
        else:  # pragma: no cover
            raise RuntimeError("no prime found")
    x = Q ** q + Q ** (q - 1) * R * W
    h = Q ** (q - 1) * R // 2
    return {"sets": sets, "primes": primes, "x": x, "h": h, "Q": Q, "R": R}


def build_instance_theory(q: int, sets: Sequence[Tuple[int, int, int]], rng: random.Random):
    """The Las Vegas reduction of Prop XS as analysed: R = K^2 (Huxley prime supply), unbounded retries.

    Primality here is sympy.isprime (BPSW above 2^64; no known counterexample, but not a proof). The theorem uses a
    deterministic primality test (AKS) in this step, which makes the reduction zero-error.
    """
    sets = sorted(set(tuple(sorted(s)) for s in sets))
    b = q + 1
    w = [sum(b ** e for e in s) for s in sets]
    W = sum(b ** e for e in range(3 * q))
    Bmax = b ** (3 * q)
    K = (1 << (q + 2)) * (Bmax + 1) ** 2
    R = K * K
    Rp = R // (4 * q)
    Amax = (Bmax + 1) * R
    Q = (1 << (q + 2)) * Amax * Amax // R + 2
    primes = []
    for wj in w:
        lo = Q + wj * R
        while True:
            c = lo + rng.randrange(Rp)
            if sympy.isprime(c):
                primes.append(c)
                break
    return {"sets": sets, "primes": primes, "x": Q ** q + Q ** (q - 1) * R * W, "h": Q ** (q - 1) * R // 2,
            "Q": Q, "R": R}


def count_S_smooth_in_window(primes: Sequence[int], x: int, h: int) -> int:
    """Exhaustive: number of integers z in [x, x+h] all of whose prime factors lie in `primes`.

    Independent of the reduction's logic: enumerates every exponent vector with product <= x + h.
    """
    ps = sorted(primes)
    top = x + h
    count = 0

    def rec(i: int, prod: int):
        nonlocal count
        if i == len(ps):
            if prod >= x:
                count += 1
            return
        p = ps[i]
        v = prod
        while v <= top:
            rec(i + 1, v)
            v *= p

    rec(0, 1)
    return count


def count_exact_covers(q: int, sets: Sequence[Tuple[int, int, int]]) -> int:
    sets = sorted(set(tuple(sorted(s)) for s in sets))
    U = set(range(3 * q))
    c = 0
    for comb in itertools.combinations(sets, q):
        cov = [e for s in comb for e in s]
        if len(cov) == len(set(cov)) and set(cov) == U:
            c += 1
    return c


def random_x3c(q: int, m: int, rng: random.Random, plant: bool) -> List[Tuple[int, int, int]]:
    U = list(range(3 * q))
    sets = []
    if plant:
        perm = U[:]
        rng.shuffle(perm)
        sets = [tuple(sorted(perm[3 * i: 3 * i + 3])) for i in range(q)]
    while len(set(sets)) < m:
        sets.append(tuple(sorted(rng.sample(U, 3))))
    return sorted(set(sets))
