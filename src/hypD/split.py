"""Split counting (Theorem S of notes/hypD_results.md).

OPERATIONAL part (`split_in_S`): given N, a uniform square R and a membership oracle for a
set S of positive integers <= X (X^2 < N), decide whether R = u/v (mod N) with u, v in S.
All lattice vectors (u, v) with |u|, |v| <= X and u = R v (mod N) are integer multiples of
the primitive shortest vector (u0, v0) (two independent ones would have |det| < N), so it
suffices to test (j*u0, j*v0) for j <= X / max(|u0|, |v0|).  It never factors N.

HARNESS part (`theorem_s_bound`, `sample_rate`, `make_S`) builds S from a prime set G
(products of exactly k distinct primes of G) and estimates Pr_R[E_S] by sampling.
The bound needs omega(N), which is a harness-only quantity.
"""
from __future__ import annotations

import random
from math import gcd, log, isqrt
from typing import Callable, Dict, List, Sequence

from ..lemmaM.descent import gauss_short


class FactorFound(Exception):
    def __init__(self, g):
        super().__init__(g)
        self.g = int(g)


def split_in_S(R: int, N: int, X: int, inS: Callable[[int], bool]) -> bool:
    """OPERATIONAL. True iff some (u, v) with 1 <= u, v <= X, u = R v (mod N), u, v in S."""
    assert X * X < N
    u0, v0 = gauss_short(R, N)
    if v0 < 0:
        u0, v0 = -u0, -v0
    if u0 == 0 or v0 == 0:
        return False
    if u0 < 0:                      # S consists of positive integers; (u, v) both > 0 needed
        return False
    j = 1
    while j * max(u0, v0) <= X:
        u, v = j * u0, j * v0
        if (u - R * v) % N == 0 and inS(u) and inS(v):
            return True
        j += 1
    return False


def k_prime_member(G: Sequence[int], k: int) -> Callable[[int], bool]:
    """Membership in S_k(G) = {products of exactly k distinct primes of G}."""
    Gs = sorted(set(G))
    Gset = set(Gs)
    P = 1
    for g in Gs:
        P *= g

    def inS(u: int) -> bool:
        if gcd(u, P) == 1:
            return False
        cnt = 0
        for g in Gs:
            if u % g == 0:
                u //= g
                if u % g == 0:
                    return False
                cnt += 1
                if cnt > k:
                    return False
                if u == 1:
                    break
        return u == 1 and cnt == k
    return inS


def count_S(G: Sequence[int], k: int, X: int) -> int:
    """|S_k(G) ∩ [1, X]| by exact enumeration (small k)."""
    Gs = sorted(set(G))
    tot = 0

    def rec(start: int, depth: int, prod: int):
        nonlocal tot
        if depth == k:
            tot += 1
            return
        for i in range(start, len(Gs)):
            q = prod * Gs[i]
            # remaining k-depth-1 factors are each >= Gs[i+1] > Gs[i]
            if q * (Gs[i] ** (k - depth - 1)) > X:
                break
            rec(i + 1, depth + 1, q)
    rec(0, 0, 1)
    return tot


def theorem_s_bound(S_size: int, Y: int, N: int, omega: int) -> float:
    """HARNESS (needs omega). Theorem S: if S is a set of integers in [Y/2, Y] coprime to N,
    Y^2 < N, and R is uniform on the squares of (Z/N)^x, then
    Pr_R[E_S] >= |S|^4 / (2 * 2^omega * phi(N) * Y^2 * (1 + log Y)).  phi(N) <= N is used."""
    return S_size ** 4 / (2.0 * (2 ** omega) * N * Y * Y * (1.0 + log(Y)))


def count_S_window(G: Sequence[int], k: int, Y: int) -> int:
    """|S_k(G) ∩ [Y/2, Y]| (HARNESS helper)."""
    return count_S(G, k, Y) - count_S(G, k, Y // 2)


def sample_rate(N: int, X: int, inS: Callable[[int], bool], trials: int, seed: int) -> Dict:
    rng = random.Random(f"hypD-split|{N}|{seed}")
    hits = 0
    for _ in range(trials):
        r = rng.randrange(2, N - 1)
        g = gcd(r, N)
        if g != 1:
            raise FactorFound(g)          # ACCIDENTAL; probability ~ 2/p
        R = r * r % N
        if split_in_S(R, N, X, inS):
            hits += 1
    return dict(trials=trials, hits=hits)
