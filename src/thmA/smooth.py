"""Batch smoothness (Bernstein smooth parts; notes/lv6_results.md §6.2 Lemma BS).  OPERATIONAL.

smooth_parts(xs, P, e): product tree of the xs, P mod (product) pushed down a remainder tree, then
s_i = gcd(x_i, (P mod x_i)^(2^e) mod x_i), with 2^e >= log2(max x).  x_i is B-smooth iff s_i = x_i.
Every B-smooth value is recognised; the test is deterministic.
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List, Sequence


def primorial(primes: Sequence[int]) -> int:
    level = list(primes) or [1]
    while len(level) > 1:
        level = [level[i] * level[i + 1] if i + 1 < len(level) else level[i] for i in range(0, len(level), 2)]
    return level[0]


def smooth_parts(xs: Sequence[int], P: int) -> List[int]:
    xs = [abs(x) for x in xs]
    if not xs:
        return []
    assert all(x > 0 for x in xs)
    # product tree
    tree = [xs]
    while len(tree[-1]) > 1:
        cur = tree[-1]
        tree.append([cur[i] * cur[i + 1] if i + 1 < len(cur) else cur[i] for i in range(0, len(cur), 2)])
    rems = [P % tree[-1][0]]
    for lvl in range(len(tree) - 2, -1, -1):
        cur = tree[lvl]
        rems = [rems[i // 2] % cur[i] for i in range(len(cur))]
    e = max(1, max(x.bit_length() for x in xs).bit_length())
    out = []
    for x, r in zip(xs, rems):
        y = r
        for _ in range(e):
            y = y * y % x
        out.append(gcd(x, y) if x > 1 else 1)
    return out


def factor_smooth(x: int, primes: Sequence[int]) -> Dict[int, int]:
    """Factorisation of a KNOWN B-smooth |x| by trial division over the primes <= B."""
    x = abs(x)
    fac: Dict[int, int] = {}
    for p in primes:
        if x == 1:
            break
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            fac[p] = e
    assert x == 1, "factor_smooth called on a non-smooth value"
    return fac
