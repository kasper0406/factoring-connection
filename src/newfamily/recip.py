"""Reciprocity-computable invariants of (a, N).  OPERATIONAL: every function here sees only integers it is
given (N and public arguments); none imports the instance generator, a harness, or a factoring routine.

Exact integer arithmetic only.

* `dedekind6(a, c)`  = 6c * s(a, c), an integer (Rademacher-Grosswald), via the reciprocity law
      s(a,c) + s(c,a) = -1/4 + (a^2 + c^2 + 1) / (12 a c)      (a, c > 0 coprime),
  i.e. 2a*D(a,c) + 2c*D(c,a) = a^2 + c^2 + 1 - 3ac with D(a,c) = 6c s(a,c), and s(a,c) = s(a mod c, c).
  O(log c) Euclid steps.
* `cf(a, c)`  partial quotients of a/c (0 < a < c): a/c = [0; a_1, ..., a_r].
* `cf_stats`  length, #(a_i = 1), max bitlen, sum bitlen, alternating sum (the Rademacher/Barkan-Hickerson-Knuth
  part of 12 s(a,c)), and Euclid-remainder residue counts (class R2 of notes/newfamily_results.md).
* `jacobi(a, n)`  Jacobi symbol by reciprocity (gmpy2).
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List

import gmpy2


def dedekind6(a: int, c: int) -> int:
    """6c * s(a, c) for c >= 1, gcd(a, c) = 1.  Exact integer; iterative reciprocity."""
    if c <= 0:
        raise ValueError("c >= 1")
    a %= c
    if gcd(a, c) != 1:
        raise ValueError("gcd(a, c) must be 1")
    # Unroll: D(a,c) = (a^2 + c^2 + 1 - 3ac - 2c D(c mod a, a)) / (2a), D(0 or anything, 1) = 0.
    stack = []
    while c > 1:
        stack.append((a, c))
        a, c = c % a, a
    D = 0  # D(., 1) = 6 * 1 * s(., 1) = 0
    while stack:
        a, c = stack.pop()
        num = a * a + c * c + 1 - 3 * a * c - 2 * c * D
        assert num % (2 * a) == 0
        D = num // (2 * a)
    return D


def cf(a: int, c: int) -> List[int]:
    """Partial quotients a_1..a_r of a/c for 0 < a < c, gcd(a,c)=1 (last a_r >= 2 unless a/c = 1/1)."""
    out = []
    x, y = c, a % c
    while y:
        out.append(x // y)
        x, y = y, x % y
    return out


def euclid_remainders(c: int, a: int) -> List[int]:
    """r_0 = c, r_1 = a mod c, r_{i+1} = r_{i-1} mod r_i, stopping before 0."""
    rs = [c, a % c]
    while rs[-1]:
        rs.append(rs[-2] % rs[-1])
    return rs[:-1]


def cf_stats(a: int, c: int) -> Dict[str, int]:
    q = cf(a, c)
    rs = euclid_remainders(c, a)
    return {
        "len": len(q),
        "n_a1": sum(1 for x in q if x == 1),
        "maxbits": max(x.bit_length() for x in q),
        "sumbits": sum(x.bit_length() for x in q),
        "altsum": sum(x if i % 2 == 0 else -x for i, x in enumerate(q)),
        # Euclid-automaton statistics (remainder sequence mod 4 and mod 3), class R2
        "rem1mod4": sum(1 for r in rs if r % 4 == 1),
        "rem0mod3": sum(1 for r in rs if r % 3 == 0),
    }


def jacobi(a: int, n: int) -> int:
    return int(gmpy2.jacobi(a % n, n))
