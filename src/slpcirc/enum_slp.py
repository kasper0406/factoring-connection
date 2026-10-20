"""Exact enumeration of constant-free division-free SLPs over Z (inputs 1 and x).

tau(F) = least number of {+,-,*} steps computing F from {1, x}. Polynomials are coefficient tuples
(low degree first, no trailing zeros). Breadth-first search over *sets* of computed values; the last
level only collects polynomials (states are not stored), so tau <= 6 needs ~120 MB.
"""
from __future__ import annotations

ONE = (1,)
X = (0, 1)


def padd(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i, c in enumerate(a):
        r[i] += c
    for i, c in enumerate(b):
        r[i] += c
    while r and r[-1] == 0:
        r.pop()
    return tuple(r)


def pneg(a):
    return tuple(-c for c in a)


def pmul(a, b):
    if not a or not b:
        return ()
    r = [0] * (len(a) + len(b) - 1)
    for i, c in enumerate(a):
        if c:
            for j, d in enumerate(b):
                r[i + j] += c * d
    return tuple(r)


def enumerate_tau(smax: int):
    """Return dict poly -> tau for all nonzero polys with tau <= smax (1 and x have tau 0)."""
    base = (ONE, X)
    tau = {ONE: 0, X: 0}
    states = {()}
    for depth in range(1, smax + 1):
        last = depth == smax
        new = set()
        for st in states:
            vals = base + st
            S = set(vals)
            n = len(vals)
            for i in range(n):
                a = vals[i]
                for j in range(i, n):
                    b = vals[j]
                    for c in (padd(a, b), pmul(a, b), padd(a, pneg(b)), padd(b, pneg(a))):
                        if c and c not in S:
                            if c not in tau:
                                tau[c] = depth
                            if not last:
                                new.add(tuple(sorted(st + (c,))))
        states = new
    return tau
