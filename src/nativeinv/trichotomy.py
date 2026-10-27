"""OPERATIONAL. Quasi-polynomial local data and the AGG / FIBER / FACT trichotomy (Prop TRI, notes/nativeinv_results.md).

A multiplicative function f with quasi-polynomial local data has f(r) = F(r) = sum_i c_i(r mod m) r^i at primes r.
On N = pq (p != q primes, coprime to m) f(N) = F(p)F(q).  `decode` sees only N, V = f(N) and the public
description (m, c); it never sees p, q and never calls an integer-factorisation routine.  Exact arithmetic only
(integers and Fractions); real roots are isolated with Sturm sequences.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------- exact univariate polynomials (low -> high)

def _trim(P: List) -> List:
    P = list(P)
    while P and P[-1] == 0:
        P.pop()
    return P


def _ev(P: Sequence, x):
    acc = 0
    for c in reversed(P):
        acc = acc * x + c
    return acc


def _deriv(P: Sequence) -> List:
    return _trim([i * P[i] for i in range(1, len(P))])


def _divmod(A: List, B: List) -> Tuple[List, List]:
    A = [Fraction(a) for a in _trim(A)]
    B = [Fraction(b) for b in _trim(B)]
    if not B:
        raise ZeroDivisionError
    Q = [Fraction(0)] * max(len(A) - len(B) + 1, 0)
    while len(A) >= len(B) and A:
        c = A[-1] / B[-1]
        k = len(A) - len(B)
        Q[k] = c
        for i, b in enumerate(B):
            A[i + k] -= c * b
        A = _trim(A)
    return Q, A


def _gcd(A: List, B: List) -> List:
    A, B = _trim(A), _trim(B)
    while B:
        _, R = _divmod(A, B)
        A, B = B, R
    return A


def _squarefree(P: List) -> List:
    g = _gcd(P, _deriv(P))
    if len(g) <= 1:
        return [Fraction(c) for c in _trim(P)]
    Q, R = _divmod(P, g)
    assert not R
    return _trim(Q)


def _sturm(P: List) -> List[List]:
    seq = [P, _deriv(P)]
    while seq[-1]:
        _, R = _divmod(seq[-2], seq[-1])
        R = [-c for c in R]
        if not R:
            break
        seq.append(R)
    return [s for s in seq if s]


def _sign_var(seq: List[List], x) -> int:
    v, prev = 0, 0
    for s in seq:
        y = _ev(s, x)
        if y == 0:
            continue
        sg = 1 if y > 0 else -1
        if prev and sg != prev:
            v += 1
        prev = sg
    return v


def integer_roots(P: Sequence[int], lo: int, hi: int) -> List[int]:
    """All integer roots of the nonzero integer polynomial P in [lo, hi] (exact)."""
    P = _trim(list(P))
    if not P:
        raise ValueError("zero polynomial")
    if len(P) == 1:
        return []
    S = _squarefree(P)
    seq = _sturm(S)
    out = set()

    def nudge(t: Fraction) -> Fraction:
        k = 3
        while _ev(S, t) == 0:          # t is a non-integer rational root; move inside the same unit cell
            t += Fraction(1, k)
            k = k * 7 + 1
        return t

    stack = [(nudge(Fraction(2 * lo - 1, 2)), nudge(Fraction(2 * hi + 1, 2)))]
    while stack:
        a, b = stack.pop()
        cnt = _sign_var(seq, a) - _sign_var(seq, b)
        if cnt <= 0:
            continue
        if b - a <= 1:
            for x in range(int(a) - 1, int(b) + 2):
                if lo <= x <= hi and _ev(P, x) == 0:
                    out.add(x)
            continue
        if b - a > 2:
            mid = nudge(Fraction((a + b) // 2) + Fraction(1, 2))
        else:
            mid = nudge((a + b) / 2)
        stack.append((a, mid))
        stack.append((mid, b))
    return sorted(out)


# ---------------------------------------------------------------- quasi-polynomial local data

@dataclass
class QuasiPoly:
    """F(r) = sum_i coeffs[r mod m][i] * r^i for primes r coprime to m."""
    name: str
    m: int
    coeffs: Dict[int, Tuple[int, ...]]

    def units(self) -> List[int]:
        if self.m == 1:
            return [0]
        return [u for u in range(self.m) if gcd(u, self.m) == 1]

    def at(self, u: int) -> Tuple[int, ...]:
        return tuple(_trim(list(self.coeffs[u % self.m])))

    def value(self, r: int) -> int:     # harness convenience; F at an explicit prime
        return _ev(self.at(r % self.m), r)


def degrees(c: Tuple[int, ...]) -> Optional[Tuple[int, int]]:
    nz = [i for i, x in enumerate(c) if x != 0]
    return (min(nz), max(nz)) if nz else None


def pair_tier(F: QuasiPoly, u: int, v: int) -> str:
    """Tier of the residue pair (u, v) = (p mod m, q mod m): ZERO, FIBER (equal monomials) or FACT."""
    du, dv = degrees(F.at(u)), degrees(F.at(v))
    if du is None or dv is None:
        return "ZERO"
    if du[0] == du[1] == dv[0] == dv[1]:
        return "FIBER"
    return "FACT"


def global_tier(F: QuasiPoly) -> str:
    """AGG if F = c * chi(r) * r^k with chi a character mod m on all units (then f(N) = c^2 chi(N) N^k);
    else FIBER if every pair with the same product class is FIBER/ZERO, else FACT (some pair is FACT)."""
    U = F.units()
    tiers = {(u, v): pair_tier(F, u, v) for u in U for v in U}
    if any(t == "FACT" for t in tiers.values()):
        return "FACT"
    # all FIBER/ZERO: AGG iff c(u)c(v) depends only on uv (mod m) on each degree class
    vals: Dict[int, set] = {}
    for (u, v), t in tiers.items():
        w = (u * v) % F.m if F.m > 1 else 0
        if t == "ZERO":
            key = ("Z",)
        else:
            k = degrees(F.at(u))[0]
            key = (k, F.at(u)[k] * F.at(v)[k])
        vals.setdefault(w, set()).add(key)
    return "AGG" if all(len(s) == 1 for s in vals.values()) else "FIBER"


@dataclass
class DecodeResult:
    split: Optional[Tuple[int, int]]
    tiers: List[str] = field(default_factory=list)
    candidates_tested: int = 0


def decode(N: int, V: int, F: QuasiPoly) -> DecodeResult:
    """Recover {p, q} from N and V = F(p)F(q) whenever the true residue pair is FACT-tier (Prop TRI(c))."""
    res = DecodeResult(None)
    U = F.units()
    w = N % F.m if F.m > 1 else 0
    for u in U:
        for v in U:
            if F.m > 1 and (u * v) % F.m != w:
                continue
            t = pair_tier(F, u, v)
            res.tiers.append(t)
            if t != "FACT":
                continue
            cu, cv = F.at(u), F.at(v)
            Dv = len(cv) - 1
            # x = p (the factor = u mod m), q = N/x:  sum_{i,j} cu_i cv_j N^j x^(i-j+Dv) - V x^Dv = 0
            deg = (len(cu) - 1) + Dv
            G = [0] * (deg + 1)
            for i, a in enumerate(cu):
                if a == 0:
                    continue
                for j, b in enumerate(cv):
                    if b == 0:
                        continue
                    G[i - j + Dv] += a * b * N ** j
            G[Dv] -= V
            if not _trim(G):
                continue
            for x in integer_roots(G, 2, N - 1):
                res.candidates_tested += 1
                g = gcd(x, N)
                if 1 < g < N:
                    res.split = tuple(sorted((g, N // g)))
                    return res
    return res
