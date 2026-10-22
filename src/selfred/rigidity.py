"""OPERATIONAL. Theorem SR2 (rigidity): an exact algebraic self-reduction is already a factoring algorithm.

Input: N, an integer N2 and integer polynomials A, B (coefficient lists, A[i] = coeff of x^i) with the promise
N2 = A(p) * B(q) for the two prime factors p, q of N (in some order).  We form

    P(x) = A(x) * x^b B(N/x) - N2 * x^b,   b = deg B,      P in Z[x],   P(p) = 0,

and return the nontrivial divisors of N among the integer roots of P (and of the swapped polynomial).
P is identically zero iff A = a x^j and B = b x^j (the multiplier class); then nothing can be learned.

Integer roots are found exactly by recursive sign-change bracketing (monotone pieces between the brackets of
the derivative's sign changes) applied to the squarefree part.  No integer factoring routine is called.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import List, Optional, Sequence, Tuple

Poly = List[int]


# ----------------------------------------------------------------------------- integer polynomial helpers
def trim(P: Sequence[int]) -> Poly:
    P = list(P)
    while P and P[-1] == 0:
        P.pop()
    return P


def peval(P: Sequence[int], x: int) -> int:
    r = 0
    for c in reversed(P):
        r = r * x + c
    return r


def pmul(A: Sequence[int], B: Sequence[int]) -> Poly:
    if not A or not B:
        return []
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                out[i + j] += a * b
    return trim(out)


def psub(A: Sequence[int], B: Sequence[int]) -> Poly:
    n = max(len(A), len(B))
    return trim([(A[i] if i < len(A) else 0) - (B[i] if i < len(B) else 0) for i in range(n)])


def pderiv(P: Sequence[int]) -> Poly:
    return trim([i * P[i] for i in range(1, len(P))])


def content(P: Sequence[int]) -> int:
    g = 0
    for c in P:
        g = gcd(g, c)
    return g


def primitive(P: Sequence[int]) -> Poly:
    P = trim(P)
    if not P:
        return []
    g = content(P)
    if P[-1] < 0:
        g = -g
    return [c // g for c in P]


def _qdivmod(A: List[Fraction], B: List[Fraction]) -> Tuple[List[Fraction], List[Fraction]]:
    A = list(A)
    q = [Fraction(0)] * max(len(A) - len(B) + 1, 1)
    while len(A) >= len(B) and any(A):
        c = A[-1] / B[-1]
        k = len(A) - len(B)
        q[k] = c
        for i, b in enumerate(B):
            A[i + k] -= c * b
        while A and A[-1] == 0:
            A.pop()
    return q, A


def _to_int_poly(F: List[Fraction]) -> Poly:
    den = 1
    for c in F:
        den = den * c.denominator // gcd(den, c.denominator)
    return primitive([int(c * den) for c in F])


def pgcd(A: Sequence[int], B: Sequence[int]) -> Poly:
    """Primitive gcd over Q (Euclid with Fractions; fine for the small degrees used here)."""
    a = [Fraction(c) for c in trim(A)]
    b = [Fraction(c) for c in trim(B)]
    while b and any(b):
        _, r = _qdivmod(a, b)
        a, b = b, r
    return _to_int_poly(a) if a else []


def squarefree_part(P: Sequence[int]) -> Poly:
    P = primitive(P)
    if len(P) <= 2:
        return P
    g = pgcd(P, pderiv(P))
    if len(g) <= 1:
        return P
    q, r = _qdivmod([Fraction(c) for c in P], [Fraction(c) for c in g])
    assert not any(r)
    return _to_int_poly(q)


def _sgn(v: int) -> int:
    return (v > 0) - (v < 0)


def sign_change_brackets(P: Sequence[int], lo: int, hi: int) -> List[int]:
    """Integers c in [lo, hi-1] with sgn P(c) * sgn P(c+1) <= 0 and P(c) != 0 or c a root.

    Every real root of odd multiplicity of P in [lo, hi] lies in some returned [c, c+1].  Exact.
    """
    P = trim(P)
    if len(P) <= 1 or lo >= hi:
        return []
    if len(P) == 2:  # a0 + a1 x
        a0, a1 = P
        # root at -a0/a1
        num, den = -a0, a1
        if den < 0:
            num, den = -num, -den
        c = num // den  # floor
        out = []
        for cc in (c - 1, c, c + 1):
            if lo <= cc < hi and _sgn(peval(P, cc)) * _sgn(peval(P, cc + 1)) <= 0:
                out.append(cc)
        return sorted(set(out))
    crit = sign_change_brackets(pderiv(P), lo, hi)
    # split points: P is monotone on each [a, b] between consecutive points
    pts = [lo]
    for c in crit:
        pts += [c, c + 1]
    pts.append(hi)
    pts = sorted(set(x for x in pts if lo <= x <= hi))
    out = set()
    for a, b in zip(pts, pts[1:]):
        sa, sb = _sgn(peval(P, a)), _sgn(peval(P, b))
        if sa == 0:
            out.add(a if a < hi else a - 1)
        if sb == 0:
            out.add(b - 1 if b - 1 >= lo else b)
        if sa * sb < 0:
            # monotone on [a, b]: binary search for c with a sign change on [c, c+1]
            l, r = a, b  # invariant: sgn P(l) = sa, sgn P(r) = sb
            while r - l > 1:
                m = (l + r) // 2
                sm = _sgn(peval(P, m))
                if sm == 0:
                    l = r = m
                    break
                if sm == sa:
                    l = m
                else:
                    r = m
            out.add(min(l, hi - 1))
    return sorted(c for c in out if lo <= c < hi)


def integer_roots(P: Sequence[int], lo: int, hi: int) -> List[int]:
    """All integer roots of P in [lo, hi] (P not identically zero). Exact."""
    S = squarefree_part(P)
    if not S:
        raise ValueError("zero polynomial")
    roots = set()
    for c in sign_change_brackets(S, lo, hi):
        for x in (c, c + 1):
            if lo <= x <= hi and peval(S, x) == 0:
                roots.add(x)
    for x in (lo, hi):
        if peval(S, x) == 0:
            roots.add(x)
    return sorted(roots)


# ----------------------------------------------------------------------------- Theorem SR2
def rigidity_polynomial(N: int, N2: int, A: Sequence[int], B: Sequence[int]) -> Poly:
    """P(x) = A(x) * x^b B(N/x) - N2 x^b, with x^b B(N/x) = sum_j B_j N^j x^(b-j)."""
    B = trim(B)
    b = len(B) - 1
    Bt = [0] * (b + 1)
    for j, Bj in enumerate(B):
        Bt[b - j] = Bj * N ** j
    xb = [0] * b + [N2]
    return psub(pmul(A, Bt), xb)


def is_degenerate(N: int, A: Sequence[int], B: Sequence[int]) -> bool:
    """True iff A(x) B(N/x) is constant in x (then P == 0 for the matching N2): the multiplier class."""
    P = rigidity_polynomial(N, 0, A, B)  # A(x) x^b B(N/x)
    nz = [i for i, c in enumerate(P) if c]
    return len(nz) <= 1 and (not nz or nz[0] == len(trim(B)) - 1)


def solve_rigid(N: int, N2: int, A: Sequence[int], B: Sequence[int]) -> Optional[int]:
    """Return a nontrivial factor of N, or None (degenerate map or broken promise)."""
    for (AA, BB) in ((A, B), (B, A)):  # the map may attach A to either factor
        P = rigidity_polynomial(N, N2, AA, BB)
        if not P:
            continue
        for r in integer_roots(P, 2, N - 1):
            if N % r == 0 and 1 < r < N:
                return r
    return None
