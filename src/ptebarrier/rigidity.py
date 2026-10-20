"""Lemma R (rigidity) and companions.  Exact integer arithmetic only.

Polynomials are lists of Python ints, lowest degree first, with no trailing zeros
(the zero polynomial is []).

Lemma R.  P, P' monic of degree M, split over Z, roots in [-B, B]; Q, Q' nonzero in Z[X] of
degree <= j; G = Q P - Q' P' with deg G <= K-1 and ||G||_inf <= h.  Put
Delta = 2^(K-1) sqrt(K) h B^(K-1), L = floor(log2 Delta).
  j = 0 (Q = Q' = 1):  M >= K + 2L + 3        ==>  P = P'.
  general j:           M >= K + 2L + 4j + 3   ==>  Q P = Q' P'.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, isqrt
from typing import Iterable, Sequence

Poly = list


# ---------------------------------------------------------------- basic poly arithmetic
def trim(a: Sequence[int]) -> Poly:
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def deg(a: Sequence[int]) -> int:
    a = trim(a)
    return len(a) - 1 if a else -1


def add(a, b) -> Poly:
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def neg(a) -> Poly:
    return [-c for c in a]


def sub(a, b) -> Poly:
    return add(a, neg(b))


def mul(a, b) -> Poly:
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for k, y in enumerate(b):
                out[i + k] += x * y
    return trim(out)


def evalp(a, x):
    v = 0
    for c in reversed(a):
        v = v * x + c
    return v


def height(a) -> int:
    return max((abs(c) for c in a), default=0)


def from_roots(roots: Iterable[int]) -> Poly:
    p = [1]
    for r in roots:
        p = mul(p, [-r, 1])
    return p


def divmod_monic(a, f):
    """Exact division with remainder by a monic f over Z."""
    a = trim(a)
    f = trim(f)
    assert f and f[-1] == 1
    df = len(f) - 1
    q = [0] * max(len(a) - df, 0)
    r = list(a)
    for i in range(len(a) - 1, df - 1, -1):
        c = r[i]
        if c:
            q[i - df] = c
            for k in range(df + 1):
                r[i - df + k] -= c * f[k]
    return trim(q), trim(r[:df] if df > 0 else [])


# ---------------------------------------------------------------- thresholds
def floor_log2_delta(B: int, K: int, h: int) -> int:
    """floor(log2(2^(K-1) * sqrt(K) * h * B^(K-1))), exactly (via the square)."""
    assert B >= 1 and K >= 1 and h >= 1
    sq = 4 ** (K - 1) * K * h * h * B ** (2 * (K - 1))
    e = sq.bit_length() - 1  # floor(log2 sq)
    return e // 2


def regime_ok(M: int, B: int, K: int, h: int, j: int = 0) -> bool:
    L = floor_log2_delta(B, K, h)
    return M >= K + 2 * L + 4 * j + 3


def h_R(B: int, M: int, K: int, j: int = 0) -> int:
    """Largest h >= 1 for which the Lemma-R hypothesis on M holds (0 if none)."""
    if not regime_ok(M, B, K, 1, j):
        return 0
    lo, hi = 1, 2
    while regime_ok(M, B, K, hi, j):
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if regime_ok(M, B, K, mid, j):
            lo = mid
        else:
            hi = mid
    return lo


# ---------------------------------------------------------------- multiset helpers
def multiset_of(roots: Iterable[int]) -> dict:
    d: dict = {}
    for r in roots:
        d[r] = d.get(r, 0) + 1
    return d


def split_common(A: Sequence[int], A2: Sequence[int]):
    """Return (common, U, U2) with A = common + U, A2 = common + U2 and U, U2 disjoint."""
    ma, mb = multiset_of(A), multiset_of(A2)
    common, U, U2 = [], [], []
    for v in set(ma) | set(mb):
        c = min(ma.get(v, 0), mb.get(v, 0))
        common += [v] * c
        U += [v] * (ma.get(v, 0) - c)
        U2 += [v] * (mb.get(v, 0) - c)
    return sorted(common), sorted(U), sorted(U2)


def int_roots(q: Sequence[int]) -> set:
    """Integer roots of a nonzero integer polynomial (brute force over divisors of the
    lowest nonzero coefficient)."""
    q = trim(q)
    if not q:
        raise ValueError("zero polynomial")
    out = set()
    i0 = 0
    while q[i0] == 0:
        i0 += 1
    if i0 > 0:
        out.add(0)
    c = abs(q[i0])
    for dvs in range(1, isqrt(c) + 1):
        if c % dvs == 0:
            for t in (dvs, c // dvs):
                for s in (t, -t):
                    if evalp(q, s) == 0:
                        out.add(s)
    return out


# ---------------------------------------------------------------- proof certificate
def lemma_R_certificate(A, A2, Q=(1,), Q2=(1,), B=None):
    """Recompute every quantity in the proof of Lemma R for one pair and check each claim.

    Returns a dict with the proof's intermediate quantities and booleans for each step.
    The final field 'consistent' is True iff no step of the proof is contradicted (if the
    regime hypothesis holds and QP != Q'P', the proof says this cannot happen, so then
    'regime' and 'distinct' cannot both be True).
    """
    A, A2 = sorted(A), sorted(A2)
    if B is None:
        B = max(1, max(abs(a) for a in A + A2))
    P, P2 = from_roots(A), from_roots(A2)
    Q, Q2 = trim(Q), trim(Q2)
    G = sub(mul(Q, P), mul(Q2, P2))
    M = len(A)
    j = max(deg(Q), deg(Q2), 0)
    out = {"M": M, "B": B, "j": j, "G": G, "distinct": bool(G)}
    if not G:
        out.update(consistent=True, regime=None)
        return out
    K = deg(G) + 1
    h = height(G)
    L = floor_log2_delta(B, K, h)
    out.update(K=K, h=h, L=L, regime=regime_ok(M, B, K, h, j))
    common, U, U2 = split_common(A, A2)
    g = from_roots(common)
    Uq, U2q = from_roots(U), from_roots(U2)
    Gt = sub(mul(Q, Uq), mul(Q2, U2q))
    # step 1: G = g * Gt exactly
    s1 = mul(g, Gt) == G
    # step 2: |Gt(a)| <= Delta for |a| <= B  (check at the roots of U, U2; compare squares)
    Dsq = 4 ** (K - 1) * K * h * h * B ** (2 * (K - 1))
    s2 = all(evalp(Gt, a) ** 2 <= Dsq for a in set(U) | set(U2))
    # step 3: at good roots a of U (Q2(a) != 0): 0 < |U2(a)| <= |Gt(a)|
    s3 = True
    far_U2 = None
    for a in set(U):
        if evalp(Q2, a) != 0:
            v = evalp(U2q, a)
            s3 &= v != 0 and abs(v) <= abs(evalp(Gt, a))
            nf = sum(1 for b in U2 if abs(b - a) >= 2)
            s3 &= nf <= L
            far_U2 = nf if far_U2 is None else max(far_U2, nf)
    for a in set(U2):
        if evalp(Q, a) != 0:
            v = evalp(Uq, a)
            s3 &= v != 0 and abs(v) <= abs(evalp(Gt, a))
            s3 &= sum(1 for b in U if abs(b - a) >= 2) <= L
    # step 4: Vandermonde/PTE: union support (incl. complex roots of Q,Q2, counted by degree)
    #         must be >= J + 2 where J = number of agreeing top coefficients of QU and Q2U2.
    QU, Q2U2 = mul(Q, Uq), mul(Q2, U2q)
    supp_bound = len(set(U) | set(U2)) + max(deg(Q), 0) + max(deg(Q2), 0)
    if deg(QU) == deg(Q2U2) and QU[-1] == Q2U2[-1]:
        J = deg(QU) - 1 - deg(Gt)
        s4 = supp_bound >= J + 2
    else:
        # degree/leading-coefficient branch of the proof: deg Gt = max(deg QU, deg Q2U2) >= r,
        # which forces M <= K - 1 (so the regime cannot hold)
        J = None
        s4 = deg(Gt) == max(deg(QU), deg(Q2U2)) and deg(Gt) >= len(U) and M <= K - 1
    out.update(J=J, supp_bound=supp_bound, steps=(s1, s2, s3, s4))
    out["consistent"] = s1 and s2 and s3 and s4 and not out["regime"]
    return out


# ---------------------------------------------------------------- Theorem W invariant
def psi_invariant(h, P, r):
    """The rational function h*P/r, as a reduced pair (num, den) over Q, normalised so that
    den is primitive with positive leading coefficient.  Theorem W: in the regime, all members
    of a single-presentation family have the same invariant."""
    from sympy import Poly as SPoly, symbols, cancel
    X = symbols("X")
    num = SPoly(list(reversed(mul(h, P))), X).as_expr()
    den = SPoly(list(reversed(trim(r))), X).as_expr()
    return cancel(num / den)


# ---------------------------------------------------------------- Prop S (View A support)
def view_a_small_value_set(F, m: int, Y: int, rho: int, H: int):
    """S = {a in [m-Y, m+Y] : |F(a)| <= (rho+1) H Bmax^rho} where Bmax = |m| + Y.

    Prop S: if (rho+1) H Bmax^rho (d+1) Y^d < |F(m)| then |S| <= d = deg F."""
    Bmax = abs(m) + Y
    thr = (rho + 1) * H * max(1, Bmax) ** rho
    return [a for a in range(m - Y, m + Y + 1) if abs(evalp(F, a)) <= thr]


def prop_S_hypothesis(F, m: int, Y: int, rho: int, H: int) -> bool:
    d = deg(F)
    Bmax = abs(m) + Y
    Fm = abs(evalp(F, m))
    return Fm > 0 and (rho + 1) * H * max(1, Bmax) ** rho * (d + 1) * Y ** d < Fm


# ---------------------------------------------------------------- bounds used in the note
def injection_bound(q: int, k: int) -> int:
    """max fiber of (e_1..e_s), s = M-k, over multisets from q consecutive integers."""
    return comb(q + k - 1, k)


def mahler_l1_bound(K: int, h: int) -> Fraction:
    """Upper bound 2^(K-1) sqrt(K) h on ||G/g||_1 (used in step 2), returned squared."""
    return Fraction(4 ** (K - 1) * K * h * h)
