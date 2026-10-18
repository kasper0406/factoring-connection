"""Constructive covering search (MODEL CODE; N is only a size parameter).

For a support I (exponents of p; negative k means q^{|k|} = N^{|k|} p^k) and
window W, sweep the balanced range p in [ceil(sqrt(N/c)), floor(sqrt(N))]
left to right.  At the current uncovered integer x, look for the longest run
[x, x+L] and a lattice target g (integer combination of the basis functions,
optionally the integer-valued lattice of a residue class) such that every
integer p in the run (in the class, if any) satisfies
      ceil(min g) <= g(p) <= ceil(min g) + W - 1,
i.e. the window {U, ..., U+W-1}, U = ceil(min g over the run), contains g(p).
That is exactly what a one-sided collision search needs: if p is the hidden
factor, g(p) is an integer in the window and alpha^{g(p)} has a computable
handle.

Lattice search: LLL (flint) on the exact, rationally-rounded Taylor images
(c_n h^n)_{n=1..K} at the run centre -- a heuristic CHOICE of candidate; the
coverage of every accepted target is then VERIFIED EXACTLY (rational
arithmetic; critical points located with rigorous arb root balls).  The
result is therefore a valid covering (an upper bound on the optimum cost).
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, floor, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

import flint

from .calibration import ceil_sqrt_ratio

Term = Dict[int, Fraction]          # exponent of p -> coefficient


# --------------------------------------------------------------------------
# basis functions
# --------------------------------------------------------------------------
def monomial_basis(I: Sequence[int], N: int) -> List[Term]:
    out = []
    for k in I:
        out.append({k: Fraction(1)} if k > 0 else {k: Fraction(N ** (-k))})
    return out


def class_basis(I: Sequence[int], N: int, rows: List[List[int]], scale: int) -> List[Term]:
    """Basis of the integer-valued lattice (rows in monomial coordinates,
    divided by `scale`)."""
    mono = monomial_basis(I, N)
    out = []
    for row in rows:
        t: Term = {}
        for a, f in zip(row, mono):
            if a:
                for k, c in f.items():
                    t[k] = t.get(k, Fraction(0)) + Fraction(a, scale) * c
        out.append(t)
    return out


def combine(basis: List[Term], a: Sequence[int]) -> Term:
    t: Term = {}
    for ai, f in zip(a, basis):
        if ai:
            for k, c in f.items():
                t[k] = t.get(k, Fraction(0)) + ai * c
    return {k: c for k, c in t.items() if c}


def gen_binom(k: int, n: int) -> Fraction:
    num = 1
    for j in range(n):
        num *= (k - j)
    return Fraction(num, 1) / Fraction(_fact(n))


def _fact(n: int) -> int:
    r = 1
    for j in range(2, n + 1):
        r *= j
    return r


def taylor(t: Term, c: Fraction, nmax: int) -> List[Fraction]:
    """Taylor coefficients g^{(n)}(c)/n!, n = 0..nmax."""
    out = []
    for n in range(nmax + 1):
        s = Fraction(0)
        for k, co in t.items():
            s += co * gen_binom(k, n) * c ** (k - n)
        out.append(s)
    return out


def evaluate(t: Term, p: int) -> Fraction:
    s = Fraction(0)
    for k, c in t.items():
        s += c * (Fraction(p) ** k)
    return s


# --------------------------------------------------------------------------
# exact oscillation of g over integers of [x0, x1] (optionally a class)
# --------------------------------------------------------------------------
def _arb_to_interval(a) -> Tuple[Fraction, Fraction]:
    m, e = a.mid().man_exp()
    mid = Fraction(int(m)) * (Fraction(2) ** int(e))
    rr = a.rad()
    m2, e2 = rr.upper().mid().man_exp()
    rad = Fraction(int(m2)) * (Fraction(2) ** int(e2))
    return mid - rad, mid + rad


def critical_points(t: Term) -> List[Tuple[Fraction, Fraction]]:
    """Rigorous enclosures of the real parts of all complex roots of g'."""
    ks = [k for k in t if k != 0]
    if not ks:
        return []
    emin = min(min(ks) - 1, 0)
    # p^{-emin} * g'(p) = sum k c_k p^{k-1-emin}
    coeffs: Dict[int, Fraction] = {}
    for k, c in t.items():
        if k:
            coeffs[k - 1 - emin] = coeffs.get(k - 1 - emin, Fraction(0)) + k * c
    deg = max(coeffs)
    den = 1
    for c in coeffs.values():
        den = den * c.denominator // _gcd(den, c.denominator)
    ints = [int(coeffs.get(j, Fraction(0)) * den) for j in range(deg + 1)]
    poly = flint.fmpz_poly(ints)
    if poly.degree() <= 0:
        return []
    old = flint.ctx.prec
    flint.ctx.prec = 256
    try:
        roots = poly.complex_roots()
    finally:
        flint.ctx.prec = old
    return [_arb_to_interval(r.real) for r, _mult in roots]


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def class_points(x0: int, x1: int, M: int, r: int) -> Tuple[int, int]:
    """first and last integers in [x0,x1] congruent to r mod M (or None)."""
    a = x0 + ((r - x0) % M)
    if a > x1:
        return None
    b = x1 - ((x1 - r) % M)
    return a, b


def oscillation_ok(t: Term, x0: int, x1: int, W: int, crit=None,
                   M: int = 1, r: int = 0) -> Tuple[bool, Optional[int]]:
    """Exact test: all class integers p in [x0,x1] have g(p) in
    {U..U+W-1}, U = ceil(min g).  Returns (ok, U)."""
    cp = class_points(x0, x1, M, r)
    if cp is None:
        return True, None
    a, b = cp
    cands = {a, b}
    if crit is None:
        crit = critical_points(t)
    for lo, hi in crit:
        if hi < a - M or lo > b + M:
            continue
        flo, chi = floor(lo), -floor(-hi)
        for z in (flo, chi):
            # nearest class points around z
            base = z - ((z - r) % M)
            for y in (base - M, base, base + M, base + 2 * M):
                if a <= y <= b:
                    cands.add(y)
        if chi - flo > 4 * M:          # very wide ball: be safe
            return False, None
    vals = [evaluate(t, y) for y in cands]
    mn, mx = min(vals), max(vals)
    U = -((-mn.numerator) // mn.denominator)       # ceil
    top = mx.numerator // mx.denominator           # floor
    return (top - U + 1 <= W), U


# --------------------------------------------------------------------------
# LLL candidate search
# --------------------------------------------------------------------------
def _round(fr: Fraction) -> int:
    return (2 * fr.numerator + fr.denominator) // (2 * fr.denominator)


def lll_candidates(basis: List[Term], x0: int, x1: int, W: int, extra: int = 1) -> List[List[int]]:
    K = len(basis)
    c = Fraction(x0 + x1, 2)
    h = Fraction(max(x1 - x0, 1), 2)
    S = Fraction(1 << 48, max(W, 1))
    rows = []
    for j, f in enumerate(basis):
        tc = taylor(f, c, K + extra)
        img = [_round(S * tc[n] * h ** n) for n in range(1, K + 1 + extra)]
        rows.append([1 if i == j else 0 for i in range(K)] + img)
    red = flint.fmpz_mat(rows).lll().tolist()
    vecs = [[int(v) for v in row[:K]] for row in red]
    cands = []
    for v in vecs:
        if any(v):
            cands.append(v)
    if len(vecs) >= 2:
        for s1 in (1, -1):
            w = [a + s1 * b for a, b in zip(vecs[0], vecs[1])]
            if any(w):
                cands.append(w)
    return cands


def best_target(basis, x0, x1, W, M=1, r=0):
    for v in lll_candidates(basis, x0, x1, W):
        t = combine(basis, v)
        if not t:
            continue
        ok, U = oscillation_ok(t, x0, x1, W, M=M, r=r)
        if ok:
            return v, U
    return None, None


# --------------------------------------------------------------------------
# greedy sweep
# --------------------------------------------------------------------------
def sweep(basis: List[Term], plo: int, phi: int, W: int, M: int = 1, r: int = 0,
          max_targets: int = None, record: bool = False) -> Dict:
    """Greedy left-to-right covering of class points of [plo, phi]."""
    x = plo + ((r - plo) % M) if M > 1 else plo
    n_targets, L_prev, checks = 0, M, 0
    targets = []
    while x <= phi:
        # gallop: find largest L (multiple of M) with a successful target
        lo_ok = 0                        # L = 0 always works (single point)
        v_ok, U_ok = None, None
        L = max(M, L_prev)
        hi_bad = None
        while True:
            x1 = min(phi, x + L)
            v, U = best_target(basis, x, x1, W, M, r)
            checks += 1
            if v is not None:
                lo_ok, v_ok, U_ok = x1 - x, v, U
                if x1 == phi:
                    break
                L *= 2
            else:
                hi_bad = L
                break
        if hi_bad is not None:
            lo, hi = lo_ok, hi_bad
            while hi - lo > max(M, (lo // 64) // M * M):
                mid = (lo + hi) // 2
                mid -= mid % M
                if mid <= lo:
                    break
                v, U = best_target(basis, x, x + mid, W, M, r)
                checks += 1
                if v is not None:
                    lo, v_ok, U_ok = mid, v, U
                else:
                    hi = mid
            lo_ok = lo
        if v_ok is None:
            # single point always coverable by any nonzero vector with W >= 1
            v_ok = [1] + [0] * (len(basis) - 1)
        n_targets += 1
        if record:
            targets.append((x, x + lo_ok, v_ok, U_ok))
        L_prev = max(M, lo_ok)
        x = x + lo_ok + M
        if max_targets and n_targets >= max_targets:
            return {"n_targets": n_targets, "covered_to": x, "checks": checks,
                    "complete": False, "targets": targets}
    return {"n_targets": n_targets, "covered_to": x, "checks": checks,
            "complete": True, "targets": targets}


def balanced_range(N: int, c_num: int = 2, c_den: int = 1) -> Tuple[int, int]:
    return ceil_sqrt_ratio(N * c_den, c_num), isqrt(N)


def harvey_cost(n_targets: int, W: int) -> Dict:
    """min over baby-list size m of m + n_targets * ceil(W/m)."""
    best = None
    cands = {1, W, max(1, isqrt(n_targets * W))}
    m = 1
    while m <= W:
        cands.add(m)
        m *= 2
    for m in cands:
        if m < 1 or m > W:
            continue
        c = m + n_targets * (-(-W // m))
        if best is None or c < best[0]:
            best = (c, m)
    return {"cost": best[0], "m": best[1]}
