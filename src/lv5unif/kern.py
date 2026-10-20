"""Exact fibre counts for the pair map Phi(c) = (<c, v(P)>, <c, v(P')>) on C_H = [-H, H]^d (Lemma Kern-U).

Convention (pcest / preprint): R_c(x, y) = sum_i c_i x^(d-1-i) y^i, v(P) = (a^(d-1-i) b^i)_{0<=i<d}.

Lemma Kern-U (notes/lv5unif_results.md §3): ker Phi has the basis e_m = coefficient vector of
Q * x^(d-3-m) y^m, Q = (b x - a y)(b' x - a' y), m = 0..d-3.  e_m has bb' at index m, -(ab' + a'b) at m+1,
aa' at m+2 (lower triangular on the first d-2 coordinates).  Hence every fibre has at most
(floor(2H/(bb')) + 1)^(d-2) points.

Two independent counting methods:
- `fibre_histogram_brute`: one pass over all of C_H (independent of the lemma); for small (2H+1)^d.
- `fibre_count_lattice`: exact count of lambda in Z^(d-2) with c0 + sum lambda_m e_m in C_H, checking ALL d
  box constraints (uses Kern(ii), the kernel basis, which pcest/preprint prove; cross-checked against brute force).
No floating point.
"""
from __future__ import annotations

from itertools import product
from typing import Dict, List, Sequence, Tuple


def vvec(d: int, a: int, b: int) -> List[int]:
    return [a ** (d - 1 - i) * b ** i for i in range(d)]


def phi(c: Sequence[int], v: Sequence[int], vp: Sequence[int]) -> Tuple[int, int]:
    return (sum(ci * vi for ci, vi in zip(c, v)), sum(ci * vi for ci, vi in zip(c, vp)))


def kernel_basis(d: int, P: Tuple[int, int], Pp: Tuple[int, int]) -> List[List[int]]:
    a, b = P
    ap, bp = Pp
    top, mid, low = b * bp, -(a * bp + ap * b), a * ap
    basis = []
    for m in range(d - 2):
        e = [0] * d
        e[m] += top
        e[m + 1] += mid
        e[m + 2] += low
        basis.append(e)
    return basis


def triangular_bound(d: int, b: int, bp: int, H: int) -> int:
    return (2 * H // (b * bp) + 1) ** (d - 2)


def fibre_histogram_brute(d: int, P, Pp, H: int) -> Dict[Tuple[int, int], int]:
    v, vp = vvec(d, *P), vvec(d, *Pp)
    hist: Dict[Tuple[int, int], int] = {}
    rng = range(-H, H + 1)
    for c in product(rng, repeat=d):
        z = phi(c, v, vp)
        hist[z] = hist.get(z, 0) + 1
    return hist


def _ceil_div(p: int, q: int) -> int:
    return -((-p) // q)


def fibre_count_lattice(d: int, P, Pp, H: int, c0: Sequence[int]) -> int:
    """#{c in C_H : Phi(c) = Phi(c0)} by exact enumeration of the kernel coordinates lambda.

    Level i (0 <= i <= d-3): c_i = c0_i + top*lam_i + mid*lam_{i-1} + low*lam_{i-2} must lie in [-H, H];
    since top = bb' >= 1 this gives an exact integer interval for lam_i.  At the end the remaining two
    coordinates c_{d-2} = c0 + mid*lam_{d-3} + low*lam_{d-4}, c_{d-1} = c0 + low*lam_{d-3} are checked.
    """
    a, b = P
    ap, bp = Pp
    top, mid, low = b * bp, -(a * bp + ap * b), a * ap
    assert top >= 1
    n = d - 2
    lam = [0] * n
    count = 0

    def coord(i: int) -> int:
        s = c0[i]
        if i < n:
            s += top * lam[i]
        if 0 <= i - 1 < n:
            s += mid * lam[i - 1]
        if 0 <= i - 2 < n:
            s += low * lam[i - 2]
        return s

    def rec(i: int) -> None:
        nonlocal count
        if i == n:
            if all(-H <= coord(j) <= H for j in range(d)):
                count += 1
            return
        base = c0[i]
        if i - 1 >= 0:
            base += mid * lam[i - 1]
        if i - 2 >= 0:
            base += low * lam[i - 2]
        lo = _ceil_div(-H - base, top)
        hi = (H - base) // top
        for t in range(lo, hi + 1):
            lam[i] = t
            rec(i + 1)
        lam[i] = 0

    rec(0)
    return count
