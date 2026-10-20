"""Exact size-model floors for in-field multi-level descent with a component budget (Prop R3(C)).  DIAGNOSTIC.

Exponent bookkeeping in powers of D = log N / log log N (log log factors dropped, recurse §3.3-3.4).  All
log-sizes are D^x; we track x.  t = log N ~ D^1, Lambda = log B ~ D^alpha, U ~ D^alpha.  A value of log-size D^x is
D^y-smooth at cost L[alpha] iff x <= alpha + y (u = D^{x-y} <= U).  A value of log-size D^x has about D^{x-y}
children of log-size D^y (typical-factorisation model).

Chain sub-model (one field, degree d = D^delta, log H_f = log H_g = log N/(d+1) ~ D^{1-delta}):
  root   : rational target zeta_0 = 1 - alpha (forced by the Dixon split);
  level j: algebraic target zeta_j, value v_j = max(1 - delta, delta + zeta_j), children zeta_{j+1}:
           v_j <= alpha + zeta_{j+1},  zeta_{j+1} <= zeta_j;
  terminal (j = r): v_r <= 2 alpha (cofactor B-smooth);
  rational sides B-smooth: 1 - delta <= 2 alpha and zeta_j <= 2 alpha for j >= 1;
  component size: sum_{j<r} max(0, v_j - zeta_{j+1})  (points in the algebraic component, as a D-exponent);
  budget C: component size <= C * delta   (i.e. at most d^C points).
Feasibility is an exact LP (Fractions, phase-1 simplex, Bland's rule).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import List, Optional, Sequence, Tuple


# ------------------------------------------------------------------ exact phase-1 simplex
def lp_feasible(A: Sequence[Sequence[Fr]], b: Sequence[Fr]) -> bool:
    """Is {x >= 0 : A x <= b} non-empty?  Exact (Fractions), Bland's rule."""
    m, n = len(A), len(A[0])
    # rows: A x + s = b ; flip rows with b < 0 and add an artificial variable there
    rows, rhs, art_rows = [], [], []
    for i in range(m):
        row = [Fr(a) for a in A[i]] + [Fr(1) if k == i else Fr(0) for k in range(m)]
        bi = Fr(b[i])
        if bi < 0:
            row = [-a for a in row]
            bi = -bi
            art_rows.append(i)
        rows.append(row)
        rhs.append(bi)
    na = len(art_rows)
    ncol = n + m + na
    T = []
    basis = []
    for i in range(m):
        art = [Fr(0)] * na
        if i in art_rows:
            art[art_rows.index(i)] = Fr(1)
            basis.append(n + m + art_rows.index(i))
        else:
            basis.append(n + i)
        T.append(rows[i] + art + [rhs[i]])
    if na == 0:
        return True
    # objective: minimise sum of artificials  <=> maximise -sum; reduced costs row
    obj = [Fr(0)] * (ncol + 1)
    for k in range(na):
        obj[n + m + k] = Fr(1)
    for i in range(m):
        if basis[i] >= n + m:
            obj = [o - t for o, t in zip(obj, T[i])]
    while True:
        enter = next((j for j in range(ncol) if obj[j] < 0), None)
        if enter is None:
            break
        best, leave = None, None
        for i in range(m):
            if T[i][enter] > 0:
                ratio = T[i][-1] / T[i][enter]
                if best is None or ratio < best or (ratio == best and basis[i] < basis[leave]):
                    best, leave = ratio, i
        if leave is None:  # unbounded (cannot happen for phase 1)
            return True
        piv = T[leave][enter]
        T[leave] = [x / piv for x in T[leave]]
        for i in range(m):
            if i != leave and T[i][enter] != 0:
                f = T[i][enter]
                T[i] = [x - f * y for x, y in zip(T[i], T[leave])]
        if obj[enter] != 0:
            f = obj[enter]
            obj = [x - f * y for x, y in zip(obj, T[leave])]
        basis[leave] = enter
    return obj[-1] == 0  # minimum of sum of artificials is -obj[-1]; feasible iff 0


# ------------------------------------------------------------------ the chain sub-model
def chain_feasible(alpha: Fr, r: int, C: Optional[Fr]) -> bool:
    """Feasibility of the chain sub-model with r algebraic levels and budget d^C (C=None: no budget)."""
    # variable order: delta, zeta_1..zeta_r, v_0..v_r, e_0..e_{r-1}
    nz, nv, ne = r, r + 1, r
    n = 1 + nz + nv + ne
    D_ = 0
    Z = lambda j: 1 + (j - 1)          # zeta_j, j >= 1
    V = lambda j: 1 + nz + j            # v_j, j >= 0
    E = lambda j: 1 + nz + nv + j       # e_j
    A: List[List[Fr]] = []
    b: List[Fr] = []
    z0 = 1 - alpha

    def row(coefs, rhs):
        a = [Fr(0)] * n
        for k, c in coefs:
            a[k] += Fr(c)
        A.append(a)
        b.append(Fr(rhs))

    row([(D_, 1)], 1)                    # delta <= 1
    row([(D_, -1)], 2 * alpha - 1)       # 1 - delta <= 2 alpha
    for j in range(r + 1):
        row([(V(j), -1), (D_, -1)], -1)  # v_j >= 1 - delta
        if j == 0:
            row([(D_, 1), (V(0), -1)], -z0)          # v_0 >= delta + zeta_0
        else:
            row([(D_, 1), (Z(j), 1), (V(j), -1)], 0)  # v_j >= delta + zeta_j
            row([(Z(j), 1)], 2 * alpha)               # rational side at an algebraic node
    for j in range(r):
        row([(V(j), 1), (Z(j + 1), -1)], alpha)      # v_j <= alpha + zeta_{j+1}
        if j == 0:
            row([(Z(1), 1)], z0)                      # zeta_1 <= zeta_0
        else:
            row([(Z(j + 1), 1), (Z(j), -1)], 0)       # zeta_{j+1} <= zeta_j
        row([(V(j), 1), (Z(j + 1), -1), (E(j), -1)], 0)  # e_j >= v_j - zeta_{j+1}
    row([(V(r), 1)], 2 * alpha)                       # terminal
    if C is not None and r > 0:
        row([(E(j), 1) for j in range(r)] + [(D_, -C)], 0)  # budget
    return lp_feasible(A, b)


def floor_alpha(r: int, C: Optional[Fr], tol: Fr = Fr(1, 10 ** 10)) -> Fr:
    """Least alpha (to tol) for which the chain sub-model with r levels and budget C is feasible."""
    lo, hi = Fr(1, 3), Fr(1, 2)
    if not chain_feasible(hi, r, C):
        return Fr(1)
    if chain_feasible(lo, r, C):
        return lo
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if chain_feasible(mid, r, C):
            hi = mid
        else:
            lo = mid
    return hi


def best_floor(C: Optional[Fr], rmax: int = 12) -> Tuple[Fr, int]:
    best, arg = Fr(1), -1
    for r in range(rmax + 1):
        a = floor_alpha(r, C)
        if a < best:
            best, arg = a, r
    return best, arg


# ------------------------------------------------------------------ closed forms (to be compared)
def closed_no_budget(r: int) -> Fr:
    return Fr(r + 2, 3 * r + 5)


def closed_budget(r: int, C: Fr) -> Fr:
    """max((r+2)/(3r+5), (r+2-C)/(2r+5-2C)) when r+1 > C; else (r+2)/(3r+5).  (Derivation in the results note.)"""
    base = closed_no_budget(r)
    if r + 1 <= C:
        return base
    return max(base, Fr(r + 2 - C) / Fr(2 * r + 5 - 2 * C))
