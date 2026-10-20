"""Frame-Hoelder (FH) tools for linked modules, degenerate or not (agent rankdegen).

DIAGNOSTIC / harness code. It takes polynomials, module generators and N; it never factors N inside any algorithm and
implements no factoring method. Reals appear only as certified arb/acb balls (python-flint); every integer reported is
exact, and every inequality reported as "holds" is certified by ball arithmetic.

Setting (as src/rank3/forms.py). f, g monic squarefree integer polynomials; K-side generators u_j = P_j(alpha),
K'-side u'_j = Q_j(beta); rows L_sigma = (P_j(alpha_sigma))_j, then (Q_j(beta_tau))_j (tagged union, even if f == g).

Theorem FH (notes/rankdegen_results.md section 3). For a probability distribution w on frames S (r+1 rows, every r of
them independent) with point marginals P_w(i in S) = (r+1) y_i / r,
    J_y := (1/r) int_{S^{r-1}} prod_i |L_i|^{-y_i}  <=  2^r C_r prod_S Phi_S^{-w_S/(r+1)},
    Phi_S = prod_{i in S} |det L_{S minus i}|.
With linking (N squarefree, linking primes unramified) and w invariant under the configuration's side-preserving
automorphisms: sum_S w_S log Phi_S >= (r-1) * min_{cross pairs} P_w(pair in S) * log N.
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction as Fr
from typing import Dict, List, Sequence, Tuple

from flint import acb, arb, ctx

from src.rank3.forms import _det, _ev, acb_roots, bean_thunder_C, poly_comb, resultant

Poly = Sequence[int]


# ----------------------------------------------------------------------------- rows and certified minors


def rows_at(f: Poly, P: Sequence[Poly], g: Poly, Q: Sequence[Poly], prec: int) -> Tuple[List[List[acb]], int]:
    old = ctx.prec
    ctx.prec = prec
    try:
        if f[-1] != 1 or g[-1] != 1:
            raise ValueError("monic f, g required")
        rows = [[_ev(p, z) for p in P] for z in acb_roots(f)] + [[_ev(q, z) for q in Q] for z in acb_roots(g)]
        return rows, len(f) - 1
    finally:
        ctx.prec = old


def galois_degree_bound(f: Poly, g: Poly) -> int:
    """[L:Q] <= d! d'! for the compositum of the splitting fields (d! if f == g)."""
    d, d2 = len(f) - 1, len(g) - 1
    if list(f) == list(g):
        return math.factorial(d)
    return math.factorial(d) * math.factorial(d2)


def minor_table(f: Poly, P: Sequence[Poly], g: Poly, Q: Sequence[Poly], r: int,
                prec0: int = 256) -> Dict[str, object]:
    """Certified zero / nonzero status of every r x r minor (unordered r-subsets of rows), and certified log|det|.

    Zero test (rigorous): a minor is an algebraic integer of the Galois closure L; its conjugates are +-(other minors),
    all bounded by Hmax := max Hadamard bound. If it is nonzero then |det| >= Hmax^{-(n-1)}, n = [L:Q] (norm >= 1).
    So at precision p with the ball upper bound below Hmax^{-(n-1)} the minor is exactly 0."""
    n = galois_degree_bound(f, g)
    prec = prec0
    while True:
        rows, d = rows_at(f, P, g, Q, prec)
        D = len(rows)
        old = ctx.prec
        ctx.prec = prec
        try:
            norms = []
            for row in rows:
                s = arb(0)
                for z in row:
                    s += abs(z) ** 2
                norms.append(s.sqrt())
            subsets = list(itertools.combinations(range(D), r))
            hmax = arb(1)
            for I in subsets:
                h = arb(1)
                for i in I:
                    h *= norms[i]
                hmax = hmax.max(h)
            log_sep = -(n - 1) * hmax.log()  # log of the separation bound
            status, logabs, undecided = {}, {}, 0
            for I in subsets:
                dt = _det([rows[i] for i in I])
                a = abs(dt)
                if a.upper() < log_sep.exp().lower():
                    status[I] = 0
                elif a.lower() > 0:
                    status[I] = 1
                    logabs[I] = a.log()
                else:
                    undecided += 1
                    break
        finally:
            ctx.prec = old
        if undecided == 0:
            return {"D": D, "d": d, "r": r, "status": status, "logabs": logabs, "prec": prec, "n": n}
        prec *= 2
        if prec > (1 << 18):
            raise RuntimeError("could not certify minors")


# ----------------------------------------------------------------------------- frames, marginals, coverage


def frames(table: Dict[str, object]) -> List[Tuple[int, ...]]:
    D, r, st = table["D"], table["r"], table["status"]
    out = []
    for S in itertools.combinations(range(D), r + 1):
        if all(st[tuple(x for x in S if x != i)] == 1 for i in S):
            out.append(S)
    return out


def log_phi(table, S) -> arb:
    t = arb(0)
    for i in S:
        t += table["logabs"][tuple(x for x in S if x != i)]
    return t


def cross_pairs(table) -> List[Tuple[int, int]]:
    d, D = table["d"], table["D"]
    return [(i, j) for i in range(d) for j in range(d, D)]


def distribution_stats(table, fr: List[Tuple[int, ...]], w: Sequence[Fr]) -> Dict[str, object]:
    """Point marginals P(i in S), per-cross-pair coverage P(i,j in S), and min coverage (exact Fractions)."""
    D = table["D"]
    marg = [Fr(0)] * D
    cov = {p: Fr(0) for p in cross_pairs(table)}
    for S, ws in zip(fr, w):
        if ws == 0:
            continue
        for i in S:
            marg[i] += ws
        Sset = set(S)
        for (i, j) in cov:
            if i in Sset and j in Sset:
                cov[(i, j)] += ws
    return {"marg": marg, "cov": cov, "mincov": min(cov.values())}


def pure_s(table, fr, s: int) -> List[Fr]:
    """Uniform distribution on the frames with exactly s K-rows (invariant under every side-preserving symmetry of
    the configuration, in particular under Galois). Empty -> all zeros."""
    d = table["d"]
    sel = [1 if sum(1 for i in S if i < d) == s else 0 for S in fr]
    k = sum(sel)
    return [Fr(x, k) if k else Fr(0) for x in sel]


def fh_log_bound(table, fr, w) -> arb:
    """log of 2^r C_r prod_S Phi_S^{-w_S/(r+1)} (certified ball)."""
    r = table["r"]
    old = ctx.prec
    ctx.prec = 256
    try:
        t = arb(2).log() * r + bean_thunder_C(r).log()
        for S, ws in zip(fr, w):
            if ws:
                t -= arb(ws.numerator) / ws.denominator * log_phi(table, S) / (r + 1)
        return t
    finally:
        ctx.prec = old


# ----------------------------------------------------------------------------- exact LP (dense simplex, Fractions)


def lp_max(c: List[Fr], A_eq: List[List[Fr]], b_eq: List[Fr], A_ub: List[List[Fr]], b_ub: List[Fr]):
    """maximise c.x s.t. A_eq x = b_eq, A_ub x <= b_ub, x >= 0. Two-phase simplex, Bland's rule, exact.
    Returns (value, x) or (None, None) if infeasible. Small dense problems only."""
    n = len(c)
    rows, rhs, basis = [], [], []
    # slack for <=, artificial for = (and for <= rows with negative rhs, not needed here: we require b_ub >= 0)
    m_ub, m_eq = len(A_ub), len(A_eq)
    ncol = n + m_ub + m_eq
    for k, (a, b) in enumerate(zip(A_ub, b_ub)):
        if b < 0:
            raise ValueError("b_ub must be >= 0")
        row = list(a) + [Fr(0)] * (m_ub + m_eq)
        row[n + k] = Fr(1)
        rows.append(row); rhs.append(Fr(b)); basis.append(n + k)
    for k, (a, b) in enumerate(zip(A_eq, b_eq)):
        a = list(a); b = Fr(b)
        if b < 0:
            a = [-x for x in a]; b = -b
        row = a + [Fr(0)] * (m_ub + m_eq)
        row[n + m_ub + k] = Fr(1)
        rows.append(row); rhs.append(b); basis.append(n + m_ub + k)
    art = set(range(n + m_ub, ncol))

    def pivot(pr, pc):
        pv = rows[pr][pc]
        rows[pr] = [x / pv for x in rows[pr]]; rhs[pr] /= pv
        for i in range(len(rows)):
            if i != pr and rows[i][pc] != 0:
                f_ = rows[i][pc]
                rows[i] = [x - f_ * y for x, y in zip(rows[i], rows[pr])]
                rhs[i] -= f_ * rhs[pr]
        basis[pr] = pc

    def run(obj, allowed):
        while True:
            # reduced costs
            cb = [obj[b] for b in basis]
            best = None
            for j in range(ncol):
                if j not in allowed or j in basis:
                    continue
                rc = obj[j] - sum(cb[i] * rows[i][j] for i in range(len(rows)))
                if rc > 0:
                    best = j
                    break
            if best is None:
                return
            ratios = [(rhs[i] / rows[i][best], basis[i], i) for i in range(len(rows)) if rows[i][best] > 0]
            if not ratios:
                raise RuntimeError("unbounded")
            _, _, pr = min(ratios)
            pivot(pr, best)

    obj1 = [Fr(0)] * ncol
    for j in art:
        obj1[j] = Fr(-1)
    run(obj1, set(range(ncol)))
    if any(basis[i] in art and rhs[i] != 0 for i in range(len(rows))):
        return None, None
    # drive artificials out of basis where possible
    for i in range(len(rows)):
        if basis[i] in art:
            for j in range(n + m_ub):
                if rows[i][j] != 0:
                    pivot(i, j)
                    break
    obj2 = list(c) + [Fr(0)] * (m_ub + m_eq)
    run(obj2, set(range(n + m_ub)))
    x = [Fr(0)] * ncol
    for i, b in enumerate(basis):
        x[b] = rhs[i]
    return sum(ci * xi for ci, xi in zip(c, x[:n])), x[:n]


def maxmin_coverage(table, fr, target_marg: Sequence[Fr]):
    """max_w min_{cross pairs} P_w(pair in S) subject to P_w(i in S) = target_marg[i], w a distribution on frames.
    Variables (w_1..w_F, t). The optimum is attained by a symmetry-invariant w (average an optimal w over the
    automorphism group: constraints and the concave objective are invariant), so the value is a valid FH exponent.
    Returns (t*, w*) or (None, None) if the marginals are infeasible."""
    D = table["D"]
    F = len(fr)
    pairs = cross_pairs(table)
    c = [Fr(0)] * F + [Fr(1)]
    A_eq, b_eq = [], []
    A_eq.append([Fr(1)] * F + [Fr(0)]); b_eq.append(Fr(1))
    for i in range(D):
        A_eq.append([Fr(1) if i in S else Fr(0) for S in fr] + [Fr(0)]); b_eq.append(Fr(target_marg[i]))
    A_ub, b_ub = [], []
    for (i, j) in pairs:
        A_ub.append([Fr(-1) if (i in S and j in S) else Fr(0) for S in fr] + [Fr(1)]); b_ub.append(Fr(0))
    val, x = lp_max(c, A_eq, b_eq, A_ub, b_ub)
    if val is None:
        return None, None
    return val, x[:F]


# ----------------------------------------------------------------------------- exact box counts


def box_norm_pairs(f, P, g, Q, r: int, E: int) -> List[Tuple[int, int]]:
    """(|W(a)|, |W'(a)|) for all a in [-E,E]^r with both nonzero (exact resultants)."""
    out = []
    for a in itertools.product(range(-E, E + 1), repeat=r):
        AK = poly_comb(P, a)
        AK2 = poly_comb(Q, a)
        if all(x == 0 for x in AK) or all(x == 0 for x in AK2):
            continue
        w1, w2 = abs(resultant(f, AK)), abs(resultant(g, AK2))
        if w1 and w2:
            out.append((w1, w2))
    return out
