"""Ledgers for degenerate and split-side linked modules (agent rankdegen). Analysis code: floats only in the numeric
cross-checks (labelled), closed forms in exact sympy/Fraction arithmetic.

Units (as rank3/resbarrier/heurnfs): x = log N, t = x^{1/3}(log x)^{2/3}; side degrees d = delta_1 (x/log x)^{1/3},
d' = delta_2 (x/log x)^{1/3}; side log-norms l_i = eta_i x^{2/3}(log x)^{1/3}; log B = beta t; count entropy gamma t.

Frame program (Theorem FH + linking, pure-s frame distributions; notes section 4):
  for each allowed s:  (r/(r+1)) (s p1 + (r+1-s) p2) >= gamma + c_r s (r+1-s)/(delta_1 delta_2),  c_r = (r-1)/(r+1),
  with p_i = eta_i/delta_i; projection constraints rho p1 >= gamma (rank-deficient K side, Lemma PR), and the
  one-sided frames s = 0, r+1 give r p2 >= gamma, r p1 >= gamma.
  H (both sides charged):  gamma >= beta + (eta_1+eta_2)/(3 beta);  cost >= max(gamma, 2 beta).
Reduction: q = p/gamma, R = gamma delta_1 delta_2; Phi := inf q1 q2 R over the constraints; then
  9 beta^2 (gamma-beta)^2 >= 4 gamma Phi, formal c^3 = 64 Phi/9, scoped (gamma <= r beta/(r-1)) c = 2(4r(r-1)Phi/9)^{1/3}.
Proposition PHI (proved, notes): Phi >= (r-1)/r whenever the allowed s form a contiguous range [s_lo, s_hi] in [1, r]
with s_hi q1 >= 1 and (r+1-s_lo) q2 >= 1 enforced; equality in the rational-side limit.
"""
from __future__ import annotations

import math
from fractions import Fraction as Fr
from typing import Iterable, Sequence

import sympy as sp


def c_r(r: int) -> Fr:
    return Fr(r - 1, r + 1)


def phi_closed(r: int) -> Fr:
    return Fr(r - 1, r)


def formal_constant(phi) -> sp.Expr:
    return sp.cbrt(sp.Rational(64) * sp.nsimplify(phi) / 9)


def scoped_constant(phi, r: int) -> sp.Expr:
    return 2 * sp.cbrt(sp.Rational(4 * r * (r - 1)) * sp.nsimplify(phi) / 9)


def phi_constraint_R(q1: float, q2: float, r: int, S: Iterable[int]) -> float:
    """Smallest R allowed by the pure-s constraints at (q1, q2) (inf if infeasible)."""
    cr = (r - 1) / (r + 1)
    R = 0.0
    for s in S:
        den = (r / (r + 1)) * (s * q1 + (r + 1 - s) * q2) - 1
        num = cr * s * (r + 1 - s)
        if num == 0:
            if den < -1e-12:
                return math.inf
            continue
        if den <= 0:
            return math.inf
        R = max(R, num / den)
    return R


def phi_numeric(r: int, S: Sequence[int], lo1: float = 0.0, lo2: float = 0.0, grid: int = 400,
                rounds: int = 6) -> float:
    """NUMERIC (floating point, analysis only): inf_{q1>=lo1, q2>=lo2} q1 q2 R(q). Log-grid + zoom."""
    import numpy as np
    cr = (r - 1) / (r + 1)
    a1, b1 = math.log(max(lo1, 1e-3)), math.log(1e5)
    a2, b2 = math.log(max(lo2, 1e-3)), math.log(1e5)
    best = math.inf
    arg = None
    for _ in range(rounds):
        q1 = np.exp(np.linspace(a1, b1, grid))[:, None]
        q2 = np.exp(np.linspace(a2, b2, grid))[None, :]
        R = np.zeros((grid, grid))
        bad = (q1 < lo1 - 1e-12) | (q2 < lo2 - 1e-12)
        for s in S:
            den = (r / (r + 1)) * (s * q1 + (r + 1 - s) * q2) - 1
            num = cr * s * (r + 1 - s)
            if num == 0:
                bad = bad | (den < -1e-12)
                continue
            bad = bad | (den <= 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                R = np.maximum(R, np.where(den > 0, num / den, np.inf))
        V = np.where(bad, np.inf, q1 * q2 * R)
        i, j = np.unravel_index(np.argmin(V), V.shape)
        if V[i, j] < best:
            best, arg = float(V[i, j]), (float(q1[i, 0]), float(q2[0, j]))
        w1, w2 = (b1 - a1) / grid * 4, (b2 - a2) / grid * 4
        l1, l2 = math.log(arg[0]), math.log(arg[1])
        a1, b1 = max(math.log(max(lo1, 1e-3)), l1 - w1), l1 + w1
        a2, b2 = max(math.log(max(lo2, 1e-3)), l2 - w2), l2 + w2
    return best


def phi_crossing_certificate(r: int, s_lo: int, s_hi: int, q1: Fr, q2: Fr) -> bool:
    """Exact check of the crossing lemma at a point: some s in [s_lo, s_hi] has (s q1 - 1)((r+1-s) q2 - 1) >= -1/r
    (which is equivalent to q1 q2 g_s(q) >= (r-1)/r for that s)."""
    return any((s * q1 - 1) * ((r + 1 - s) * q2 - 1) >= Fr(-1, r) for s in range(s_lo, s_hi + 1))


# ----------------------------------------------------------------------------- split side (ptesplit O3)


def split_side_charge(r: int, d: int, d2: int, l1: float, x: float, s: int) -> float:
    """N-part of the pure-s lower bound on the charged log-norm l2 (Corollary SPL):
    l2 >= (s/d)[(r-1)x/r - d' l1/(r+1-s)] + ((r+1) d'/(r (r+1-s))) (lambda - log(2^r C_r)).  Returns the N-part."""
    return (s / d) * ((r - 1) * x / r - d2 * l1 / (r + 1 - s))


def split_threshold(r: int, d2: int, x: float, s: int = 1) -> float:
    """l1 above which the pure-s bound no longer charges the N-part: (r-1)(r+1-s) x/(r d'); s = 1: (r-1)x/d'."""
    return (r - 1) * (r + 1 - s) * x / (r * d2)


def r5_threshold(r: int, D: int, x: float) -> float:
    return (r - 1) * x / (D - 1)


def hsplit_floor(r: int, scoped: bool = False, seed: int = 1, iters: int = 6000) -> float:
    """NUMERIC (analysis only). Floor of the O3 program under H_split (split-by-construction candidates occur among
    box points of given sizes with density rho(u_1)): count entropy gamma >= beta + (eta1+eta2)/(3 beta) (box points
    needed), cost = max(beta + eta2/(3 beta), 2 beta) (only split candidates are examined), GP two-sided frame
    constraints s = 0..r+1, scoped: gamma <= r beta/(r-1). Random search + local refinement over
    (beta, eta1, eta2, delta1, delta2); returns the smallest cost found (an upper estimate of the program's infimum)."""
    import random
    rnd = random.Random(seed)
    cr = (r - 1) / (r + 1)

    def feasible(beta, e1, e2, d1, d2):
        gamma = beta + (e1 + e2) / (3 * beta)
        if scoped and gamma > r * beta / (r - 1) + 1e-12:
            return False
        p1, p2 = e1 / d1, e2 / d2
        for s in range(0, r + 2):
            lhs = (r / (r + 1)) * (s * p1 + (r + 1 - s) * p2)
            if lhs < gamma + cr * s * (r + 1 - s) / (d1 * d2) - 1e-12:
                return False
        return True

    def cost(beta, e1, e2):
        return max(beta + e2 / (3 * beta), 2 * beta)

    best, x = math.inf, None
    for _ in range(iters):
        beta = rnd.uniform(0.05, 2.0)
        e1, e2 = math.exp(rnd.uniform(-4, 3)), math.exp(rnd.uniform(-4, 3))
        d1, d2 = math.exp(rnd.uniform(-4, 4)), math.exp(rnd.uniform(-4, 4))
        if feasible(beta, e1, e2, d1, d2):
            c = cost(beta, e1, e2)
            if c < best:
                best, x = c, [beta, e1, e2, d1, d2]
    if x is None:
        return math.inf
    step = 0.3
    for it in range(iters * 3):
        y = [v * math.exp(rnd.gauss(0, step)) for v in x]
        if feasible(*y):
            c = cost(y[0], y[1], y[2])
            if c < best:
                best, x = c, y
        if it % 1000 == 999:
            step *= 0.7
    return best
