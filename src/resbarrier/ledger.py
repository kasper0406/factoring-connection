"""Exact ledgers, LPs and exponent floors for the generalized resultant barrier (agent resbarrier).

Units: x = log N; heights h are in units of x. Constants are L_N[1/3, .] exponents (single N, RAM, omega = 2).
Floating point appears only in the numeric MNFS minimisation reused from src/heurnfs (analysis only).
See notes/resbarrier_results.md for statements, labels and proofs.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from itertools import combinations
from typing import Dict, List, Sequence, Tuple

import sympy as sp

from src.heurnfs.models import gnfs_closed, height_lp, mnfs


def rank2_constant(k) -> sp.Expr:
    """Cor RB2: rank-2 floor with N^k | disc(F): ledger log X >= k x/(2(D-1)) + D C/2, i.e. kappa = k/2 in
    the GNFS program, constant (32 kappa/9)^{1/3} = (16 k/9)^{1/3}.
    k = 2 (two-map / NFS pair): (32/9)^{1/3} = 1.5263 (= SNFS); k = 1 (ramified, Kirchner-type): 1.2114;
    k = 4 corresponds to the GNFS ledger 2x/D (heuristic G or Theorem GAE): (64/9)^{1/3} = 1.9230."""
    return gnfs_closed(sp.Rational(k, 2))


def alpha_floor(a: Fr, b: Fr) -> Fr:
    """Ledger log X >~ x^a C^b with C = x^alpha (up to logs). Smoothing L[alpha] needs (a + b alpha) - alpha
    <= alpha, i.e. alpha >= a/(2 - b). Generic (a,b) = (1/2,1/2): 1/3. Aggregate (1/2, 0): 1/4.
    Dixon-type (1, 0): 1/2."""
    return Fr(a) / (2 - Fr(b))


def lp_min(cons: Sequence[Tuple[Sequence[Fr], Fr]], n: int = 2) -> Fr:
    """Exact min sum(h) over h >= 0 subject to sum_i c_i h_i >= rhs for each (c, rhs), by vertex enumeration."""
    rows: List[Tuple[List[Fr], Fr]] = [([Fr(c) for c in cc], Fr(r)) for cc, r in cons]
    for i in range(n):  # h_i >= 0
        e = [Fr(0)] * n
        e[i] = Fr(1)
        rows.append((e, Fr(0)))
    best = None
    for sub in combinations(rows, n):
        M = sp.Matrix([[c for c in r[0]] for r in sub])
        if M.det() == 0:
            continue
        sol = M.LUsolve(sp.Matrix([r[1] for r in sub]))
        h = [Fr(int(sp.fraction(s)[0]), int(sp.fraction(s)[1])) for s in sol]
        if all(sum(c * hi for c, hi in zip(r[0], h)) >= r[1] for r in rows):
            v = sum(h)
            best = v if best is None or v < best else best
    return best


def res_only_min(d1: int, d2: int) -> Fr:
    """Lemma RES alone: d2 h1 + d1 h2 >= 1  =>  min h1 + h2 = 1/max(d1, d2) (SNFS attains it)."""
    return lp_min([((d2, d1), 1)])


def gae_min(d1: int, d2: int) -> Fr:
    """Theorem GAE (almost all N) + RES: adds (d1+1) h1 + (d2+1) h2 >= 1; equals heurnfs.height_lp."""
    return lp_min([((d2, d1), 1), ((d1 + 1, d2 + 1), 1)])


def gr_min(r: int, d1: int, d2: int) -> Fr:
    """HEURISTIC G_r (rank-r linked modules): Schmidt-type count M^{max(r, d+1)} of rank-r norm forms of degree
    d and height M, against r - 1 independent linking conditions mod N:
        max(r, d1+1) h1 + max(r, d2+1) h2 >= r - 1,
    together with the rigorous 2-plane resultant constraint d2 h1 + d1 h2 >= 1."""
    return lp_min([((max(r, d1 + 1), max(r, d2 + 1)), r - 1), ((d2, d1), 1)])


def kl_min(r: int, d2: int) -> Fr:
    """Lemma KL (rigorous; rational side W1 = sum a_j v_j, nonlinear side of degree d2, rank r <= d2 + 1):
        h2 + (d2/(r-1)) h1 >= 1   =>   min h1 + h2 = min(1, (r-1)/d2)."""
    return lp_min([((Fr(d2, r - 1), 1), 1)])


def mnfs_res_only_closed() -> sp.Expr:
    """Exact optimum of the RES-only MNFS program (kr = 1, ka = 1/2): c^3 = min_t 32(1/2+t)/(9 t^2 (2-t)^2)
    = (36 + 16 sqrt 6)/15 at t = beta_a/beta = sqrt(2/3) (codex resbarrier_2 #9). Numerically 1.7114326...
    This coincides with the KNOWN medium-prime MexTNFS constant (Kim-Barbulescu 2016; Sarkar-Singh 2016):
    the same Coppersmith-type program with a halved algebraic height term. Not new."""
    return sp.cbrt((36 + 16 * sp.sqrt(6)) / 15)


def mnfs_floors(seed: int = 0) -> Dict[str, float]:
    """MNFS constants (numeric, heurnfs program) under (i) Lemma RES only: rational multiplier kr = 1 and,
    from pairwise resultants among the V algebraic sides, ka >= 1/2; (ii) G / Theorem GAE: kr = ka = 1."""
    return {"res_only": mnfs(kr=1.0, ka=0.5, seed=seed)["c"], "gae": mnfs(kr=1.0, ka=1.0, seed=seed)["c"]}


def gr_expected_log_count_upper(r: int, d: int, h_total: Fr) -> Fr:
    """Same with Schmidt's UPPER count M^{d+r-1}: (d+r-1) h_total - (r-1) (the exponent relevant for nonexistence)."""
    return (d + r - 1) * Fr(h_total) - (r - 1)


def gr_expected_log_count(r: int, d: int, h_total: Fr) -> Fr:
    """log_N of the G_r expected number of linked rank-r pairs (both sides degree d) at total height h_total:
    max(r, d+1) * h_total - (r - 1). Aggregate families need h_total ~ 1/d: value (d+1)/d - (r-1) < 0, r >= 3."""
    return max(r, d + 1) * Fr(h_total) - (r - 1)
