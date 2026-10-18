"""Exact (Fraction) exponent bookkeeping for notes/ecmavg_results.md.

Units: an exponent x means N^x up to N^{o(1)} (polylog factors and constants c(b) dropped).
delta = N^{-theta} is the density of the subset A of S_n handed to a fallback.

Proved inputs (cited, not re-proved here):
  * Theorem R (notes/dchworst_results.md §3): DCL-wc with D''s design at height H = N^h
    factors every pair covered by some block; cost N^{1/6+o(1)} + list cost N^{1/6+5h/3}.
  * sieve2d Theorem 2: uncovered fraction <= C n^4/H^4 + C n^6 H^{1/3} N^{-1/6}, valid for
    2 <= H <= N^{1/15}/C0; Proposition 4: <= C n/H^2.  So the failure exponent of one
    DCL-wc run is  fail_R(h) = min(-2h, max(-4h, h/3 - 1/6)),  0 < h <= 1/15.
  * Harvey--Hittmeir (HH): worst case N^{1/5}.
  * Theorem AV ladder (D' + detector) with sieve2d's n^4/H^4 for E_b: per-rung E_c and
    detector-coincidence sets of exponent max(1/6 + 14h/3 - t, 3t - 1, 10h/3 - 1/6) that pay
    the worst case W = N^{1/5}; h <= 1/50 (Proposition U).

Nothing here is an experiment.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from typing import Callable, Dict, List, Optional, Sequence, Tuple

SIXTH, FIFTH = F(1, 6), F(1, 5)
H_MAX_R = F(1, 15)      # sieve2d Theorem 2 range (H <= N^{1/15}/C0)
H_MAX_AV = F(1, 50)     # avgtime Proposition U range


# ---------------------------------------------------------------------------
# 1. One DCL-wc run (Theorem R) and the fallback on an arbitrary subset A
# ---------------------------------------------------------------------------
def fail_R(h: F) -> F:
    """Exponent of the fraction of S_n failing DCL-wc at height N^h (Theorem R + sieve2d)."""
    assert 0 < h <= H_MAX_R
    return min(-2 * h, max(-4 * h, h / 3 - SIXTH))


def cost_R(h: F) -> F:
    """Cost exponent of one DCL-wc run at height N^h (base selection 1/6 is dominated)."""
    return SIXTH + F(5, 3) * h


def ladder_R(theta: F, h0: F, hmax: Optional[F] = None, hh_after: bool = True) -> F:
    """Exponent bound for E[1_A * cost] of: geometric ladder of DCL-wc runs at heights
    h0 <= h <= hmax (ratio O(1) in H, O(log n) rungs), then HH, on a subset A of density
    N^{-theta}.  Rung at height h is reached by <= min(delta, fail(previous rung)) of S_n; the
    previous rung has the same exponent h up to o(1) (constant ratio), so the contribution
    is cost_R(h) - max(theta, -fail_R(h)) for h > h0, and cost_R(h0) - theta for the first.
    All terms piecewise linear: the sup over (h0, hmax] is taken at breakpoints."""
    if hmax is None:
        hmax = h0
    assert 0 < h0 <= hmax <= H_MAX_R
    terms = [cost_R(h0) - theta]
    pts = {hmax} | {x for x in (F(1, 26), F(1, 14), theta / 4, theta / 2, F(1, 2) - 3 * theta)
                     if h0 < x < hmax}
    for h in pts:
        if h > h0:
            terms.append(cost_R(h) - max(theta, -fail_R(h)))
    if hh_after:
        terms.append(FIFTH - max(theta, -fail_R(hmax)))
    return max(terms)


def fallback_R(theta: F) -> F:
    """Closed form of the best fallback exponent on a density-N^{-theta} subset A:
         min( 1/5 - theta ,  max(1/6 - 7 theta/12, 1/5 - 2/13) )   for 0 <= theta <= 2/13,
         1/5 - theta                                               for theta >= 2/13,
    floored at 0.  The first option is HH alone on A; the second is the DCL-wc ladder started
    at H0 = N^{theta/4} (so its first rung already has fail <= delta) up to H = N^{1/26}
    (where the sieve2d bound bottoms out at N^{-2/13}), then HH.  (Proved given the cited
    theorems; minimality is within this family of ladders: `fallback_R_search`.)"""
    theta = F(theta)
    hh = FIFTH - theta
    if theta <= F(2, 13):
        lad = max(SIXTH - F(7, 12) * theta, FIFTH - F(2, 13))
        return max(F(0), min(hh, lad))
    return max(F(0), hh)


def fallback_R_search(theta: F, grid: int = 78) -> F:
    """Brute-force min over ladders (h0, hmax) on a rational grid, and HH alone."""
    best = max(F(0), FIFTH - F(theta))
    hs = [H_MAX_R * F(i, grid) for i in range(1, grid + 1)] + [F(theta) / 4, F(1, 26)]
    hs = sorted({h for h in hs if 0 < h <= H_MAX_R})
    for h0 in hs:
        for hmax in hs:
            if hmax >= h0:
                best = min(best, ladder_R(F(theta), h0, hmax))
    return best


# ---------------------------------------------------------------------------
# 2. codex r18 (*) : the AV ladder (D' + detector) with sieve2d's n^4/H^4
# ---------------------------------------------------------------------------
def fallback_AV_floor() -> Tuple[F, F, F]:
    """Smallest z such that the AV-ladder top terms are <= z:
       HH on E_b(H_max):       1/5 - 4h <= z
       E_c u C(T) pay N^{1/5}: 1/5 + max(1/6 + 14h/3 - t, 3t - 1, 10h/3 - 1/6) <= z
    Solve: from 3t - 1 <= z - 1/5 and 1/6 + 14h/3 - t <= z - 1/5 we need
    z >= 3/40 + 7h/2; with h >= (1/5 - z)/4 this gives z >= 2/15 at h = 1/60, t = 14/45.
    Returns (z, h, t)."""
    z = F(2, 15)
    h = (FIFTH - z) / 4
    t = (z + F(4, 5)) / 3
    assert FIFTH - 4 * h <= z
    assert FIFTH + max(SIXTH + F(14, 3) * h - t, 3 * t - 1, F(10, 3) * h - SIXTH) <= z
    assert h <= H_MAX_AV
    return z, h, t


def fallback_AV(theta: F) -> F:
    """codex r18 (*): max(1/6 - 7 theta/12, 2/15) (ladder started at N^{theta/4}); the
    detector on A costs t/2 - theta = 7/45 - theta, binding only for theta < 1/45."""
    theta = F(theta)
    z, h, t = fallback_AV_floor()
    return max(SIXTH - F(7, 12) * theta, z, t / 2 - theta)


# ---------------------------------------------------------------------------
# 3. ECM prefilter cost models and total exponent
# ---------------------------------------------------------------------------
def ecm_cost(b: F, theta: F, kappa: F = F(1), stopping: bool = True) -> F:
    """Prefilter cost exponent with B = N^b, J_max = N^theta curves, under a failure law
    P(first j curves fail) <= j^{-kappa} N^{o(1)} for every j <= J_max.
      stopping=False (codex r18): every input runs all J curves: b + theta.
      stopping=True: run curves in order, stop at the first proper divisor.  Expected cost
      B * sum_{j<=J} P(fail after j-1) = B * sum_j min(1, j^{-kappa}):
        kappa < 1: b + (1-kappa) theta;  kappa >= 1: b (kappa = 1 costs a log factor)."""
    if not stopping:
        return F(b) + F(theta)
    return F(b) + max(F(0), 1 - F(kappa)) * F(theta)


def total(b: F, theta: F, kappa: F = F(1), stopping: bool = True,
          fallback: Callable[[F], F] = fallback_R) -> F:
    """Average exponent: max(prefilter cost, fallback on the failure set of density
    N^{-kappa theta})."""
    return max(ecm_cost(b, theta, kappa, stopping), fallback(F(kappa) * F(theta)))


def _argmin_piecewise(g: Callable[[F], F], lines: Sequence[Tuple[F, F]],
                      lo: F, hi: F) -> Tuple[F, F]:
    """Minimise a continuous piecewise-linear g on [lo, hi] whose pieces are among `lines`
    (a + c x): the min is at an endpoint or a pairwise intersection."""
    cands = {F(lo), F(hi)}
    for (a1, c1), (a2, c2) in combinations(lines, 2):
        if c1 != c2:
            x = (a2 - a1) / (c1 - c2)
            if lo <= x <= hi:
                cands.add(x)
    best = min(cands, key=lambda x: (g(x), x))
    return g(best), best


def optimise(b: F, kappa: F = F(1), stopping: bool = True, theta_max: F = FIFTH,
             fallback: str = "R") -> Tuple[F, F]:
    """min over theta in [0, theta_max] of total(b, theta, ...). Returns (exponent, theta)."""
    b, kappa = F(b), F(kappa)
    fb = fallback_R if fallback == "R" else fallback_AV
    g = lambda th: total(b, th, kappa, stopping, fb)
    slope = 1 if not stopping else max(F(0), 1 - kappa)
    lines = [(b, F(slope)), (FIFTH, -kappa), (SIXTH, -F(7, 12) * kappa),
             (FIFTH - F(2, 13), F(0)), (F(0), F(0)),
             (F(2, 15), F(0)), (F(7, 45), -kappa)]
    return _argmin_piecewise(g, lines, F(0), F(theta_max))


def codex_point() -> Dict[str, F]:
    """r18 §5's design: B = N^{1/8}, J = N^{1/38}, no stopping, AV fallback."""
    b, th = F(1, 8), F(1, 38)
    return {"ecm": ecm_cost(b, th, stopping=False), "fallback": fallback_AV(th),
            "total": total(b, th, stopping=False, fallback=fallback_AV)}


def weak_kappa_curve(kappa: F, theta_max: F = FIFTH, b: F = F(0)) -> Tuple[F, F]:
    """Best exponent (b -> 0 limit when b = 0) under failure <= J^{-kappa}, stopping, R fallback."""
    return optimise(b, kappa, True, theta_max, "R")


def kappa_closed_form(kappa: F) -> F:
    """b -> 0, theta unconstrained (<= 1/5), stopping, R fallback, 0 < kappa < 1:
    balance (1-kappa) theta against the fallback F_R(kappa theta).
      ladder branch  (1-kappa) th = 1/6 - 7 kappa th/12  ->  2(1-kappa)/(12 - 5 kappa)
      HH branch      (1-kappa) th = 1/5 - kappa th       ->  (1-kappa)/5
    The fallback is the min of the two branches, so the value is
    min(2(1-kappa)/(12-5kappa), (1-kappa)/5) floored by 3/65 when the ladder branch is used."""
    kappa = F(kappa)
    lad = 2 * (1 - kappa) / (12 - 5 * kappa)
    th_l = lad / (1 - kappa)
    if kappa * th_l <= F(2, 13):
        lad = max(lad, FIFTH - F(2, 13))
    else:
        lad = F(1)
    hh = (1 - kappa) / 5
    return min(lad, hh)
