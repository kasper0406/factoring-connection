"""Cost formulas for agent `siconstruct` (exact sympy where possible; mpmath only to evaluate closed forms).

Notation: x = ln N, F_a = x^a (ln x)^{1-a}, L[a, c] = exp((c + o(1)) F_a).  Factor base y = L[1/3, kappa],
window h = L[2/3, nu], per-output constructor cost L[1/3, tau].

* `si_cost(tau)`          : min_kappa max(kappa + tau + 2/(9 kappa^2), 2 kappa)   (Prop SI with an L[1/3, tau] oracle)
* `si_cost_slack(t)`      : min_kappa max(kappa + 2/(9 kappa^2) + t/(3 kappa), 2 kappa)
                            (residue inflated to h * L[2/3, t]: a CVP oracle with slack t, Prop CV(d))
* `mitm_per_output(...)`  : amortised MITM per-output cost 2 sqrt(N / (R (h+1))) (Lemma MA)
* `cvp_bounds(...)`       : the exact inequalities of Prop CV, evaluated at high precision for tests
"""
from __future__ import annotations

from math import isqrt
from typing import Dict, Sequence

import mpmath as mp
import sympy as sp

K = sp.Symbol("kappa", positive=True)


def _min_max(trial: sp.Expr) -> Dict[str, sp.Expr]:
    """min over kappa > 0 of max(trial(kappa), 2 kappa), trial convex decreasing-then-increasing."""
    crit = [r for r in sp.solve(sp.Eq(sp.diff(trial, K), 0), K) if r.is_real and r > 0]
    kmin = crit[0]
    inter = [r for r in sp.solve(sp.Eq(trial, 2 * K), K) if r.is_real and r > 0]
    # the max is minimised at the intersection if it lies left of the trial minimiser, else at the minimiser
    left = [r for r in inter if sp.N(r - kmin) <= 0]
    if left:
        k0 = max(left, key=lambda r: sp.N(r))
        return {"kappa": k0, "c": 2 * k0, "regime": "intersection"}
    return {"kappa": kmin, "c": trial.subs(K, kmin), "regime": "trial-minimum"}


def si_cost(tau) -> Dict[str, sp.Expr]:
    tau = sp.nsimplify(tau)
    return _min_max(K + tau + sp.Rational(2, 9) / K ** 2)


def si_cost_slack(t) -> Dict[str, sp.Expr]:
    t = sp.nsimplify(t)
    return _min_max(K + sp.Rational(2, 9) / K ** 2 + t / (3 * K))


def tau_threshold(target) -> sp.Expr:
    """Largest per-output tau with si_cost(tau) < target (trial-minimum regime): (3/2)(4/9)^{1/3} + tau = target."""
    return sp.nsimplify(target) - sp.Rational(3, 2) * sp.Rational(4, 9) ** sp.Rational(1, 3)


def linear_oracle_cost() -> sp.Expr:
    """tau = kappa (per-output cost linear in the factor-base size): min_kappa 2 kappa + 2/(9 kappa^2)."""
    tr = 2 * K + sp.Rational(2, 9) / K ** 2
    k0 = [r for r in sp.solve(sp.diff(tr, K), K) if r.is_real and r > 0][0]
    return sp.simplify(tr.subs(K, k0))


# ------------------------------------------------------------------------------------------- Lemma MA (amortised MITM)
def mitm_expected_hits(nA: int, nB: int, h: int, N: int) -> sp.Rational:
    """Heuristic H_eq: a*b mod N uniform for N-oblivious lists, for every product range [N, KN]."""
    return sp.Rational(nA * nB * (h + 1), N)


def mitm_per_output(N: int, h: int, R: int) -> float:
    """Optimal amortised MITM per-output cost (|A| = |B| = sqrt(R N/(h+1))): 2 sqrt(N/(R (h+1)))."""
    return float(2 * mp.sqrt(mp.mpf(N) / (R * (h + 1))))


def mitm_list_size(N: int, h: int, R: int) -> int:
    return isqrt(R * N // (h + 1)) + 1


# ------------------------------------------------------------------------------------------- Prop CV (one-sided CVP)
def cvp_bounds(e: Sequence[int], P: Sequence[int], N: int, k: int, h: int, dps: int = 80) -> Dict[str, mp.mpf]:
    """For the Adleman-weight lattice b(e) = (sqrt(l_i) e_i, C <e, l>), target t_k = (0, C ln kN), C = kN/h:
    returns Q, Delta, dist2 = Q + C^2 Delta^2, the lower bound ln kN - 1/(4C^2) (Prop CV(a)) and the witness
    bound ln kN + h/kN + 1 (Prop CV(b), valid for squarefree u in I_k)."""
    with mp.workdps(dps):
        l = [mp.log(p) for p in P]
        C = mp.mpf(k * N) / h
        dot = mp.fsum(ei * li for ei, li in zip(e, l))
        Delta = dot - mp.log(k * N)
        Q = mp.fsum(ei * ei * li for ei, li in zip(e, l))
        return {"Q": Q, "Delta": Delta, "dist2": Q + C * C * Delta * Delta,
                "lower": mp.log(k * N) - 1 / (4 * C * C),
                "witness_upper": mp.log(k * N) + mp.mpf(h) / (k * N) + 1}


def cvp_center_excess(e: Sequence[int], P: Sequence[int], X: int, h: int, dps: int = 60) -> mp.mpf:
    """Prop CV(c*) (repaired converse, after codex siconstruct_2): lattice with C = X/h and target C·ln T,
    ln T = m + 1/(2C^2), m = (ln X + ln(X+h))/2.  Returns D_T^2 - (m + 1/(4C^2)), which equals the exact identity
    2 ln v + E_rep + C^2 (ln(u/v) - m)^2 with E_rep = sum (e_i^2 - |e_i|) ln p_i >= 0."""
    with mp.workdps(dps):
        l = [mp.log(p) for p in P]
        C = mp.mpf(X) / h
        m = (mp.log(X) + mp.log(X + h)) / 2
        lnT = m + 1 / (2 * C * C)
        dot = mp.fsum(ei * li for ei, li in zip(e, l))
        Q = mp.fsum(ei * ei * li for ei, li in zip(e, l))
        return Q + C * C * (dot - lnT) ** 2 - (m + 1 / (4 * C * C))


def cvp_center_exact(P: Sequence[int], X: int, h: int):
    """Exact CVP for the centred target by exhaustive enumeration (tiny P only): every vector with D_T^2 <= m + 2
    has |e_i| <= sqrt((m + 2)/ln p_i), so the box below is complete.  Returns (min excess, minimising e)."""
    import itertools
    m = (mp.log(X) + mp.log(X + h)) / 2
    bounds = [int(mp.sqrt((m + 2) / mp.log(p))) for p in P]
    best = None
    for e in itertools.product(*[range(-b, b + 1) for b in bounds]):
        ex = cvp_center_excess(e, P, X, h, dps=30)
        if best is None or ex < best[0]:
            best = (ex, e)
    return best
