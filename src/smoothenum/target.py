"""Task 1: which enumeration cost gives which L[1/3] constant (agent smoothenum).

Units as in src/heurnfs/models.py: x = log N; every cost is an L_N[1/3, .] exponent; d = delta (x/log x)^{1/3};
box E = L[e]; rational FB B = L[beta]; algebraic FB (each of V = L[v] fields) L[beta_a].
  r = 1/(3 beta delta)                rational smoothness penalty (rational norm ~ N^{1/d}, size 1/delta)
  a = (delta e + 1/delta)/(3 beta_a)  algebraic smoothness penalty (norm ~ E^d N^{1/d})
  s = max(beta, beta_a + v)           matrix columns
Relation count (HEURISTIC H): 2e + v - r - a >= s.   LA: 2s (block Wiedemann, omega = 2).

Enumeration modes (the only thing this module varies):
  'box'      : sieve the box, cost 2e - sigma (sigma = an absolute L-exponent gain; known tricks give 0)
  'rat'      : list the rational-smooth pairs at cost 2e - theta*r  (theta = 1: output-sensitive)
  'alg'      : list, per field, the algebraic-smooth pairs at output cost 2e - a (V fields: + v),
               then test the rational side on each
  'both'     : list the doubly-smooth pairs at output cost 2e + v - r - a
Candidate tests: modes box/rat test every rational-smooth pair in V fields (2e - r + v); 'alg' tests the
rational side of every algebraic-smooth pair (2e + v - a); 'both' has no test term.

All numeric optima: seeded Nelder-Mead (src/heurnfs/optim.py) = analysis only (AGENT_RULES rule 1).
"""
from __future__ import annotations

from typing import Dict

import numpy as np
import sympy as sp

from src.heurnfs.optim import max_violation, minimise

C_MNFS = float(sp.cbrt((92 + 26 * sp.sqrt(13)) / sp.Integer(27)))
C_GNFS = float(sp.cbrt(sp.Rational(64, 9)))


def factory_closed() -> sp.Expr:
    """Per-N factory constant = MNFS program with output-sensitive rational enumeration (v* = 0).

    Derivation (model-PROVED). Binding: 2e - r = 2 beta = c, a = beta, beta_a = beta, v = 0.
    e = beta + 1/(6 beta delta); a <= beta  <=>  delta beta + 1/(6 beta) + 1/delta <= 3 beta^2; the left side is
    minimal at delta = beta^{-1/2}: 2 sqrt(beta) + 1/(6 beta) = 3 beta^2. With y = beta^{3/2}:
    18 y^2 - 12 y - 1 = 0, y = (2 + sqrt 6)/6.  c = 2 ((2+sqrt6)/6)^{2/3},  c^3 = (20 + 8 sqrt 6)/9.
    Precomputation 2e = 2 beta + 1/(3 sqrt beta)."""
    y = (2 + sp.sqrt(6)) / 6
    return 2 * y ** sp.Rational(2, 3)


def factory_precomp_closed() -> sp.Expr:
    b = ((2 + sp.sqrt(6)) / 6) ** sp.Rational(2, 3)
    return 2 * b + 1 / (3 * sp.sqrt(b))


def alg_output_infimum() -> sp.Expr:
    """HYPOTHETICAL algebraic-side output-sensitive listing: infimum 2 * 18^{-1/3} = 0.7631 (not attained).

    Cost max(2s, s + r) with r <= s, r = 1/(3 beta delta), and growth of the algebraic-smooth count in e
    needs delta < 6 beta_a (2e - a = 2e(1 - delta/(6 beta_a)) - ...). Hence 1/(18 beta beta_a) < s <= ...,
    s^3 > 1/18, c = 2s -> 2/18^{1/3} as delta -> 6 beta, e -> infinity. V is irrelevant (v* = 0)."""
    return 2 / sp.cbrt(18)


def _parts(x):
    e, b, ba, v, d = x
    r = 1 / (3 * b * d)
    a = (d * e + 1 / d) / (3 * ba)
    s = max(b, ba + v)
    return e, b, ba, v, d, r, a, s


def costs(x, mode: str = "box", theta: float = 0.0, sigma: float = 0.0) -> Dict[str, float]:
    e, b, ba, v, d, r, a, s = _parts(x)
    if mode == "box":
        enum, test = 2 * e - sigma, 2 * e - r + v
    elif mode == "rat":
        enum, test = 2 * e - theta * r, 2 * e - r + v
    elif mode == "alg":
        enum, test = 2 * e + v - a, 2 * e + v - a
    elif mode == "both":
        enum, test = 2 * e + v - r - a, -1e9
    else:
        raise ValueError(mode)
    return {"enum": enum, "test": test, "la": 2 * s, "yield_slack": 2 * e + v - r - a - s, "r": r, "a": a}


def program(mode: str = "box", theta: float = 0.0, sigma: float = 0.0, emax: float = 1.4,
            vmax: float = 0.6, restarts: int = 60, seed: int = 0) -> Dict[str, float]:
    """min max(enum, test, LA) s.t. relation count, over (e, beta, beta_a, v, delta) in a box (e <= emax)."""
    def obj(x):
        c = costs(x, mode, theta, sigma)
        return max(c["enum"], c["test"], c["la"])

    def cons(x):
        return [costs(x, mode, theta, sigma)["yield_slack"], sigma - 0.0 if mode != "box" else
                _parts(x)[5] - sigma]  # in 'box' mode the gain sigma cannot exceed the penalty r
    x, val = minimise(obj, cons, [.3, .3, .3, .001, .3], [emax, 1.3, 1.3, vmax, 3.0],
                      restarts=restarts, seed=seed)
    x = np.minimum(x, [emax, 9, 9, vmax, 9])
    out = {k: float(t) for k, t in zip(("e", "beta", "beta_a", "v", "delta"), x)}
    out.update({k: float(t) for k, t in costs(x, mode, theta, sigma).items()})
    out["c"] = float(val)
    out["viol"] = max_violation(cons, x)
    return out


def theta_curve(thetas=(0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0), seed: int = 0) -> Dict[float, float]:
    """c(theta): rational listing at cost E^2 L[-theta r]. MNFS <= GNFS and c(theta) monotone are enforced by
    taking running minima (Nelder-Mead artifacts only ever overestimate)."""
    out, best = {}, C_MNFS
    for t in thetas:
        c = min(program("rat", theta=t, seed=seed)["c"], best)
        out[t] = c
        best = c
    return out


def sigma_curve(sigmas=(0.0, 0.01, 0.02, 0.05, 0.1, 0.2), seed: int = 0) -> Dict[float, float]:
    """c(sigma): box enumeration cost 2e - sigma (an absolute L[1/3, sigma] saving on the sieve)."""
    out, best = {}, C_MNFS
    for s_ in sigmas:
        c = min(program("box", sigma=s_, seed=seed)["c"], best)
        out[s_] = c
        best = c
    return out


def degree_one_both_sides_polytime(kappa: float = 2.0, K: float = 1.0) -> Dict[str, float]:
    """HYPOTHETICAL doubly-smooth output-sensitive listing for two LINEAR forms of height N^{1/2} (Gaussian
    heuristic) with B = (log N)^kappa and box E = N^K: log(#doubly smooth)/x ~ 2K - (2K+1)/kappa (to leading
    order, log u ~ log x). Positive for K > 1/(2(kappa-1)), so B^2 = (log N)^{2 kappa} relations exist and the
    algorithm would be polynomial time. Returns the leading-order exponent (units of x = log N)."""
    return {"log_count_over_x": 2 * K - (2 * K + 1) / kappa, "needs_K_above": 1 / (2 * (kappa - 1))}


def table(seed: int = 0) -> Dict[str, Dict]:
    out = {
        "mnfs_box": program("box", seed=seed),
        "rat_output_theta1": program("rat", theta=1.0, seed=seed),
        "alg_output_emax3": program("alg", emax=3.0, seed=seed),
        "alg_output_emax10": program("alg", emax=10.0, seed=seed),
        "both_output_emax3": program("both", emax=3.0, seed=seed),
    }
    return out


if __name__ == "__main__":  # pragma: no cover
    import json
    print("closed: MNFS", C_MNFS, "factory", float(factory_closed()), "precomp", float(factory_precomp_closed()),
          "alg-inf", float(alg_output_infimum()))
    for k, v in table().items():
        print(f"{k:22s}", json.dumps({a: round(b, 4) for a, b in v.items()}))
    print("theta curve", {k: round(v, 4) for k, v in theta_curve().items()})
    print("sigma curve", {k: round(v, 4) for k, v in sigma_curve().items()})
