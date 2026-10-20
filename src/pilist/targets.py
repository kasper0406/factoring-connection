"""Cost targets for the canonical product-listing problem Pi (agent pilist, Task 1).

ANALYSIS ONLY (AGENT_RULES rule 1): floating point is used to evaluate and minimise L-constant programs.

Two families of targets.

(A) Congruence form at modulus N (Dixon-shaped; covers nonexact Pi_beta and Pi-split).
    Level alpha, factor base K = L[alpha, kappa], window |v| <= V = L[2 alpha, nu]
    (i.e. log V = nu x^{2 alpha} (log x)^{1 - 2 alpha}), lister cost L[alpha, tau] per output.
    One-sided outputs: v is K-smooth with probability L[alpha, -alpha nu / kappa]  (u log u with
    u = (nu/kappa)(x/log x)^alpha, log u ~ alpha log x).  Total exponent:
        one-sided:     max(kappa + tau + alpha nu / kappa, 2 kappa)
        doubly split:  max(kappa + tau, 2 kappa)
    Sanity check: alpha = 1/2, nu = 1 (V = N), tau = 0 is Dixon, L[1/2, sqrt 2].

(B) NFS-shaped listers (heuristic H of notes/heurnfs_results.md; program of src/smoothenum/target.py):
    partial listers with cost E^2 L[-theta * penalty] on the rational side, the algebraic side, or both.
"""
from __future__ import annotations

from math import sqrt
from typing import Dict

import numpy as np

from src.heurnfs.optim import max_violation, minimise
from src.smoothenum.target import C_MNFS, _parts, c_rat_theta_closed


# ----------------------------------------------------------------------------------------------- (A)
def congruence_cost(alpha: float, kappa: float, nu: float, tau: float, doubly: bool = False) -> float:
    """L[alpha, .] exponent of a congruence-of-squares algorithm fed by a Pi lister (see module doc)."""
    coll = kappa + tau + (0.0 if doubly else alpha * nu / kappa)
    return max(coll, 2 * kappa)


def tau_max(alpha: float, c: float, nu: float, doubly: bool = False) -> float:
    """Largest per-output lister exponent tau that still gives total exponent <= c at level alpha.

    One-sided: max over 0 < kappa <= c/2 of c - kappa - alpha nu/kappa. The unconstrained maximiser is
    kappa* = sqrt(alpha nu); if kappa* <= c/2 the value is c - 2 sqrt(alpha nu), else (kappa = c/2)
    c/2 - 2 alpha nu / c.  Doubly split: c - kappa at the smallest useful kappa ... but kappa also has to
    pay LA, so tau_max = c - kappa for any kappa <= c/2: the sup is c (kappa -> 0), approached, not attained.
    """
    if doubly:
        return c
    ks = sqrt(alpha * nu)
    if ks <= c / 2:
        return c - 2 * ks
    return c / 2 - 2 * alpha * nu / c


def dixon_check() -> float:
    """min_kappa congruence_cost(1/2, kappa, 1, 0) = sqrt 2 (Dixon); returns the numeric minimum."""
    ks = np.linspace(0.05, 2.0, 20001)
    return float(min(congruence_cost(0.5, k, 1.0, 0.0) for k in ks))


def blind_cost_exponent_at_level(nu: float) -> str:
    """A blind (residue-uniform) search pays N/V per output: L[1, 1 - o(1)], i.e. exponential in L[alpha]
    units for every alpha < 1.  Returned as text for the tables."""
    return "N/V = L[1, 1-o(1)] per output"


# ----------------------------------------------------------------------------------------------- (B)
def nfs_costs(x, mode: str, theta: float) -> Dict[str, float]:
    """Partial listers in the MNFS program. mode in {'rat', 'alg', 'both'}; theta in [0, 1].

    rat : enumerate the rational-smooth pairs at E^2 V^0 L[-theta r]; test them in V fields (2e - r + v).
    alg : enumerate the algebraic-smooth pairs (all V fields) at E^2 V L[-theta a]; test the rational side
          of each (2e + v - a).
    both: enumerate the doubly-smooth pairs at E^2 V L[-theta (r + a)]; no test.
    theta = 0 is the box sieve (with the respective test term); theta = 1 is output-sensitive."""
    e, b, ba, v, d, r, a, s = _parts(x)
    if mode == "rat":
        enum, test = 2 * e - theta * r, 2 * e - r + v
    elif mode == "alg":
        enum, test = 2 * e + v - theta * a, 2 * e + v - a
    elif mode == "both":
        enum, test = 2 * e + v - theta * (r + a), -1e9
    else:
        raise ValueError(mode)
    return {"enum": enum, "test": test, "la": 2 * s, "slack": 2 * e + v - r - a - s, "r": r, "a": a}


def nfs_program(mode: str, theta: float, emax: float = 3.0, vmax: float = 0.6, restarts: int = 60,
                seed: int = 0) -> Dict[str, float]:
    def obj(x):
        c = nfs_costs(x, mode, theta)
        return max(c["enum"], c["test"], c["la"])

    def cons(x):
        return [nfs_costs(x, mode, theta)["slack"]]
    x, val = minimise(obj, cons, [.3, .3, .3, .001, .3], [emax, 1.3, 1.3, vmax, 3.0], restarts=restarts,
                      seed=seed)
    out = {k: float(t) for k, t in zip(("e", "beta", "beta_a", "v", "delta"), x)}
    out.update({k: float(t) for k, t in nfs_costs(x, mode, theta).items()})
    out["c"] = float(val)
    out["viol"] = max_violation(cons, x)
    return out


def nfs_theta_curve(mode: str, thetas=(0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0), emax: float = 3.0,
                    seed: int = 0) -> Dict[float, float]:
    """c(theta) for a partial lister; running minimum with C_MNFS (Nelder-Mead only overestimates)."""
    out, best = {}, C_MNFS
    for t in thetas:
        if mode == "rat":
            c = c_rat_theta_closed(t)
        else:
            c = nfs_program(mode, t, emax=emax, seed=seed)["c"]
        c = min(c, best)
        out[t] = c
        best = c
    return out


def doubly_per_output(tau: float, h: float = 1e-3) -> float:
    """Doubly-smooth lister with per-output cost L[1/3, tau] (any box): cost max(s + tau, 2 s), and the
    codex family beta = beta_a = h, delta = 3h, e = h + 2/(9h^2) makes the box hold L[h] doubly-smooth pairs
    (slack 0, s = h). So c = max(h + tau, 2h) -> tau as h -> 0."""
    from fractions import Fraction as Fr
    H = Fr(h).limit_denominator(10 ** 12)
    e, b, ba, d = H + Fr(2) / (9 * H * H), H, H, 3 * H
    r = 1 / (3 * b * d)
    a = (d * e + 1 / d) / (3 * ba)
    s = max(b, ba)
    assert 2 * e - r - a - s == 0          # exact: the box holds L[h] doubly-smooth pairs
    return max(float(s) + tau, 2 * float(s))


if __name__ == "__main__":  # pragma: no cover
    import json
    print("Dixon check", dixon_check())
    for al in (1 / 3, 1 / 4):
        for nu in (0.25, 0.5, 1.0, 2.0):
            print("alpha", round(al, 3), "nu", nu, "tau_max(1.9019)" if al > .3 else "tau_max(c=1)",
                  round(tau_max(al, C_MNFS if al > .3 else 1.0, nu), 4))
    for m in ("rat", "alg", "both"):
        print(m, json.dumps({k: round(v, 4) for k, v in nfs_theta_curve(m).items()}))
