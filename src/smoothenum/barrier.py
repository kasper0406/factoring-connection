"""Barriers for output-sensitive enumeration (agent smoothenum). ANALYSIS ONLY (floats allowed: rho, fits).

Model LQ (lattice-query). The algorithm's only access to the arithmetic of the box is:
  (1) set up a lattice Lambda_{M,R} = {(a, b) : a == R b (mod M)}, M a product of factor-base primes/ideals
      with known factorization (cost >= 1 per setup, even when the lattice has no box point);
  (2) enumerate its box points (cost 1 per point);  (3) test a point (cost L^{o(1)}).
Heuristic H_LQ: for a point of Lambda_M the cofactor G/M behaves like a uniform integer of size |G|/M,
independently of how M was chosen. Then a lattice of index M yields, per unit cost,
      Y(M) <= min(1, E^2/M) * P(y/M)        (P = smoothness probability, y = norm size),
and the best ratio Y(M)/P(y) over all M is the "LQ gain" g*. Proposition LQ: g* = L^{o(1)}, so in the L[1/3]
program the LQ model has theta = 0 (no change from box sieving).

Model GEN (bottom-up generation). The algorithm creates integers s with known factorization (cost >= 1 each)
and maps each to the unique box pair with G(a, b) = s (if any). A residue-blind generator hits the box with
probability (#box) / (2 * max|G|) ~ E/m (exact density statement: G is injective on the box once |Y0| > 2A).
A generator that fixes a factor s1 and lets s2 range over the residues making s1*s2 land in the box is
*exactly* an LQ query with M = s1 (the admissible s2 are the box points of Lambda_{s1}).
"""
from __future__ import annotations

from functools import lru_cache
from math import exp, log, log2
from typing import Dict, List

import numpy as np


@lru_cache(maxsize=None)
def _rho_table(umax: float = 60.0, h: float = 1e-3):
    """Dickman rho from the integral equation u rho(u) = int_{u-1}^{u} rho(t) dt (trapezoid, running sum)."""
    k = int(round(1 / h))
    n = int(umax / h) + 1
    us = np.arange(n) * h
    r = np.ones(n)
    for i in range(k + 1, n):  # exact re-summation each step (no cancellation when rho is tiny)
        S = float(np.sum(r[i - k + 1:i]))
        r[i] = h * (0.5 * r[i - k] + S) / (us[i] - 0.5 * h)
    return us, np.log(r)


def log_rho(u: float) -> float:
    """log of the Dickman function (ANALYSIS). Uses a table to u = 60, then the de Bruijn asymptotic."""
    if u <= 1:
        return 0.0
    us, lr = _rho_table()
    if u < us[-1]:
        return float(np.interp(u, us, lr))
    # rho(u) = exp(-u (log u + log log u - 1 + o(1)))
    return -u * (log(u) + log(log(u)) - 1 + (log(log(u)) - 1) / log(u))


def rho(u: float) -> float:
    return exp(log_rho(u))


def lq_gain(u: float, v: float) -> float:
    """rho(u - v)/rho(u): the per-point yield gain of a lattice of index B^v (v = log M/log B)."""
    return exp(log_rho(max(u - v, 0.0)) - log_rho(u))


def lq_best_at(u: float, v_box: float, logB: float) -> Dict[str, float]:
    """Best yield-per-cost gain over all lattice indices M = B^v (natural logs):
    gain(v) = min(1, exp((v_box - v) logB)) * rho(u - v)/rho(u). Returns argmax v and gain."""
    vs = np.union1d(np.linspace(0.0, min(u, v_box + 50.0), 5001), [min(v_box, u)])  # include v_box exactly
    best = (0.0, 1.0)
    for v in vs:
        g = log_rho(max(u - v, 0)) - log_rho(u) + min(0.0, (v_box - v) * logB)
        if g > log(best[1]):
            best = (float(v), exp(g))
    return {"v": best[0], "gain": best[1], "gain_at_box": lq_gain(u, v_box)}


def mnfs_point(bits: int, e: float = 0.95094, beta: float = 0.95094, delta: float = 1.40175,
               beta_a: float = 0.82591) -> Dict[str, float]:
    """Concrete sizes at the (asymptotic) MNFS optimum for an N of `bits` bits, natural logs, o(1) terms dropped
    and continuous d: ILLUSTRATIVE, not tuned parameters. lam = x^{1/3}(log x)^{2/3}."""
    x = bits * log(2)
    lam = x ** (1 / 3) * log(x) ** (2 / 3)
    d = max(1.0, delta * (x / log(x)) ** (1 / 3))
    logB, logBa, logE = beta * lam, beta_a * lam, e * lam
    log_rat = x / (d + 1) + logE        # rational norm ~ m E
    log_alg = x / (d + 1) + d * logE    # algebraic norm ~ m E^d
    return {"bits": bits, "x": x, "lam": lam, "d": d, "logB": logB, "logE": logE, "u_rat": log_rat / logB,
            "v_box": 2 * logE / logB, "logBa": logBa, "u_alg": log_alg / logBa, "v_box_alg": 2 * logE / logBa}


def lq_gain_table(bits_list=(512, 768, 1024, 2048, 4096, 16384, 65536, 2 ** 20)) -> List[Dict[str, float]]:
    """Proposition LQ, numerically: the LQ upper bound g* (with multiplicity) on each side separately at the MNFS
    point, and sigma_eff = log g*/lam, which tends to 0 like (log x)^{1/3}/x^{1/3}
    (log g* ~ v_box log(u log u) ~ log x, lam = x^{1/3} (log x)^{2/3})."""
    rows = []
    for n in bits_list:
        p = mnfs_point(n)
        b = lq_best_at(p["u_rat"], p["v_box"], p["logB"])
        ba = lq_best_at(p["u_alg"], p["v_box_alg"], p["logBa"])
        rows.append({**p, "v": b["v"], "gain": b["gain"], "gain_alg": ba["gain"],
                     "sigma_eff": log(b["gain"]) / p["lam"], "sigma_eff_alg": log(ba["gain"]) / p["lam"],
                     "frac_rat_dickman_penalty": log(b["gain"]) / -log_rho(p["u_rat"])})
    return rows


def wagner_max_k(log_m: float, log_E: float, kmax: int = 2 ** 20) -> int:
    """k-tree (Wagner 2002) for 'product of k listed smooth numbers lands in the window |s mod m| <= E with
    0 < s <= E m'. The k-tree needs bottom lists of size S >= (m/E)^{1/(1 + lg k)} (window of m/E residues),
    while k factors with product <= E m allow at most S <= (E m)^{1/k} per list. Feasible iff
        (1 + lg k)/k >= log(m/E)/log(E m).
    Returns the largest feasible power of two k. As log E/log m -> 0 (L[1/3] vs L[2/3]) only k <= 2 survives,
    and k = 2 (plain meet-in-the-middle) needs lists of size (m/E)^{1/2} = L[2/3] >> E^2."""
    need = (log_m - log_E) / (log_m + log_E)
    best, k = 1, 2
    while k <= kmax:
        if (1 + log2(k)) / k >= need:
            best = k
        k *= 2
    return best


def wagner_table(bits_list=(512, 1024, 4096, 2 ** 16, 2 ** 20)) -> List[Dict[str, float]]:
    rows = []
    for n in bits_list:
        p = mnfs_point(n)
        log_m = p["x"] / (p["d"] + 1)
        rows.append({"bits": n, "log_m": log_m, "log_E": p["logE"], "max_k": wagner_max_k(log_m, p["logE"]),
                     "log_list_k2": 0.5 * (log_m - p["logE"]), "log_box": 2 * p["logE"]})
    return rows


if __name__ == "__main__":  # pragma: no cover
    for r in lq_gain_table():
        print({k: round(v, 4) if isinstance(v, float) else v for k, v in r.items()})
    for r in wagner_table():
        print({k: round(v, 3) if isinstance(v, float) else v for k, v in r.items()})
