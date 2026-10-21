"""L[1/3] constants for Prop SI, its prime-field DL transfer (Prop DLI), and the model collapses.

Floating point is used for constants only (closed forms are given in the docstrings and checked in tests).
"""
from __future__ import annotations

import math

MNFS = 1.9019


def prop_si_cost(kappa: float, nu: float, tau: float = 0.0) -> float:
    """L[1/3] constant of the Prop SI program: max(kappa + tau + nu/(3 kappa), 2 kappa) (needs nu > 2/(3 kappa))."""
    return max(kappa + tau + nu / (3 * kappa), 2 * kappa)


def prop_si_optimum(tau: float = 0.0, grid: int = 200000):
    """inf over kappa of prop_si_cost(kappa, 2/(3 kappa), tau). tau = 0: kappa^3 = 2/9, c^3 = 16/9 (closed form)."""
    if tau == 0.0:
        k = (2 / 9) ** (1 / 3)
        return k, 2 * k
    best = (None, float("inf"))
    for i in range(1, grid):
        k = 3.0 * i / grid
        c = prop_si_cost(k, 2 / (3 * k), tau)
        if c < best[1]:
            best = (k, c)
    return best


def tau_threshold_mnfs() -> float:
    """Largest per-call constructor cost tau with Prop SI below MNFS: 1.9019 - (3/2)(4/9)^(1/3) (suffmap m4(d))."""
    return MNFS - 1.5 * (4 / 9) ** (1 / 3)


def dl_booting_ecm() -> float:
    """Classical ECM booting for individual logs in F_p: min_c 1/(3c) + sqrt(4c/3) = 3^(1/3) at c = 3^(-1/3)."""
    c = 3 ** (-1 / 3)
    return 1 / (3 * c) + math.sqrt(4 * c / 3)


def dl_descent_cost(c: float, kappa: float, nu: float | None = None) -> float:
    """Descent of one medium prime q = L[2/3, c] with the interval oracle on modulus M = floor(jp/q):
    v = q u - kjp has |v| <= q (h + k) = L[2/3, c + nu]; y-smooth with prob L[1/3, -(c+nu)/(3 kappa)]."""
    if nu is None:
        nu = 2 / (3 * kappa)
    return (c + nu) / (3 * kappa)


def dl_descent_optimum(kappa: float, grid: int = 200000):
    """Ignoring the ECM cost of booting (smoothness test at L^{o(1)}): min_c max(1/(3c), descent(c))."""
    best = (None, float("inf"))
    for i in range(1, grid):
        c = 3.0 * i / grid
        v = max(1 / (3 * c), dl_descent_cost(c, kappa))
        if v < best[1]:
            best = (c, v)
    return best


def dl_individual_log_with_ecm(kappa: float, grid: int = 200000):
    """With ECM booting cost: min_c max(1/(3c) + sqrt(4c/3), descent(c)). Returns (c, cost)."""
    best = (None, float("inf"))
    for i in range(1, grid):
        c = 3.0 * i / grid
        v = max(1 / (3 * c) + math.sqrt(4 * c / 3), dl_descent_cost(c, kappa))
        if v < best[1]:
            best = (c, v)
    return best


def johnson_nu_needed_if_half_power(kappa: float) -> float:
    """Superseded crude bound (E <= V^{1/2} rho): would need nu >= 4/(3 kappa). Kept to document that the
    proved Theorem J gives the much stronger E <= L[2/3, -2/(3 kappa) + o(1)] for every nu."""
    return 4 / (3 * kappa)


def collapse_table(kappa: float | None = None):
    """Per-output constructor cost in each model at Prop SI scale (exponent class, constant)."""
    if kappa is None:
        kappa = (2 / 9) ** (1 / 3)
    g = 2 / (3 * kappa)
    return {
        "S (sampling/sieving)": ("L[2/3]", g),
        "J (Johnson certificates, oblivious modulus)": ("L[2/3]", g),
        "M2 (oblivious two-list merge)": ("exp(n/2)", 0.5),
        "BKZ/GSA Schnorr log-lattice": ("exp(L[1/3])", None),
        "whole relation program in model S": ("L[1/2]", math.sqrt(2)),
    }


if __name__ == "__main__":
    print("Prop SI optimum", prop_si_optimum())
    print("tau threshold", tau_threshold_mnfs())
    k = (2 / 9) ** (1 / 3)
    print("DL booting ECM", dl_booting_ecm())
    print("DL descent (no ECM cost)", dl_descent_optimum(k))
    print("DL individual log with ECM booting", dl_individual_log_with_ecm(k))
    print(collapse_table())
