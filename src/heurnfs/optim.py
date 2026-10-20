"""Seeded Nelder-Mead with penalty constraints (agent heurnfs).

ANALYSIS ONLY: floating point is used to *minimise L[1/3]-constant programs*, i.e. as a fit, never
inside an operational factoring routine (AGENT_RULES rule 1). Every program whose optimum has a closed
form is cross-checked against the exact sympy value in tests/test_heurnfs.py.
"""
from __future__ import annotations

from typing import Callable, Sequence, Tuple

import numpy as np


def _nm(f: Callable, x0: np.ndarray, iters: int, scale: float) -> Tuple[np.ndarray, float]:
    n = len(x0)
    S = [np.array(x0, float)]
    for i in range(n):
        y = np.array(x0, float)
        y[i] += scale
        S.append(y)
    F = [f(s) for s in S]
    for _ in range(iters):
        o = np.argsort(F)
        S = [S[i] for i in o]
        F = [F[i] for i in o]
        c = np.mean(S[:-1], 0)
        xr = c + (c - S[-1]); fr = f(xr)
        if fr < F[0]:
            xe = c + 2 * (c - S[-1]); fe = f(xe)
            S[-1], F[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < F[-2]:
            S[-1], F[-1] = xr, fr
        else:
            xc = c + 0.5 * (S[-1] - c); fc = f(xc)
            if fc < F[-1]:
                S[-1], F[-1] = xc, fc
            else:
                S = [S[0] + 0.5 * (s - S[0]) for s in S]
                F = [f(s) for s in S]
    i = int(np.argmin(F))
    return S[i], F[i]


def minimise(obj: Callable[[np.ndarray], float], cons: Callable[[np.ndarray], Sequence[float]],
             lo: Sequence[float], hi: Sequence[float], restarts: int = 40, seed: int = 0,
             pen: float = 1e3) -> Tuple[np.ndarray, float]:
    """min obj(x) s.t. every g in cons(x) >= 0 and x > 0 (exact penalty, several restarts).

    Returns (x*, obj(x*)); the caller checks feasibility with `max_violation`.
    """
    rng = np.random.default_rng(seed)
    lo = np.asarray(lo, float); hi = np.asarray(hi, float)

    def f(x):
        if np.any(x <= 1e-9):
            return 1e9
        return obj(x) + pen * sum(max(0.0, -g) for g in cons(x))

    best = None
    for _ in range(restarts):
        x, v = _nm(f, rng.uniform(lo, hi), 3000, 0.2)
        x, v = _nm(f, x, 3000, 0.01)
        x, v = _nm(f, x, 2000, 0.0005)
        if best is None or v < best[1]:
            best = (x, v)
    x = best[0]
    return x, float(obj(x))


def max_violation(cons, x) -> float:
    return max([0.0] + [-g for g in cons(x)])
