"""Prop DC (decoupling) bookkeeping: the bound of notes/pc3_results.md §2, evaluated on a finite joint law.

For a label with parent count Z >= 0 and child good-point count N >= 0 (joint law p on finitely many atoms):
    E[Z 1{N=0}] <= E[Z] * (K*phi + (sig^2 + sig)/(K - 1)),   K > 1,
phi = Var N / (E N)^2, sig^2 = Var Z / (E Z)^2.  Optimising K: E[Z] * (phi + 2 sqrt(phi (sig^2 + sig))).
"""
from __future__ import annotations

import numpy as np


def _moments(p, Z, N):
    p = np.asarray(p, float)
    Z = np.asarray(Z, float)
    N = np.asarray(N, float)
    ez = (p * Z).sum()
    en = (p * N).sum()
    phi = (p * (N - en) ** 2).sum() / en ** 2
    sig = np.sqrt((p * (Z - ez) ** 2).sum()) / ez
    return ez, phi, sig


def dc_bound(p, Z, N, K: float) -> float:
    ez, phi, sig = _moments(p, Z, N)
    return ez * (K * phi + (sig * sig + sig) / (K - 1))


def dc_bound_opt(p, Z, N) -> float:
    ez, phi, sig = _moments(p, Z, N)
    return ez * (phi + 2 * np.sqrt(phi * (sig * sig + sig)))


def required_margin(eta: float, phi: float, sig: float) -> float:
    """Ratio (phi + 2 sqrt(phi (sig^2+sig))) / (eta/2): Prop UB needs <= 1 on incidence-weighted average."""
    return (phi + 2 * np.sqrt(phi * (sig * sig + sig))) / (eta / 2)
