"""Linear-statistics barrier (notes/coverage38_results.md §3).  DIAGNOSTIC / HARNESS.

Reciprocal model.  For a pair (p, q) and band primes l, x_l = p * l^{-1} mod q.
l | Bp - Aq with |A| <~ H  <=>  B x_l mod q is ~ a multiple of l of size <~ H q / Y,
so coverage is "some x_l falls in the Farey-arc set
    S = {x in Z/q : exists 1 <= B <= H with |B x mods q| <= k0}."
K(h) = sum_l e(h x_l / q).  For fixed q, K_p(h) = K_1(h p mod q).

Exact parts: the point sets, the target set, hit counts (integers).
Floating point (numpy FFT) is used ONLY for the magnitudes |K(h)|, |S^(h)| -- the
analytic statistics whose size is the object of study; nothing is a factoring output.
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np


def farey_arc_set(q: int, H: int, k0: int) -> np.ndarray:
    S = np.zeros(q, dtype=bool)
    r = np.arange(-k0, k0 + 1, dtype=np.int64) % q
    for B in range(1, H + 1):
        Binv = pow(B, -1, q)
        S[(r * Binv) % q] = True
    return S


def points(p: int, q: int, band: Sequence[int]) -> np.ndarray:
    return np.array([(p * pow(l, -1, q)) % q for l in band], dtype=np.int64)


def fourier_abs(ind: np.ndarray) -> np.ndarray:
    """|sum_x ind[x] e(h x / q)| for all h (float)."""
    return np.abs(np.fft.fft(ind.astype(np.float64)))


def stats_for_p(absK1: np.ndarray, absS: np.ndarray, p: int, q: int, Q0: int) -> Dict:
    """Uses |K_p(h)| = |K_1(h p)|.  E_abs = (1/q) sum_{h != 0} |S^(h)| |K_p(h)|
    (the absolute-value Erdos-Turan error bound); Mlow = max_{1<=h<=Q0} |K_p(h)|."""
    h = np.arange(q, dtype=np.int64)
    idx = (h * p) % q
    Kp = absK1[idx]
    E_abs = float((absS[1:] * Kp[1:]).sum() / q)
    Mlow = float(Kp[1:Q0 + 1].max())
    Mall = float(Kp[1:].max())
    return {"E_abs": E_abs, "Mlow": Mlow, "Mall": Mall}


def random_complement(q: int, S: np.ndarray, J: int, rng: np.random.Generator) -> np.ndarray:
    comp = np.nonzero(~S)[0]
    return rng.choice(comp, size=J, replace=False)


def hoeffding_level(J: int, Q: int) -> float:
    """Deviation level 2 sqrt(J log(8Q)) of Prop. L (simultaneously for Q test sums
    with coefficients of modulus <= 1, failure probability <= 1/2)."""
    return 2.0 * np.sqrt(J * np.log(8.0 * Q))
