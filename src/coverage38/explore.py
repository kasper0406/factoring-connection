"""Exploratory (NOT pre-registered) rerun of e2 / e2s with the corrected Lehman score
L*(p,q) = min_v max(B^2, |delta_v| sqrt(B/P))  (height-dependent stationary windows;
Dirichlet gives L* <= 2 P^{2/5}).  Added after the Opus critic found that the frozen
uniform-window score L can exceed C P^{2/5}.  HARNESS code.
Run: PYTHONPATH=. .venv/bin/python -m src.coverage38.explore"""
from __future__ import annotations

import json
import time
from math import isqrt

import numpy as np

from src.dchworst import coverage as cv
from src.coverage38.experiments import sweep_snapshots, surv_curve, h_pred, _write


def lstar_vec(P: np.ndarray, Q: np.ndarray, Bmax: int) -> np.ndarray:
    Pf = np.sqrt(P.astype(np.float64) * Q.astype(np.float64))
    best = np.full(len(P), np.inf)
    for B in range(1, Bmax + 1):
        a0 = (B * P) // Q
        for A in (a0, a0 + 1):
            d = np.abs((B * P - A * Q).astype(np.float64))
            val = np.maximum(float(B * B), d * np.sqrt(B / Pf))
            val = np.where(A >= 1, val, np.inf)
            np.minimum(best, val, out=best)
    return best


def run():
    from src.dchworst.experiments import sample_pairs
    from src.dcltheorem.params import pow2_frac
    # exhaustive n = 32
    Y = 2 ** 13
    band = cv.band(Y)
    P, Q = cv.all_pairs(32)
    Xd = (2.0 ** 16) ** 0.4
    sc = lstar_vec(P, Q, int(2 * Xd ** 0.5) + 4)
    snaps = sweep_snapshots(P, Q, band, 110, True, {16, 64, 128, len(band)})
    for G in (2, 4):
        easy = sc <= Xd / G
        for J, (h, surv) in sorted(snaps.items()):
            rec = {"exp": "x_e2star", "n": 32, "J": J, "G": G, "frac_easy": float(easy.mean()),
                   "max_Lstar_over_P25": float(sc.max() / Xd),
                   "hstar_all": int(h.max()), "hstar_hard": int(h[~easy].max()),
                   "pred_hard": h_pred(int((~easy).sum()), surv)}
            _write("explore", rec)
            print(json.dumps(rec), flush=True)
    del P, Q, snaps
    for n in (60, 72):
        P, Q = sample_pairs(n, 200000, seed0=3800 + n, procs=8)
        Y = pow2_frac(n, 1, 6)
        band = cv.band(Y)
        Xd = (2.0 ** (n / 2)) ** 0.4
        sc = lstar_vec(P, Q, int(2 * Xd ** 0.5) + 4)
        h, surv = sweep_snapshots(P, Q, band, 40, True, {len(band)})[len(band)]
        tail = surv_curve(h, 40)
        Href = min(H for H in range(41) if tail[H + 1] <= 0.05 * len(P))
        unc = h > Href
        for G in (2, 4):
            easy = sc <= Xd / G
            fe, fh = float(unc[easy].mean()), float(unc[~easy].mean())
            ne, nh = int(easy.sum()), int((~easy).sum())
            se = (fe * (1 - fe) / ne) ** 0.5 / fh
            rec = {"exp": "x_e2sstar", "n": n, "G": G, "frac_easy": ne / len(P),
                   "max_Lstar_over_P25": float(sc.max() / Xd), "Href": int(Href),
                   "unc_easy": fe, "unc_hard": fh, "R": fe / fh,
                   "R_ci95": [fe / fh - 1.96 * se, fe / fh + 1.96 * se]}
            _write("explore", rec)
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    run()
