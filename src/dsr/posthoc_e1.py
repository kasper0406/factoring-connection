"""POST HOC (labelled; not pre-registered): where the Fourier mass of smooth sets sits, by character order.

For prime moduli p, characters of (Z/p)^x are indexed by k in Z/(p-1), with order (p-1)/gcd(k, p-1).
Reports, for S_can and a random control of equal size: the share of the L^1 and L^2 (non-principal) mass carried by
characters of order <= T, and A'(S) restricted to characters of order > T.
  PYTHONPATH=. .venv/bin/python -m src.dsr.posthoc_e1
"""
from __future__ import annotations

import json
import math

import numpy as np

from src.common.instances import prime_instance
from src.dsr import core
from src.dsr.harness import A_CAN, A_BAL, fourier_stats, group_array  # noqa: F401


def orders(p: int) -> np.ndarray:
    k = np.arange(p - 1, dtype=np.int64)
    return (p - 1) // np.gcd(k, p - 1)


def run() -> None:
    for bits, seed in ((20, 0), (22, 0)):
        p = prime_instance(bits, seed).N
        G = core.units_mask(p)
        ordv = orders(p)
        for name, a in (("S_can", A_CAN), ("S_bal", A_BAL)):
            B, P = core.factor_base(p, a)
            S = (core.s_can(p, P) if name == "S_can" else core.s_bal(p, P)) & G
            size = int(S.sum())
            rng = np.random.default_rng(7)
            R = np.zeros(p, dtype=bool)
            R[rng.choice(np.nonzero(G)[0], size=size, replace=False)] = True
            out = {"bits": bits, "p": p, "set": name, "B": B, "size": size}
            for tag, M in (("S", S), ("rand", R)):
                F = np.abs(np.fft.fft(group_array(M, [p])))
                F[0] = 0.0
                l1, l2 = F.sum(), (F ** 2).sum()
                for T in (2, 10, 100, 1000, 10**4):
                    sel = ordv <= T
                    out[f"{tag}_L1share_ord<={T}"] = round(float(F[sel].sum() / l1), 4)
                    out[f"{tag}_L2share_ord<={T}"] = round(float((F[sel] ** 2).sum() / l2), 4)
                    out[f"{tag}_A1_ord>{T}"] = round(float(F[~sel].sum() / (p - 1)), 2)
                top = np.argsort(F)[-5:][::-1]
                if tag == "S":
                    ind = np.zeros(p, dtype=np.int64)
                    from src.dsr.harness import _dlog_table
                    ind = _dlog_table(p)
                    Pi = np.array([ind[q] for q in P], dtype=np.int64)
                    # mean Re chi_k(q) over the factor base at the top character (chi_k(q) = e(-k ind(q)/(p-1)))
                    k = int(top[0])
                    out["S_top_meanRe_chi_on_FB"] = round(float(np.mean(np.cos(2 * np.pi * k * Pi / (p - 1)))), 4)
                    # empirical spectral parameter of the factor base: max_{k != 0} |sum_q chi_k(q)| / |P|
                    fb = np.zeros(p - 1)
                    fb[Pi] = 1.0
                    Ffb = np.abs(np.fft.fft(fb))
                    Ffb[0] = 0.0
                    out["lambda_emp_FB"] = round(float(Ffb.max() / len(P)), 4)
                    out["nP"] = len(P)
                    out["GRH_bound_log2(pB)/sqrtB"] = round(math.log(p * B) ** 2 / math.sqrt(B), 2)
                out[f"{tag}_top5_orders"] = [int(ordv[i]) for i in top]
                out[f"{tag}_top5_bias"] = [round(float(F[i] / size), 4) for i in top]
            print(json.dumps(out))


if __name__ == "__main__":
    run()
