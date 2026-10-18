"""sieve2d EXPLORATORY (not pre-registered) diagnostics, run after the frozen runs.

  PYTHONPATH=. .venv/bin/python -m src.sieve2d.explore

X1: integer-pair 4th moments with the points x == (0,0) mod some band prime (the set D of
    notes §2) removed, vs the integer model conditioned on the same event.  Purpose:
    diagnose the frozen P5 failure (heavy-tailed atom M_l = |V| with probability 1/l^2).
X2: 6th central moment of M over the same 10^6 prime pairs (block 0) vs the
    independent-residue model.  Purpose: is the proof's j = 6 obstruction (notes §4)
    visible in the data, or only a proof barrier?
Same seeds/samplers as the pre-registered E2 runs.  Raw: experiments/raw/sieve2d_explore.jsonl
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List

import numpy as np

from src.dcltheorem.params import weight_set

from .coverage import class_table, degenerate_mask, inverse_table, lam_nominal
from .experiments import CONFIG, RAW, band, thm_Y
from .sampling import sample_integer_pairs, sample_prime_pairs


def _cum_from_raw(m):
    """cumulants k1..k6 from raw moments m[1..6]."""
    m1, m2, m3, m4, m5, m6 = m[1:7]
    k1 = m1
    k2 = m2 - m1 ** 2
    k3 = m3 - 3 * m2 * m1 + 2 * m1 ** 3
    k4 = m4 - 4 * m3 * m1 - 3 * m2 ** 2 + 12 * m2 * m1 ** 2 - 6 * m1 ** 4
    k5 = (m5 - 5 * m4 * m1 - 10 * m3 * m2 + 20 * m3 * m1 ** 2 + 30 * m2 ** 2 * m1
          - 60 * m2 * m1 ** 3 + 24 * m1 ** 5)
    k6 = (m6 - 6 * m5 * m1 - 15 * m4 * m2 + 30 * m4 * m1 ** 2 - 10 * m3 ** 2
          + 120 * m3 * m2 * m1 - 120 * m3 * m1 ** 3 + 30 * m2 ** 3 - 270 * m2 ** 2 * m1 ** 2
          + 360 * m2 * m1 ** 4 - 120 * m1 ** 6)
    return [0, k1, k2, k3, k4, k5, k6]


def model_central(vals_wts_per_l) -> Dict[str, float]:
    K = [0.0] * 7
    for vals, wts in vals_wts_per_l:
        m = [float((wts * vals ** i).sum()) for i in range(7)]
        k = _cum_from_raw(m)
        for i in range(7):
            K[i] += k[i]
    mu4 = K[4] + 3 * K[2] ** 2
    mu6 = K[6] + 15 * K[4] * K[2] + 10 * K[3] ** 2 + 15 * K[2] ** 3
    return {"mean": K[1], "var": K[2], "mu4": mu4, "mu6": mu6}


def _job(args):
    n, kind, Y, Hs, k = args
    cnt = CONFIG["samp_pairs"]
    if kind == "prime":
        p, q = sample_prime_pairs(n, cnt, CONFIG["samp_seed"])
    else:
        p, q = sample_integer_pairs(n, cnt, CONFIG["samp_seed"])
    bl = band(Y, k)
    Vs = [weight_set(H) for H in Hs]
    Ms = [np.zeros(p.shape, dtype=np.int64) for _ in Vs]
    inD = np.zeros(p.shape, dtype=bool)
    models: List[List] = [[] for _ in Vs]
    for l in bl:
        inv = inverse_table(l)
        pm, qm = p % l, q % l
        qz = qm == 0
        both = qz & (pm == 0)
        inD |= both
        r = np.where(qz, l, (pm * inv[qm]) % l)
        for i, V in enumerate(Vs):
            C = class_table(l, V)
            Ms[i] += np.where(both, len(V), C[r])
            if kind == "prime":
                vals = C[1:l].astype(float)
                wts = np.full(l - 1, 1.0 / (l - 1))
            else:   # integer model conditioned on x != (0,0) mod l
                vals = C.astype(float)
                wts = np.full(l + 1, (l - 1.0) / (l * l - 1))
            models[i].append((vals, wts))
    out = []
    for H, V, M, mod in zip(Hs, Vs, Ms, models):
        lam = lam_nominal(bl, len(V))
        mc = model_central(mod)
        # off D and, for integer pairs, off the lines L (as in the frozen E2 records)
        keep = ~inD
        if kind == "int":
            keep &= ~degenerate_mask(p, q, V)
        Mk = M[keep].astype(float)
        mu = Mk.mean()
        rec = {"n": n, "kind": kind, "block": k, "H": H, "Y": Y, "lam": lam,
               "pairs": int(M.size), "nD": int(inD.sum()), "nkept": int(keep.sum()),
               "m4_offD": float(((Mk - mu) ** 4).mean()), "model_mu4_offD": mc["mu4"],
               "m6_offD": float(((Mk - mu) ** 6).mean()), "model_mu6_offD": mc["mu6"],
               "var_offD": float(Mk.var()), "model_var_offD": mc["var"]}
        if kind == "int":
            nd = ~degenerate_mask(p, q, V)          # = the frozen E2 population
            Mf = M[nd].astype(float)
            rec["m4_frozen_pop"] = float(((Mf - Mf.mean()) ** 4).mean())
            rec["dev4_share_from_D"] = float(((Mf[inD[nd]] - Mf.mean()) ** 4).sum()
                                             / ((Mf - Mf.mean()) ** 4).sum())
        out.append(rec)
    return out


def main():
    jobs = []
    for n in CONFIG["samp_n"]:
        byY: Dict[int, List[int]] = {}
        for H in CONFIG["samp_H"]:
            Y = thm_Y(n, H)
            if Y > H * H:
                byY.setdefault(Y, []).append(H)
        for Y, Hs in byY.items():
            for k in CONFIG["samp_blocks"][str(n)]:
                jobs.append((n, "int", Y, Hs, k))
                if k == 0:
                    jobs.append((n, "prime", Y, Hs, k))
    path = os.path.join(RAW, "sieve2d_explore.jsonl")
    with ProcessPoolExecutor(max_workers=8) as ex, open(path, "w") as fh:
        for res in ex.map(_job, jobs):
            for r in res:
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(res[0]["n"], res[0]["kind"], res[0]["block"], [x["H"] for x in res], flush=True)


if __name__ == "__main__":
    main()
