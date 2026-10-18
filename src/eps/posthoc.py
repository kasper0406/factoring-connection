"""POST-HOC (not pre-registered) mechanism analysis of the n = 32 CF hits.  Written after the battery ran.

Hypothesis H-SQ: the CF-group hits are carried by steps whose Q_i has a nontrivial square factor.  For odd i,
A_{i-1}^2 = -Q_i (mod N); if Q_i = k^2 m then (m|p) = eps * (stuff) -- in particular Q_i = k^2 at odd i is impossible
when eps = -1 (the PC2 mechanism).  Test: recompute the eight CF count statistics over T = 1024 steps,
(a) all steps (must reproduce the battery values), (b) only steps with Q_i squarefree,
(c) only steps with Q_i NOT squarefree.  The step filter uses only Q_i (public), so this stays operational.
Squarefreeness: valid for Q_i < 2^18 (n <= 34): every square factor l^2 | Q_i has l <= 2^9.
"""
from __future__ import annotations

import argparse
import gzip
import json
from math import gcd
from multiprocessing import Pool

import gmpy2
import numpy as np

_P = 1
for _l in range(2, 513):
    if all(_l % d for d in range(2, int(_l ** 0.5) + 1)):
        _P *= _l

KEYS = ("n_a1", "n_a2", "sumbits", "maxbits", "n_Qeven", "n_sq_even", "n_Q1mod4_odd", "n_2sq_odd")


def squarefree_small(Q: int) -> bool:
    assert Q < (1 << 18)
    # g = product of the distinct primes <= 512 dividing Q.  Since Q < 2^18, at most one prime > 512 divides Q,
    # to the first power; so Q is squarefree iff no small prime survives in Q // g.
    g = gcd(Q, _P)
    return gcd(Q // g, g) == 1


def cf_variants(N: int, T: int = 1024):
    """OPERATIONAL (N only).  Returns {variant: {key: value}} for variants all / sqfree / nonsqfree."""
    a0 = int(gmpy2.isqrt(N))
    P, Q, a = 0, 1, a0
    out = {v: {k: 0 for k in KEYS} for v in ("all", "sqfree", "nonsqfree")}
    for i in range(1, T + 1):
        P = a * Q - P
        Q = (N - P * P) // Q
        a = (a0 + P) // Q
        odd = i & 1
        sf = squarefree_small(Q)
        is_sq = bool(gmpy2.is_square(Q))
        is_2sq = (Q & 1) == 0 and bool(gmpy2.is_square(Q >> 1))
        for v in ("all", "sqfree" if sf else "nonsqfree"):
            d = out[v]
            d["n_a1"] += a == 1
            d["n_a2"] += a == 2
            d["sumbits"] += a.bit_length()
            d["maxbits"] = max(d["maxbits"], a.bit_length())
            d["n_Qeven"] += (Q & 1) == 0
            d["n_sq_even"] += (not odd) and is_sq
            d["n_Q1mod4_odd"] += bool(odd) and (Q & 3) == 1
            d["n_2sq_odd"] += bool(odd) and is_2sq
    return out


def _work(N):
    return cf_variants(N)


def main(argv=None):
    from src.eps.analysis import load, strat_z
    from src.eps.operational import MENU_NAMES
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--out", default="experiments/analysis/eps_posthoc_n32.json")
    a = ap.parse_args(argv)
    d = load(a.path)
    with Pool(30) as pool:
        V = pool.map(_work, d["N"], chunksize=500)
    res = {"path": a.path, "n_inst": len(V)}
    # consistency with the battery values
    mism = 0
    for k in KEYS:
        j = MENU_NAMES.index(f"cf1024_{k}")
        mism += int((np.array([v["all"][k] for v in V]) != d["X"][:, j]).sum())
    res["reproduce_mismatch"] = mism
    for tgt in ("eps", "qr2"):
        mask = np.ones(len(V), bool) if tgt == "eps" else d["qr2"] >= 0
        y = ((d["eps"] == 1).astype(int) if tgt == "eps" else d["qr2"])[mask]
        strata = d["Nmod"][mask] % 256
        rows = {}
        for k in KEYS:
            rows[k] = {var: strat_z(np.array([v[var][k] for v in V], float)[mask], y, strata)["z"]
                       for var in ("all", "sqfree", "nonsqfree")}
        res[tgt] = rows
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()


def cluster_robust_z(x, y, strata, cp, cq):
    """POST-HOC: z = U / sqrt(V_cr), U as in analysis.strat_z, V_cr two-way cluster-robust (clusters = p and q):
    V_cr = V_p + V_q - V_{pq}, V_c = sum_c (sum_{i in c} e_i)^2 with e_i = (x_i - xbar_h)(y_i - ybar_h)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    e = np.zeros(len(x))
    for h in np.unique(strata):
        m = strata == h
        e[m] = (x[m] - x[m].mean()) * (y[m] - y[m].mean())
    U = e.sum()

    def vc(keys):
        acc = {}
        for k, v in zip(keys, e):
            acc[k] = acc.get(k, 0.0) + v
        return sum(v * v for v in acc.values())
    V = vc(cp) + vc(cq) - vc(list(zip(cp, cq)))
    return U / np.sqrt(V) if V > 0 else 0.0


def cluster_main(path, stats, out):
    from src.eps.analysis import load
    from src.eps.operational import MENU, MENU_NAMES
    d = load(path)
    res = {}
    for tgt in ("eps", "qr2"):
        mask = np.ones(len(d["eps"]), bool) if tgt == "eps" else d["qr2"] >= 0
        y = ((d["eps"] == 1).astype(int) if tgt == "eps" else d["qr2"])[mask]
        strata = d["Nmod"][mask] % 256
        cp = [p for p, m in zip(d["p"], mask) if m]
        cq = [q for q, m in zip(d["q"], mask) if m]
        res[tgt] = {}
        names = stats if stats else [m for m, g, r in MENU if r in ("menu", "nc")]
        for s in names:
            x = d["X"][mask, MENU_NAMES.index(s)]
            if np.all(x == x[0]):
                continue
            res[tgt][s] = cluster_robust_z(x, y, strata, cp, cq)
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    return res
