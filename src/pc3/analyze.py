"""Analysis of pc3 E1/E2 (notes/pc3_prereg.md): pooled statistics, label-cluster bootstrap, verdicts.

Usage: python -m src.pc3.analyze E1 experiments/raw/pc3_e1.jsonl
       python -m src.pc3.analyze E2 experiments/raw/pc3_e2.jsonl
"""
from __future__ import annotations

import json
import sys

import numpy as np


def boot(vals, weights=None, n=2000, seed=7):
    vals = np.asarray(vals, float)
    w = np.ones_like(vals) if weights is None else np.asarray(weights, float)
    rng = np.random.default_rng(seed)
    k = len(vals)
    out = []
    for _ in range(n):
        idx = rng.integers(0, k, k)
        out.append((vals[idx] * w[idx]).sum() / w[idx].sum())
    return float((vals * w).sum() / w.sum()), float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975))


def e1(path):
    R = [json.loads(l) for l in open(path)]
    res = {}
    m, lo, hi = boot([r["mean_resid"] for r in R])
    res["P1a_mean_resid"] = (m, lo, hi, -0.01 <= m <= 0.01 and -0.02 <= lo and hi <= 0.02)
    cov = np.array([r["cov_emp_mod"] for r in R])
    var = np.array([r["var_mod"] for r in R])
    rng = np.random.default_rng(8)
    sl = [cov[i].sum() / var[i].sum() for i in (rng.integers(0, len(R), len(R)) for _ in range(2000))]
    slope = cov.sum() / var.sum()
    res["P1b_slope"] = (float(slope), float(np.quantile(sl, .025)), float(np.quantile(sl, .975)), 0.8 <= slope <= 1.3)
    exc = np.mean([r["var_emp"] - r["noise_var"] - r["var_mod"] for r in R])
    vm = np.mean(var)
    res["P1c_excess_dispersion"] = (float(exc), float(vm), bool(exc <= 0.5 * vm))
    m, lo, hi = boot([r["mean_Cemp"] for r in R])
    res["P1d_box_mean_Cemp"] = (m, lo, hi, -0.01 <= m <= 0.01)
    res["descr"] = dict(n_labels=len(R), mean_m=float(np.mean([r["m"] for r in R])),
                        mean_Cloc=float(np.mean([r["mean_Cloc"] for r in R])),
                        mean_Carch=float(np.mean([r["mean_Carch"] for r in R])),
                        mean_abs_resid=float(np.mean([r["mean_abs_resid"] for r in R])),
                        sd_Cemp_pairs=float(np.sqrt(np.mean([r["var_emp"] for r in R]))),
                        sd_Cmod_pairs=float(np.sqrt(vm)),
                        bin_mod=list(np.mean([r["bin_mod"] for r in R], 0)),
                        bin_emp=list(np.mean([r["bin_emp"] for r in R], 0)))
    return res


def e2(path):
    R = [json.loads(l) for l in open(path)]
    res = {}
    Ks = sorted(set(r["K"] for r in R))
    byK = {K: [r for r in R if r["K"] == K] for K in Ks}
    for K in Ks:
        rs = byK[K]
        d = [r["Phi_ex"] - r["phi_ex"] for r in rs]
        m, lo, hi = boot(d)
        res["P2a_K%d" % K] = (m, lo, hi, abs(m) <= 0.01)
        m, lo, hi = boot([r["sig2_ex"] for r in rs])
        res["P2c_K%d" % K] = (m, lo, hi, -0.02 <= m <= 0.05)
        m, lo, hi = boot([r["shift"] for r in rs])
        res["P2d_K%d" % K] = (m, lo, hi, abs(m) <= 0.03)
        m, lo, hi = boot([r["dc_direct"] - r["phi"] for r in rs])
        res["P2e_K%d" % K] = (m, lo, hi, abs(m) <= 0.01)
        res["descr_K%d" % K] = dict(phi_ex=float(np.mean([r["phi_ex"] for r in rs])),
                                    Phi_ex=float(np.mean([r["Phi_ex"] for r in rs])),
                                    z=float(np.mean([r["z"] for r in rs])), sig2=float(np.mean([r["sig2"] for r in rs])),
                                    m=float(np.mean([r["m"] for r in rs])),
                                    n_par_used=float(np.mean([r["n_par_used"] for r in rs])))
    a = {r["q"] * 10 ** 6 + r["r"]: r["Phi_ex"] for r in byK[Ks[-1]]}
    b = {r["q"] * 10 ** 6 + r["r"]: r["Phi_ex"] for r in byK[Ks[0]]}
    d = [a[k] - b[k] for k in a if k in b]
    m, lo, hi = boot(d)
    res["P2b_K%d_minus_K%d" % (Ks[-1], Ks[0])] = (m, lo, hi, abs(m) <= 0.01)
    return res


if __name__ == "__main__":
    which, path = sys.argv[1], sys.argv[2]
    out = e1(path) if which == "E1" else e2(path)
    out = json.loads(json.dumps(out, default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o)))
    print(json.dumps(out, indent=1))
    with open(path.replace(".jsonl", "_summary.json"), "w") as f:
        json.dump(out, f, indent=1)
