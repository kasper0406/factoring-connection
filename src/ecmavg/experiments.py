"""Pre-registered runs (notes/ecmavg_prereg.md).  HARNESS: factor-aware by design.

  .venv/bin/python -m src.ecmavg.experiments labels      # writes experiments/raw/ecmavg_labels.jsonl
  .venv/bin/python -m src.ecmavg.experiments analyze     # writes experiments/analysis/ecmavg_summary.json

Samples (frozen):
  S1  128 primes per b in {12, 14, 16}: random.Random(20261018 + b).sample(primes_b, 128)
  S2  all primes of bit length b in {12, 14, 16}
  S3  2048 primes per b in {18, 20, 22}: random.Random(20261018 + b).sample(primes_b, 2048)
Families: 'suyama' (sigma = 6..69) and 'mont' (A = 3..66, x0 = 2); J_max = 64 curves.
Labels are computed at beta = 1/2 and truncated to beta = 1/3, 1/4.
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from fractions import Fraction as F

import numpy as np
from sympy import primerange

from .curves import bound_B, stage_steps
from .labels import INF, exposure_labels, truncate

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw", "ecmavg_labels.jsonl")
SUMMARY = os.path.join(ROOT, "experiments", "analysis", "ecmavg_summary.json")
FAMILIES = ("suyama", "mont")
JMAX = 64
JS = (8, 16, 32, 64)
BETAS = (F(1, 2), F(1, 3), F(1, 4))
SEED = 20261018


def xbound(b: int) -> int:
    return 2 ** b + 2 ** ((b + 1) // 2 + 1) + 1


def samples():
    out = []
    for b in (12, 14, 16):
        pr = list(primerange(2 ** (b - 1), 2 ** b))
        out.append(("S1", b, sorted(random.Random(SEED + b).sample(pr, 128))))
        out.append(("S2", b, pr))
    for b in (18, 20, 22):
        pr = list(primerange(2 ** (b - 1), 2 ** b))
        out.append(("S3", b, sorted(random.Random(SEED + b).sample(pr, 2048))))
    return out


def run_labels():
    os.makedirs(os.path.dirname(RAW), exist_ok=True)
    done = {}
    with open(RAW, "w") as fh:
        for sname, b, primes in samples():
            B = bound_B(b, F(1, 2))
            for fam in FAMILIES:
                key = (fam, b)
                t0 = time.time()
                # S1 is a subset of S2 at the same b: reuse the S2 labels when available
                k = exposure_labels(fam, JMAX, primes, B, xbound(b))
                for i, p in enumerate(primes):
                    fh.write(json.dumps({"sample": sname, "family": fam, "b": b, "B": B,
                                         "p": p, "k": [int(x) for x in k[:, i]]}) + "\n")
                print(sname, fam, b, len(primes), f"{time.time() - t0:.1f}s", flush=True)


def load_labels():
    data = {}
    with open(RAW) as fh:
        for line in fh:
            r = json.loads(line)
            data.setdefault((r["sample"], r["family"], r["b"]), []).append((r["p"], r["k"]))
    out = {}
    for key, rows in data.items():
        rows.sort()
        out[key] = (np.array([p for p, _ in rows], dtype=np.int64),
                    np.array([k for _, k in rows], dtype=np.int16).T)
    return out


def pair_stats(k: np.ndarray, Js=JS):
    """Moments over all unordered pairs of columns of k (J, M)."""
    J0, M = k.shape
    X = (k != INF)
    res = {J: dict(s1=0.0, s2=0.0, f=0, u1=0.0, u2=0.0, fu=0, dbl0=0, same_nz=0) for J in Js}
    rows = np.array([J - 1 for J in Js])
    npairs = M * (M - 1) // 2
    for i in range(M - 1):
        ks = (k[:, i:i + 1] != k[:, i + 1:])
        xu = (X[:, i:i + 1] != X[:, i + 1:])
        cs = np.cumsum(ks, axis=0, dtype=np.int32)[rows]
        cu = np.cumsum(xu, axis=0, dtype=np.int32)[rows]
        w0 = np.cumsum(X[:, i:i + 1] | X[:, i + 1:], axis=0, dtype=np.int32)[rows]
        for t, J in enumerate(Js):
            zs, zu = cs[t].astype(np.float64), cu[t].astype(np.float64)
            r = res[J]
            r["s1"] += zs.sum(); r["s2"] += (zs * zs).sum(); r["f"] += int((zs == 0).sum())
            r["u1"] += zu.sum(); r["u2"] += (zu * zu).sum(); r["fu"] += int((zu == 0).sum())
            fail = (zs == 0)
            r["dbl0"] += int((fail & (w0[t] == 0)).sum())
            r["same_nz"] += int((fail & (w0[t] > 0)).sum())
    out = {}
    for J in Js:
        r = res[J]
        ms, mu = r["s1"] / npairs, r["u1"] / npairs
        out[J] = {"E_Zs_over_J": ms / J, "Var_Zs_over_J": (r["s2"] / npairs - ms * ms) / J,
                  "fail_s": r["f"] / npairs,
                  "E_Zu_over_J": mu / J, "Var_Zu_over_J": (r["u2"] / npairs - mu * mu) / J,
                  "fail_u": r["fu"] / npairs,
                  "fail_s_both_unexposed": r["dbl0"] / npairs,
                  "fail_s_same_nonzero_signature": r["same_nz"] / npairs}
    return out, npairs


def prime_stats(primes: np.ndarray, k: np.ndarray, Js=JS, L: int = 120):
    """Single-prime exposure count W_p(J) and its common-mode decomposition by p mod L."""
    X = (k != INF).astype(np.float64)
    out = {}
    cls = primes % L
    for J in Js:
        W = X[:J].sum(axis=0)
        var = W.var()
        # between-class variance of E[W | p mod L]
        means = {c: W[cls == c].mean() for c in np.unique(cls)}
        wts = {c: (cls == c).mean() for c in means}
        between = sum(wts[c] * (means[c] - W.mean()) ** 2 for c in means)
        C = np.cov(X[:J], bias=True)
        off = (C.sum() - np.trace(C)) / (J * (J - 1))
        out[J] = {"E_W_over_J": W.mean() / J, "Var_W_over_J": var / J,
                  "P_W0": float((W == 0).mean()), "between_class_share": between / var if var else 0.0,
                  "mean_offdiag_cov": off, "mean_diag_var": float(np.trace(C) / J)}
    sig = {}
    for i in range(k.shape[1]):
        sig.setdefault(tuple(X[:, i].astype(np.int8)), 0)
        sig[tuple(X[:, i].astype(np.int8))] += 1
    sigs = {}
    for i in range(k.shape[1]):
        sigs.setdefault(tuple(k[:, i]), 0)
        sigs[tuple(k[:, i])] += 1
    M = k.shape[1]
    npairs = M * (M - 1) // 2
    ident_u = sum(c * (c - 1) // 2 for c in sig.values()) / npairs
    ident_s = sum(c * (c - 1) // 2 for c in sigs.values()) / npairs
    zero_sig = sig.get(tuple([0] * k.shape[0]), 0) / M
    return out, {"identical_unstaged_sig_pairs": ident_u, "identical_staged_sig_pairs": ident_s,
                 "distinct_unstaged_sigs": len(sig), "all_zero_sig_primes": zero_sig}


def slope(Js, vals):
    x, y = np.log(np.array(Js, float)), np.log(np.array(vals, float))
    return float(np.polyfit(x, y, 1)[0])


def run_analyze():
    lab = load_labels()
    summary = []
    for (sname, fam, b), (primes, k) in sorted(lab.items()):
        Bmax = bound_B(b, F(1, 2))
        for beta in BETAS:
            B = bound_B(b, beta)
            kk = truncate(k, len(stage_steps(B, xbound(b))))
            t0 = time.time()
            ps, npairs = pair_stats(kk)
            pr, sig = prime_stats(primes, kk)
            vz = [ps[J]["Var_Zu_over_J"] * J for J in JS]
            vs = [ps[J]["Var_Zs_over_J"] * J for J in JS]
            vw = [pr[J]["Var_W_over_J"] * J for J in JS]
            rec = {"sample": sname, "family": fam, "b": b, "beta": str(beta), "B": B,
                   "M": int(len(primes)), "pairs": npairs, "pairs_stats": ps, "prime_stats": pr,
                   "signatures": sig,
                   "R_u": ps[64]["Var_Zu_over_J"] / ps[8]["Var_Zu_over_J"],
                   "R_s": ps[64]["Var_Zs_over_J"] / ps[8]["Var_Zs_over_J"],
                   "R_W": pr[64]["Var_W_over_J"] / pr[8]["Var_W_over_J"],
                   "slope_Var_Zu": slope(JS, vz), "slope_Var_Zs": slope(JS, vs),
                   "slope_Var_W": slope(JS, vw),
                   "failratio_s_64_8": (ps[64]["fail_s"] / ps[8]["fail_s"]) if ps[8]["fail_s"] else None}
            summary.append(rec)
            print(sname, fam, b, str(beta), f"R_u={rec['R_u']:.2f} R_s={rec['R_s']:.2f} "
                  f"R_W={rec['R_W']:.2f} fail_s8={ps[8]['fail_s']:.3g} fail_s64={ps[64]['fail_s']:.3g} "
                  f"EW/J={pr[64]['E_W_over_J']:.3f} cbar={pr[64]['mean_offdiag_cov']:.2e} "
                  f"({time.time() - t0:.1f}s)", flush=True)
    os.makedirs(os.path.dirname(SUMMARY), exist_ok=True)
    with open(SUMMARY, "w") as fh:
        json.dump(summary, fh, indent=1, default=float)


if __name__ == "__main__":
    {"labels": run_labels, "analyze": run_analyze}[sys.argv[1]]()
