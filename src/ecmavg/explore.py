"""EXPLORATORY (not pre-registered) analyses of the frozen labels.  HARNESS, factor-aware.

  .venv/bin/python -m src.ecmavg.explore      # writes experiments/analysis/ecmavg_explore.json

E1  permutation null for the mean off-diagonal covariance cbar of the exposure indicators
    X_j(p) (rows permuted independently across primes: keeps every curve's rate, kills
    cross-curve dependence).  z = (cbar - mean_null)/sd_null.
E2  the same after conditioning on p mod L (L = 120, 840): residual covariance of
    X_j(p) - mean_{p' ~ p} X_j(p') and a within-class permutation null.  If the common mode
    is local (small-prime Frobenius data), the within-class z should collapse.
    Class sets: p mod 120, p mod 840, 16 size buckets (rank of p in the sample), and
    (p mod 120) x 4 size buckets.
E3  normalised common mode cbar / s^2 and cbar / s^2 * (log2 B)^2 across b.
"""
from __future__ import annotations

import json
import os
from fractions import Fraction as F

import numpy as np

from .curves import bound_B, stage_steps
from .experiments import ROOT, load_labels, xbound
from .labels import INF, truncate

OUT = os.path.join(ROOT, "experiments", "analysis", "ecmavg_explore.json")


def cbar(X: np.ndarray) -> float:
    J = X.shape[0]
    C = np.cov(X, bias=True)
    return float((C.sum() - np.trace(C)) / (J * (J - 1)))


def residualise(X: np.ndarray, cls: np.ndarray) -> np.ndarray:
    R = X.copy()
    for c in np.unique(cls):
        m = cls == c
        R[:, m] -= X[:, m].mean(axis=1, keepdims=True)
    return R


def perm_null(X, cls, nperm, rng, within):
    vals = []
    for _ in range(nperm):
        Y = X.copy()
        for j in range(X.shape[0]):
            if within:
                for c in np.unique(cls):
                    idx = np.flatnonzero(cls == c)
                    Y[j, idx] = X[j, rng.permutation(idx)]
            else:
                Y[j] = X[j, rng.permutation(X.shape[1])]
        vals.append(cbar(residualise(Y, cls) if within else Y))
    return float(np.mean(vals)), float(np.std(vals))


def main(nperm: int = 100):
    lab = load_labels()
    rng = np.random.default_rng(20261018)
    res = []
    for (sname, fam, b), (primes, k) in sorted(lab.items()):
        if sname == "S1":
            continue
        for beta in (F(1, 3), F(1, 4)):
            B = bound_B(b, beta)
            kk = truncate(k, len(stage_steps(B, xbound(b))))
            X = (kk != INF).astype(np.float64)
            s = float(X.mean())
            c0 = cbar(X)
            m0, sd0 = perm_null(X, None, nperm, rng, False)
            rec = {"sample": sname, "family": fam, "b": b, "beta": str(beta), "B": B,
                   "M": int(len(primes)), "s": s, "cbar": c0, "z_cbar": (c0 - m0) / sd0,
                   "cbar_over_s2": c0 / s ** 2,
                   "cbar_over_s2_log2B_sq": c0 / s ** 2 * np.log2(B) ** 2}
            order = np.argsort(np.argsort(primes))
            size16 = (order * 16) // len(primes)
            size4 = (order * 4) // len(primes)
            for L, cls in (("120", primes % 120), ("840", primes % 840), ("size16", size16),
                           ("120xsize4", (primes % 120) * 4 + size4)):
                cw = cbar(residualise(X, cls))
                mw, sdw = perm_null(X, cls, nperm, rng, True)
                rec[f"cbar_within_{L}"] = cw
                rec[f"z_within_{L}"] = (cw - mw) / sdw
                rec[f"excess_within_{L}"] = cw - mw
            res.append(rec)
            print(sname, fam, b, str(beta), f"s={s:.3f} cbar={c0:.2e} z={rec['z_cbar']:.1f} "
                  f"c/s2*lg2B^2={rec['cbar_over_s2_log2B_sq']:.2f} "
                  f"z_w120={rec['z_within_120']:.1f} z_w840={rec['z_within_840']:.1f} "
                  f"exc840={rec['excess_within_840']:.1e} exc_size16={rec['excess_within_size16']:.1e} "
                  f"exc_120xs4={rec['excess_within_120xsize4']:.1e} z_120xs4={rec['z_within_120xsize4']:.1f}",
                  flush=True)
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()


def lower_tail(out_path: str = os.path.join(ROOT, "experiments", "analysis", "ecmavg_tail.json")):
    """E4: observed P(W_J = 0) versus the independent-curves model prod_j (1 - s_j), where s_j
    is curve j's empirical exposure rate; and the staged pair failure versus P(W=0)^2."""
    from .experiments import pair_stats
    lab = load_labels()
    res = []
    for (sname, fam, b), (primes, k) in sorted(lab.items()):
        if sname == "S1":
            continue
        for beta in (F(1, 3), F(1, 4)):
            B = bound_B(b, beta)
            kk = truncate(k, len(stage_steps(B, xbound(b))))
            X = (kk != INF)
            rec = {"sample": sname, "family": fam, "b": b, "beta": str(beta)}
            for J in (8, 16, 32, 64):
                obs = float((~X[:J].any(axis=0)).mean())
                ind = float(np.prod(1 - X[:J].mean(axis=1)))
                rec[f"PW0_{J}"] = obs
                rec[f"PW0_indep_{J}"] = ind
            res.append(rec)
            print(sname, fam, b, str(beta), " ".join(
                f"J={J}: {rec[f'PW0_{J}']:.2e}/{rec[f'PW0_indep_{J}']:.2e}" for J in (8, 16, 32, 64)),
                flush=True)
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=1)
