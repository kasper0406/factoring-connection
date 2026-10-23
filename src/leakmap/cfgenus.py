"""E2: does the CF prefix of sqrt(N) leak only genus information?  (test of conjecture H-GEN)

Pre-registration: notes/leakmap_prereg.md E2.
- `cf_stats(N)` (OPERATIONAL, input N only): the 16 CF-prefix menu statistics, imported read-only from
  src/eps/cf.py, plus first indices of the one-sided controls PC{l}_even / PC{l}_odd (l = 3, 5, 7):
  some even i <= T with Q_i = l k^2 (gcd(k, N) = 1) certifies (l|p) = (l|q) = +1;
  some odd i <= T with Q_i = l k^2 certifies (-l|p) = +1  (since A_{i-1}^2 = (-1)^i Q_i mod N).
- `labels(p, q)` (HARNESS): genus, non-genus and archimedean targets.
- `analyse(...)`: stratified permutation z and two-way cluster-robust z.

Usage: .venv/bin/python -m src.leakmap.cfgenus run --n 32 --count 200000 --range main
       .venv/bin/python -m src.leakmap.cfgenus analyse experiments/raw/leakmap_e2_main_n32.jsonl.gz
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from math import gcd, sqrt
from multiprocessing import Pool
from typing import Dict, List

import gmpy2
import numpy as np

from src.eps.cf import cf_features
from src.common.instances import balanced_semiprime

STATS = [f"cf{T}_{k}" for T in (64, 1024) for k in
         ("n_a1", "n_a2", "sumbits", "maxbits", "n_Qeven", "n_sq_even", "n_Q1mod4_odd", "n_2sq_odd")]
ELLS = (3, 5, 7)
T_PC = 1024
RANGES = {"main": 9_500_000_000_000, "fresh": 9_600_000_000_000}


# ----------------------------------------------------------------------------- OPERATIONAL
def pc_ell(N: int, T: int = T_PC) -> Dict[str, int]:
    a0 = int(gmpy2.isqrt(N))
    P, Q, a = 0, 1, a0
    first = {f"pc{l}_{par}": 0 for l in ELLS for par in ("even", "odd")}
    for i in range(1, T + 1):
        P = a * Q - P
        Q = (N - P * P) // Q
        a = (a0 + P) // Q
        par = "odd" if i & 1 else "even"
        for l in ELLS:
            if Q % l == 0:
                r = Q // l
                if gmpy2.is_square(r):
                    k = int(gmpy2.isqrt(r))
                    key = f"pc{l}_{par}"
                    if first[key] == 0 and gcd(k, N) == 1:
                        first[key] = i
    return first


def cf_stats(N: int) -> Dict[str, int]:
    f = cf_features(N)
    out = {k: f[k] for k in STATS}
    out.update(pc_ell(N))
    return out


# ----------------------------------------------------------------------------- HARNESS
def labels(p: int, q: int) -> Dict:
    N = p * q
    return {
        "G1": int(p % 4 == 1) if N % 4 == 1 else None,
        "G2": int(p % 3 == 1) if N % 3 == 1 else None,
        "G3": int(p % 5 in (1, 4)) if N % 5 == 1 else None,
        "G4": int(p % 7 in (1, 2, 4)) if N % 7 == 1 else None,
        "NG1": int(p % 7 in (1, 6)) if N % 7 == 1 else None,
        "NG2": int(p % 5 == 1) if (N % 5 == 1 and p % 5 in (1, 4)) else None,
        "A1": int(32 * q < 33 * p),
    }


TARGET_MOD = {"G1": 1, "G2": 3, "G3": 5, "G4": 7, "NG1": 7, "NG2": 5, "A1": 1}


def _work(args):
    n, lo, hi = args
    out = []
    for seed in range(lo, hi):
        inst = balanced_semiprime(n, seed)
        p, q = sorted(inst.labels["factors"]) if "factors" in inst.labels else (inst.labels["p"], inst.labels["q"])
        N = inst.N
        rec = {"seed": seed, "N": N, "p": p, "q": q}
        rec.update(cf_stats(N))
        out.append(rec)
    return out


def run(n: int, count: int, rng_name: str, workers: int, out: str) -> None:
    base = RANGES[rng_name] + n * 1_000_000_000
    chunk = 2000
    jobs = [(n, base + s, base + min(s + chunk, count)) for s in range(0, count, chunk)]
    with gzip.open(out, "wt") as f, Pool(workers) as pool:
        for recs in pool.imap(_work, jobs):
            for r in recs:
                f.write(json.dumps(r) + "\n")


# ----------------------------------------------------------------------------- analysis
def strat_z(x: np.ndarray, y: np.ndarray, s: np.ndarray) -> float:
    """Exact within-stratum permutation z of sum x*y (y binary)."""
    T = float(np.dot(x, y))
    E = 0.0
    V = 0.0
    order = np.argsort(s, kind="stable")
    xs, ys, ss = x[order], y[order], s[order]
    cuts = np.flatnonzero(np.diff(ss)) + 1
    for xa, ya in zip(np.split(xs, cuts), np.split(ys, cuts)):
        m = len(xa)
        if m < 2:
            continue
        n1 = float(ya.sum())
        mu = xa.mean()
        E += n1 * mu
        V += n1 * (m - n1) / (m * (m - 1)) * float(((xa - mu) ** 2).sum())
    if V <= 0:
        return 0.0
    return (T - E) / sqrt(V)


def cluster_z(x: np.ndarray, y: np.ndarray, s: np.ndarray, cp: np.ndarray, cq: np.ndarray) -> float:
    """Two-way (p, q) cluster-robust z for the within-stratum covariance of x and y."""
    xt = x.astype(float).copy()
    yt = y.astype(float).copy()
    for arr in (xt, yt):
        sums = defaultdict(float)
        cnts = defaultdict(int)
        for v, k in zip(arr, s):
            sums[k] += v
            cnts[k] += 1
        means = np.array([sums[k] / cnts[k] for k in s])
        arr -= means
    u = xt * yt
    S = u.sum()

    def csum(c):
        d = defaultdict(float)
        for v, k in zip(u, c):
            d[k] += v
        return sum(w * w for w in d.values())

    V = csum(cp) + csum(cq) - float((u * u).sum())
    if V <= 0:
        return 0.0
    return float(S / sqrt(V))


def mde_rho(neff: int, z: float = 5.0, power_z: float = 0.84) -> float:
    return (z + power_z) / sqrt(max(neff, 1))


def analyse(path: str, cluster: bool) -> Dict:
    recs = [json.loads(l) for l in gzip.open(path, "rt")]
    res = {"file": path, "count": len(recs), "targets": {}, "controls": {}}
    for tname, mod in TARGET_MOD.items():
        rows = []
        for r in recs:
            lab = labels(r["p"], r["q"])[tname]
            if lab is not None:
                rows.append((r, lab))
        y = np.array([lab for _, lab in rows], dtype=float)
        N = [r["N"] for r, _ in rows]
        if tname == "A1":
            s = np.array([(Nv % 256) * 64 + (Nv >> (Nv.bit_length() - 6)) % 64 for Nv in N])
        else:
            s = np.array([Nv % (256 * mod) for Nv in N])
        cp = np.array([r["p"] for r, _ in rows])
        cq = np.array([r["q"] for r, _ in rows])
        out = {"n_rows": len(rows), "mean_target": float(y.mean()), "mde_rho": mde_rho(len(rows)), "stats": {}}
        for st in STATS:
            x = np.array([r[st] for r, _ in rows], dtype=float)
            z = strat_z(x, y, s)
            rho = float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 else 0.0
            d = {"z": round(z, 3), "rho": round(rho, 4)}
            if cluster:
                d["z_cr"] = round(cluster_z(x, y, s, cp, cq), 3)
            out["stats"][st] = d
        res["targets"][tname] = out
    # one-sided positive controls
    for l in ELLS:
        for par in ("even", "odd"):
            key = f"pc{l}_{par}"
            fires = viol = 0
            for r in recs:
                if r[key]:
                    fires += 1
                    sym = int(gmpy2.legendre(l if par == "even" else -l, r["p"]))
                    sym2 = int(gmpy2.legendre(l if par == "even" else -l, r["q"]))
                    if sym != 1 or sym2 != 1:
                        viol += 1
            res["controls"][key] = {"fires": fires, "violations": viol}
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd")
    r = sub.add_parser("run")
    r.add_argument("--n", type=int, required=True)
    r.add_argument("--count", type=int, default=200000)
    r.add_argument("--range", default="main")
    r.add_argument("--workers", type=int, default=6)
    r.add_argument("--out", default=None)
    an = sub.add_parser("analyse")
    an.add_argument("path")
    an.add_argument("--cluster", action="store_true")
    an.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        out = a.out or f"experiments/raw/leakmap_e2_{a.range}_n{a.n}.jsonl.gz"
        run(a.n, a.count, a.range, a.workers, out)
    elif a.cmd == "analyse":
        res = analyse(a.path, a.cluster)
        s = json.dumps(res, indent=1)
        if a.out:
            open(a.out, "w").write(s)
        print(s)


if __name__ == "__main__":
    main()
