"""sihard experiments E1-E4 (pre-registered in notes/sihard_prereg.md, commits 7f9e398 + bb7fbd0).

Usage: .venv/bin/python -m src.sihard.run_experiments e1|e2|e3|e4|analyze
Raw output: experiments/raw/sihard_E{1,2,3,4}.jsonl (HARNESS; E1 uses instance labels only for bookkeeping).
"""
from __future__ import annotations

import json
import math
import os
import random
import statistics
import sys
import time

import numpy as np

from src.common.instances import balanced_semiprime
from src.sihard.core import (dickman_rho, johnson_S, johnson_counts, merge_bohr_counts,
                             merge_interval_average, primes_upto, smooth_in_range, smooth_mask)
from src.sihard.reduction import build_instance, count_exact_covers, count_S_smooth_in_window, random_x3c

RAW = "experiments/raw"


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, f"sihard_{name}.jsonl"), "w")


# ------------------------------------------------------------------ E1
def e1():
    f = _out("E1")
    for y in (1 << 10, 1 << 12):
        ps = primes_upto(y)
        for mu in (2, 8, 32):
            for seed in range(10):
                N = balanced_semiprime(48, seed).N
                u_mid = math.log(100.5 * N) / math.log(y)
                h = round(mu / dickman_rho(u_mid))
                near = [int(smooth_mask(k * N, h + 1, ps).sum()) for k in range(1, 201)]
                rng = random.Random(f"sihard-E1-{seed}-{y}-{mu}")
                ctrl = []
                for _ in range(200):
                    x = rng.randrange(N, 201 * N)
                    ctrl.append(int(smooth_mask(x, h + 1, ps).sum()))
                rec = {"y": y, "mu": mu, "seed": seed, "N": N, "h": h, "u_mid": u_mid,
                       "near": near, "ctrl": ctrl}
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print("E1", y, mu, seed, h, sum(near) / 200, sum(ctrl) / 200, flush=True)
    f.close()


# ------------------------------------------------------------------ E2
def _p_family(y, V, rng):
    ps = primes_upto(y)
    t = math.log(V)
    fam = []
    fam.append(("divisor_trick", 1 << (V.bit_length() - 1), 1.0))
    for c in (1.5, 2, 4, 8):
        P, i = 1, 0
        while math.log(P) < c * t and i < len(ps):
            P *= ps[i]
            i += 1
        if math.log(P) >= c * t:
            fam.append((f"primorial_c{c}", P, c))
    for c in (2, 4, 8, 16):
        z, P = 1, 1
        while math.log(P) < c * t and z < y:
            z += 1
            P = math.lcm(P, z)
        if math.log(P) >= c * t:
            fam.append((f"lcm_c{c}", P, c))
    for r in range(5):
        sub = ps[:]
        rng.shuffle(sub)
        P, i = 1, 0
        while math.log(P) < 4 * t:
            P *= sub[i]
            i += 1
        fam.append((f"random{r}_c4", P, 4))
    return fam


def e2():
    f = _out("E2")
    X = 1 << 40
    for cell, y, M in (("E2", 1 << 8, 1 << 16), ("E2b", 1 << 10, 1 << 20)):
        ps = primes_upto(y)
        for V in (1 << 8, 1 << 10, 1 << 12):
            sm = smooth_in_range(X, X + M + V + 1, y)
            dens = len(sm) / (M + V + 1)
            rng = random.Random(f"sihard-E2-{cell}-{V}")
            for name, P, c in _p_family(y, V, rng):
                S = johnson_S(V, P)
                r = johnson_counts(X, M, V, P, S, sm, ps)
                rec = {"cell": cell, "y": y, "M": M, "V": V, "P_name": name, "logP": math.log(P), "c": c,
                       "logS": math.log(S), "dens": dens, "n_smooth": len(sm), **r}
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print("E2", cell, V, name, r, round(dens, 6), flush=True)
    f.close()


# ------------------------------------------------------------------ E3
def e3():
    f = _out("E3")
    lo, hi = 1 << 18, int(math.isqrt(2 * (1 << 36)))  # [2^18, 2^18 sqrt2)
    pool = smooth_in_range(lo, hi, 1 << 8)
    X = 1 << 36
    for h in (1 << 6, 1 << 10):
        rng = random.Random(f"sihard-E3-{h}")
        A = rng.sample(pool, 2000)
        B = rng.sample(pool, 2000)
        meas, bound = merge_interval_average(A, B, X, h)
        rec = {"version": "interval", "h": h, "measured": meas, "bound": bound, "ratio": meas / bound}
        f.write(json.dumps(rec) + "\n")
        print("E3", rec, flush=True)
        Ns = [rng.randrange(1 << 30, 1 << 31) for _ in range(4096)]
        cnt = merge_bohr_counts(A, B, Ns, h)
        maxprod = max(A) * max(B)
        K = -(-maxprod // (1 << 30))
        bnd = 2 * (h + 1) * (1 + math.log(K)) * len(A) * len(B) / (1 << 30)
        mean = statistics.mean(cnt)
        rec = {"version": "bohr", "h": h, "measured": mean, "sd": statistics.pstdev(cnt), "bound": bnd,
               "K": K, "ratio": mean / bnd}
        f.write(json.dumps(rec) + "\n")
        print("E3", rec, flush=True)
    f.close()


# ------------------------------------------------------------------ E4
def e4(name="E4", R=1 << 20):
    f = _out(name)
    rng = random.Random("sihard-" + name)
    for i in range(300):
        q = rng.choice((2, 3, 4))
        m = q + rng.randint(1, 4)
        sets = random_x3c(q, m, rng, plant=rng.random() < 0.5)
        inst = build_instance(q, sets, rng, R=R)
        c_smooth = count_S_smooth_in_window(inst["primes"], inst["x"], inst["h"])
        c_cov = count_exact_covers(q, inst["sets"])
        rec = {"i": i, "q": q, "m": len(inst["sets"]), "smooth": c_smooth, "covers": c_cov,
               "log2x": inst["x"].bit_length(), "log2h": inst["h"].bit_length(),
               "log2y": max(inst["primes"]).bit_length()}
        f.write(json.dumps(rec) + "\n")
    f.close()


# ------------------------------------------------------------------ analysis
def _load(name):
    with open(os.path.join(RAW, f"sihard_{name}.jsonl")) as fh:
        return [json.loads(l) for l in fh]


def analyze():
    out = {}
    if os.path.exists(os.path.join(RAW, "sihard_E1.jsonl")):
        rows = _load("E1")
        cells = {}
        for r in rows:
            cells.setdefault((r["y"], r["mu"]), []).append(r)
        for (y, mu), rs in sorted(cells.items()):
            near = [c for r in rs for c in r["near"]]
            ctrl = [c for r in rs for c in r["ctrl"]]
            mn, mc = statistics.mean(near), statistics.mean(ctrl)
            idn = statistics.pvariance(near) / mn if mn else float("nan")
            idc = statistics.pvariance(ctrl) / mc if mc else float("nan")
            e_n = sum(1 for c in near if c == 0) / len(near)
            e_c = sum(1 for c in ctrl if c == 0) / len(ctrl)
            P1 = 0.9 <= mn / mc <= 1.1
            P2 = (0.7 <= idn <= 1.6) and abs(idn - idc) <= 0.3
            P3 = abs(e_n - e_c) <= 0.05 if mu in (2, 8) else None
            out[f"E1 y={y} mu={mu}"] = dict(mean_near=mn, mean_ctrl=mc, ratio=mn / mc, ID_near=idn, ID_ctrl=idc,
                                            empty_near=e_n, empty_ctrl=e_c, ctrl_over_dickman=mc / mu,
                                            P1=P1, P2=P2, P3=P3)
    if os.path.exists(os.path.join(RAW, "sihard_E2.jsonl")):
        rows = _load("E2")
        P4 = all(r["A"] <= r["B"] + 1e-12 and r["B"] <= r["Bd"] + 1e-12 for r in rows)
        best = {}
        for r in rows:
            key = (r["cell"], r["V"])
            g = r["A"] / r["dens"] if r["dens"] else 0
            if r["P_name"] == "divisor_trick":
                best.setdefault(key, {})["divisor_trick"] = g
            b = best.setdefault(key, {})
            if g > b.get("best", (-1, ""))[0]:
                b["best"] = (g, r["P_name"])
        out["E2 P4"] = P4
        out["E2 best A/dens"] = {f"{k[0]} V={k[1]}": v for k, v in best.items()}
    if os.path.exists(os.path.join(RAW, "sihard_E3.jsonl")):
        rows = _load("E3")
        out["E3"] = rows
        out["E3 P5"] = all(r["measured"] <= r["bound"] for r in rows)
        out["E3 P6"] = all(0.9 <= r["ratio"] <= 1.0 for r in rows if r["version"] == "interval")
    for nm in ("E4", "E4b"):
        if os.path.exists(os.path.join(RAW, f"sihard_{nm}.jsonl")):
            rows = _load(nm)
            out[f"{nm} P7"] = all(r["smooth"] == r["covers"] for r in rows)
            out[f"{nm} n"] = len(rows)
            out[f"{nm} with covers"] = sum(1 for r in rows if r["covers"] > 0)
            out[f"{nm} max covers"] = max(r["covers"] for r in rows)
            out[f"{nm} max log2x"] = max(r["log2x"] for r in rows)
    os.makedirs("experiments/analysis", exist_ok=True)
    with open("experiments/analysis/sihard_summary.json", "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    t0 = time.time()
    {"e1": e1, "e2": e2, "e3": e3, "e4": e4, "e4b": lambda: e4("E4b", None), "analyze": analyze}[sys.argv[1]]()
    print("elapsed", round(time.time() - t0, 1), file=sys.stderr)
