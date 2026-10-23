"""HARNESS for the newfamily battery (knows the factors).  Never imported by operational code.

Instance law D_n exactly as the eps battery (src/eps/harness.py): N = pq from balanced_semiprime(n, seed),
conditioned on N = 1 mod 4; labels eps = (-1/p), QR_2 = [2 square mod N] on N = 1 mod 8.  Seed ranges are
disjoint from every earlier run (notes/newfamily_prereg.md section 4).

Analysis reuses src/eps/analysis.py::strat_z (stratified permutation z, strata N mod 2^8).
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import time
from multiprocessing import Pool
from typing import Dict, List

import numpy as np

from src.eps.analysis import mde_rho, raw_z, strat_z
from src.eps.harness import instance, labels
from src.newfamily.battery import MENU, MENU_NAMES, S_BONFERRONI, compute_stats

SEED_BASE = {
    "pilot": 24_500_000_000_000,   # unit tests / pipeline pilot (labels of pilot not used for any decision)
    "main": 24_600_000_000_000,    # the battery
    "fresh": 24_700_000_000_000,   # replication of hits
}
CHUNK = 250
Z_THRESH = 6.0
Z_REPLICATE = 4.0


def _work(args) -> List[Dict]:
    n, s0, s1 = args
    out = []
    for seed in range(s0, s1):
        g = instance(n, seed)
        if g is None:
            continue
        t0 = time.perf_counter()
        st, aux = compute_stats(g["N"])          # operational: N only
        dt = time.perf_counter() - t0
        lab = labels(g["p"], g["q"])
        out.append({"seed": seed, "N": g["N"], "p": g["p"], "q": g["q"], "eps": lab["eps"], "qr2": lab["qr2"],
                    "x": [st[k] for k in MENU_NAMES], "t": round(dt, 6), "nonunits": aux["nonunits"]})
    return out


def run(n: int, count: int, base: int, out: str, workers: int) -> Dict:
    t0 = time.time()
    recs: List[Dict] = []
    nxt = 0
    with Pool(workers) as pool:
        while len(recs) < count:
            nch = max(workers, int((count - len(recs)) * 2.1 / CHUNK) + 1)
            jobs = [(n, base + (nxt + i) * CHUNK, base + (nxt + i + 1) * CHUNK) for i in range(nch)]
            nxt += nch
            for part in pool.imap(_work, jobs):
                recs.extend(part)
    recs = recs[:count]
    header = {"header": True, "n": n, "count": count, "seed_base": base, "seed_last": recs[-1]["seed"],
              "menu": MENU_NAMES, "wall_s": round(time.time() - t0, 1),
              "cpu_s": round(sum(r["t"] for r in recs), 1), "workers": workers}
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with gzip.open(out, "wt") as fh:
        fh.write(json.dumps(header) + "\n")
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    return header


def load(path: str) -> Dict:
    with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline())
        assert header["menu"] == MENU_NAMES, "menu mismatch"
        recs = [json.loads(line) for line in fh]
    return {"header": header, "X": np.array([r["x"] for r in recs], dtype=np.float64),
            "Nmod": np.array([r["N"] % 256 for r in recs], dtype=np.int64),
            "eps": np.array([r["eps"] for r in recs], dtype=np.int64),
            "qr2": np.array([-1 if r["qr2"] is None else r["qr2"] for r in recs], dtype=np.int64),
            "nonunits": [r["nonunits"] for r in recs], "t": np.array([r["t"] for r in recs])}


def pc_check(x: np.ndarray, y01: np.ndarray) -> Dict:
    fires = x == 1
    viol = int((fires & (y01 == 0)).sum())
    return {"fires": int(fires.sum()), "violations": viol}


def analyze(d: Dict, target: str, planted_seed: int = 0) -> Dict:
    if target == "eps":
        m = np.ones(len(d["eps"]), dtype=bool)
        y = (d["eps"] == 1).astype(np.float64)
    else:
        m = d["qr2"] >= 0
        y = d["qr2"][m].astype(np.float64)
    X = d["X"][m]
    h = d["Nmod"][m]
    res = {"target": target, "n_inst": int(m.sum()), "stats": {}}
    for j, (name, grp, role) in enumerate(MENU):
        x = X[:, j]
        if role.startswith("pc"):
            if (role == "pc_eps") == (target == "eps"):
                res["stats"][name] = {"role": role, **pc_check(x, y), **strat_z(x, y, h)}
            continue
        z = strat_z(x, y, h)
        res["stats"][name] = {"role": role, "group": grp, "z": z["z"], "degenerate": z["degenerate"],
                              "raw_z": raw_z(x, y), "hit": (not z["degenerate"]) and abs(z["z"]) >= Z_THRESH}
    # planted control: x = y with prob (1+rho)/2 flipped noise, rho = 0.05
    rng = np.random.default_rng(planted_seed)
    rho = 0.05
    flip = rng.random(len(y)) < (1 - rho) / 2
    xp = np.where(flip, 1 - y, y)
    res["planted"] = {"rho": rho, "z": strat_z(xp, y, h)["z"], "z_expected": rho * math.sqrt(len(y))}
    # within-stratum shuffle null: mean z^2 over non-degenerate menu stats with y shuffled within strata
    ys = y.copy()
    for hv in np.unique(h):
        idx = np.where(h == hv)[0]
        ys[idx] = ys[rng.permutation(idx)]
    zs = []
    for j, (name, grp, role) in enumerate(MENU):
        if role == "menu":
            z = strat_z(X[:, j], ys, h)
            if not z["degenerate"]:
                zs.append(z["z"] ** 2)
    res["shuffle_mean_z2"] = float(np.mean(zs))
    menu_z = [v["z"] for v in res["stats"].values() if v.get("role") == "menu" and not v["degenerate"]]
    res["menu_mean_z2"] = float(np.mean(np.square(menu_z)))
    res["menu_max_abs_z"] = float(np.max(np.abs(menu_z)))
    res["mde_rho"] = mde_rho(int(m.sum()) - 256, Z_THRESH)
    res["hits"] = [k for k, v in res["stats"].items() if v.get("hit")]
    res["degenerate"] = [k for k, v in res["stats"].items() if v.get("degenerate")]
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", required=True)
    ap.add_argument("--count", type=int, default=100000)
    ap.add_argument("--range", choices=sorted(SEED_BASE), required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--analyze", action="store_true")
    a = ap.parse_args(argv)
    for n in a.n:
        out = f"experiments/raw/newfamily_{a.range}_n{n}.jsonl.gz"
        if not a.analyze:
            h = run(n, a.count, SEED_BASE[a.range] + n * 10 ** 9, out, a.workers)
            print(json.dumps({k: v for k, v in h.items() if k != "menu"}), flush=True)
        d = load(out)
        summ = {"n": n, "range": a.range, "eps": analyze(d, "eps"), "qr2": analyze(d, "qr2"),
                "nonunit_instances": sum(1 for x in d["nonunits"] if x),
                "mean_ms": float(d["t"].mean() * 1000)}
        os.makedirs("experiments/analysis", exist_ok=True)
        with open(f"experiments/analysis/newfamily_{a.range}_n{n}.json", "w") as fh:
            json.dump(summ, fh, indent=1)
        for t in ("eps", "qr2"):
            r = summ[t]
            pcs = {k: (v["fires"], v["violations"], round(v["z"], 1)) for k, v in r["stats"].items()
                   if v.get("role", "").startswith("pc")}
            print(n, t, "hits", r["hits"], "max|z|", round(r["menu_max_abs_z"], 2), "mean z2",
                  round(r["menu_mean_z2"], 2), "shuffle", round(r["shuffle_mean_z2"], 2), "planted",
                  round(r["planted"]["z"], 1), "/", round(r["planted"]["z_expected"], 1), "PC", pcs,
                  "degenerate", r["degenerate"], flush=True)


if __name__ == "__main__":
    main()
