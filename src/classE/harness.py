"""HARNESS for the classE battery (knows the factors).  Never imported by operational code.

Instance law D_n exactly as eps/newfamily: N = pq from balanced_semiprime(n, seed) (p, q with the top two bits
set, p < q), conditioned on N = 1 (mod 4).  Targets (notes/classE_prereg.md section 1):
  eps = (-1/p) on all of D_n;
  tb  = [8p >= 7q] (the first non-constant bit of p/q, which lies in (3/4, 1)).
Floors: eps -> strata N mod 2^8;  tb -> strata top 12 bits of N (registered primary); tb with N mod 2^8 strata is
reported as a secondary (size-confounded) view only.
Analysis reuses src/eps/analysis.py::strat_z.
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
from src.classE.battery import MENU, MENU_NAMES, compute_stats

SEED_BASE = {
    "pilot": 25_510_000_000_000,   # label-free degeneracy pilot and unit tests
    "main": 25_520_000_000_000,    # the battery
    "fresh": 25_530_000_000_000,   # replication of hits
}
CHUNK = 250
Z_THRESH = 6.0
Z_REPLICATE = 4.0
TARGETS = ("eps", "tb")


def tb_label(p: int, q: int) -> int:
    return int(8 * p >= 7 * q)


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
        out.append({"seed": seed, "N": g["N"], "p": g["p"], "q": g["q"], "eps": lab["eps"],
                    "tb": tb_label(g["p"], g["q"]),
                    "x": [st[k] for k in MENU_NAMES], "t": round(dt, 6), "nonunits": aux["nonunits"],
                    "path_hits": aux["path_hits"]})
    return out


def run(n: int, count: int, base: int, out: str, workers: int) -> Dict:
    """Streams records to a temporary file (bounded memory), then writes header + records."""
    t0 = time.time()
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    tmp = out + ".part"
    written = 0
    cpu = 0.0
    last_seed = None
    nxt = 0
    with gzip.open(tmp, "wt") as fh, Pool(workers) as pool:
        while written < count:
            nch = max(workers, int((count - written) * 2.1 / CHUNK) + 1)
            jobs = [(n, base + (nxt + i) * CHUNK, base + (nxt + i + 1) * CHUNK) for i in range(nch)]
            nxt += nch
            for part in pool.imap(_work, jobs):
                for r in part:
                    if written >= count:
                        break
                    fh.write(json.dumps(r) + "\n")
                    written += 1
                    cpu += r["t"]
                    last_seed = r["seed"]
    header = {"header": True, "n": n, "count": count, "seed_base": base, "seed_last": last_seed,
              "menu": MENU_NAMES, "wall_s": round(time.time() - t0, 1), "cpu_s": round(cpu, 1),
              "workers": workers}
    with gzip.open(out, "wt") as fo, gzip.open(tmp, "rt") as fi:
        fo.write(json.dumps(header) + "\n")
        for line in fi:
            fo.write(line)
    os.remove(tmp)
    return header


def load(path: str) -> Dict:
    with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline())
        assert header["menu"] == MENU_NAMES, "menu mismatch"
        n = header["n"]
        X = np.zeros((header["count"], len(MENU_NAMES)), dtype=np.float64)
        Nmod, Ntop, eps, tb, t = [], [], [], [], []
        nonunits, path_hits = [], []
        i = 0
        for line in fh:
            r = json.loads(line)
            X[i] = r["x"]
            Nmod.append(r["N"] % 256)
            Ntop.append(r["N"] >> (n - 12))
            eps.append(r["eps"])
            tb.append(r["tb"])
            t.append(r["t"])
            nonunits.append(r["nonunits"])
            path_hits.append(r["path_hits"])
            i += 1
    return {"header": header, "X": X[:i], "Nmod": np.array(Nmod, dtype=np.int64),
            "Ntop": np.array(Ntop, dtype=np.int64), "eps": np.array(eps, dtype=np.int64),
            "tb": np.array(tb, dtype=np.int64), "nonunits": nonunits, "path_hits": path_hits,
            "t": np.array(t)}


def pc_check(x: np.ndarray, y01: np.ndarray) -> Dict:
    fires = x == 1
    return {"fires": int(fires.sum()), "violations": int((fires & (y01 == 0)).sum())}


def analyze(d: Dict, target: str, strata_key: str, planted_seed: int = 0) -> Dict:
    y = (d["eps"] == 1).astype(np.float64) if target == "eps" else d["tb"].astype(np.float64)
    X = d["X"]
    h = d[strata_key]
    res = {"target": target, "strata": strata_key, "n_inst": int(len(y)), "stats": {}}
    for j, (name, grp, role) in enumerate(MENU):
        x = X[:, j]
        if role.startswith("pc"):
            if (role == "pc_eps") == (target == "eps"):
                res["stats"][name] = {"role": role, **pc_check(x, y), **strat_z(x, y, h)}
            continue
        z = strat_z(x, y, h)
        res["stats"][name] = {"role": role, "group": grp, "z": z["z"], "degenerate": z["degenerate"],
                              "raw_z": raw_z(x, y),
                              "hit": role.startswith("menu") and (not z["degenerate"]) and abs(z["z"]) >= Z_THRESH}
    rng = np.random.default_rng(planted_seed)
    rho = 0.05
    flip = rng.random(len(y)) < (1 - rho) / 2
    xp = np.where(flip, 1 - y, y)
    res["planted"] = {"rho": rho, "z": strat_z(xp, y, h)["z"], "z_expected": rho * math.sqrt(len(y))}
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
    for role in ("menu", "menuG", "nc"):
        zz = [v["z"] for v in res["stats"].values() if v.get("role") == role and not v["degenerate"]]
        res[f"{role}_mean_z2"] = float(np.mean(np.square(zz))) if zz else None
        res[f"{role}_max_abs_z"] = float(np.max(np.abs(zz))) if zz else None
    res["mde_rho"] = mde_rho(int(len(y)) - len(np.unique(h)), Z_THRESH)
    res["hits"] = [k for k, v in res["stats"].items() if v.get("hit") and v["role"] == "menu"]
    res["hitsG"] = [k for k, v in res["stats"].items() if v.get("hit") and v["role"] == "menuG"]
    res["degenerate"] = [k for k, v in res["stats"].items() if v.get("degenerate")]
    return res


def summarize(n: int, rng_name: str, d: Dict) -> Dict:
    summ = {"n": n, "range": rng_name,
            "eps": analyze(d, "eps", "Nmod"),
            "tb": analyze(d, "tb", "Ntop"),
            "tb_secondary_Nmod": analyze(d, "tb", "Nmod"),
            "nonunit_instances": sum(1 for x in d["nonunits"] if x),
            "path_hit_instances": sum(1 for x in d["path_hits"] if x),
            "mean_ms": float(d["t"].mean() * 1000),
            "p_eps": float((d["eps"] == 1).mean()), "p_tb": float(d["tb"].mean())}
    return summ


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", required=True)
    ap.add_argument("--count", type=int, default=200000)
    ap.add_argument("--range", choices=sorted(SEED_BASE), required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--analyze", action="store_true")
    a = ap.parse_args(argv)
    for n in a.n:
        out = f"experiments/raw/classE_{a.range}_n{n}.jsonl.gz"
        if not a.analyze:
            h = run(n, a.count, SEED_BASE[a.range] + n * 10 ** 9, out, a.workers)
            print(json.dumps({k: v for k, v in h.items() if k != "menu"}), flush=True)
        d = load(out)
        summ = summarize(n, a.range, d)
        os.makedirs("experiments/analysis", exist_ok=True)
        with open(f"experiments/analysis/classE_{a.range}_n{n}.json", "w") as fh:
            json.dump(summ, fh, indent=1)
        for t in ("eps", "tb", "tb_secondary_Nmod"):
            r = summ[t]
            pcs = {k: (v["fires"], v["violations"], round(v["z"], 1)) for k, v in r["stats"].items()
                   if v.get("role", "").startswith("pc")}
            print(n, t, "hits", r["hits"], "hitsG", r["hitsG"], "max|z| menu", r["menu_max_abs_z"],
                  "menuG", r["menuG_max_abs_z"], "mean z2", round(r["menu_mean_z2"], 2),
                  "shuffle", round(r["shuffle_mean_z2"], 2), "planted", round(r["planted"]["z"], 1), "/",
                  round(r["planted"]["z_expected"], 1), "PC", pcs, "ndeg", len(r["degenerate"]), flush=True)


if __name__ == "__main__":
    main()
