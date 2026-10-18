"""HARNESS for the eps/QR_2 battery (knows the factors).  Never imported by operational code.

Instance distribution D_n (NEXT_STEPS.md): N = pq from src.common.instances.balanced_semiprime(n, seed),
conditioned on N = 1 (mod 4).  Seeds s = base + j, j = 0, 1, 2, ...; the first `count` seeds (in order)
whose N is 1 mod 4 form the sample.  Labels: eps = (-1/p) = (-1/q) in {+1,-1};
QR_2 = [2 is a square mod N] in {0,1} (recorded only when N = 1 mod 8, where (2/p) = (2/q)).

Seed ranges (disjoint from all earlier runs, which used seeds < 10^6; see notes/eps_prereg.md section 4):
SEED_BASE below.  The supplementary curve uses base 'supp' + n * 10^9 for size n.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from multiprocessing import Pool
from typing import Dict, Iterator, List, Optional

from src.common.instances import balanced_semiprime
from src.eps.cf import cf_features
from src.eps.operational import MENU_NAMES, compute_stats

SEED_BASE = {
    "bench": 9_100_000_000_000,   # cost benchmark only; never reused
    "main": 9_200_000_000_000,    # the battery
    "fresh": 9_300_000_000_000,   # fresh-seed re-runs of hits (leakage audit)
    "supp": 9_400_000_000_000,    # supplementary positive-control frequency curve (other n)
    "test": 9_500_000_000_000,    # unit tests
}
CHUNK = 500


def labels(p: int, q: int) -> Dict:
    N = p * q
    assert p % 4 == q % 4, "D_n requires p = q (mod 4)"
    eps = 1 if p % 4 == 1 else -1
    out = {"eps": eps, "qr2": None}
    if N % 8 == 1:
        a, b = p % 8 in (1, 7), q % 8 in (1, 7)
        assert a == b
        out["qr2"] = int(a)
    return out


def instance(n: int, seed: int) -> Optional[Dict]:
    inst = balanced_semiprime(n, seed)
    if inst.N % 4 != 1:
        return None
    p, q = inst.labels["factors"]
    return {"N": inst.N, "p": p, "q": q}


def _work(args) -> List[Dict]:
    n, s0, s1 = args
    out = []
    for seed in range(s0, s1):
        g = instance(n, seed)
        if g is None:
            continue
        N = g["N"]
        t0 = time.perf_counter()
        st, aux = compute_stats(N)          # operational: sees N only
        dt = time.perf_counter() - t0
        lab = labels(g["p"], g["q"])
        rec = {"seed": seed, "N": N, "eps": lab["eps"], "qr2": lab["qr2"], "p": g["p"], "q": g["q"],
               "x": [st[k] for k in MENU_NAMES], "t": round(dt, 6),
               "tg": {k: round(v, 6) for k, v in aux["t"].items()},
               "first": {k[6:]: v for k, v in aux.items() if k.startswith("first_")},
               "nsmooth": aux["cf_nsmooth"], "mr_split": aux["mr_split"], "cf_nonunit": aux["cf_nonunit"], "cf_fb_ok": aux["cf_fb_ok"],
               "gcd_checks": aux["gcd_checks"],
               "nonunits": [[w, d] for w, d in aux["nonunits"]]}
        out.append(rec)
    return out


def run(n: int, count: int, base: int, out: str, workers: int, chunk: int = CHUNK) -> Dict:
    t0 = time.time()
    recs: List[Dict] = []
    next_chunk = 0
    with Pool(workers) as pool:
        while len(recs) < count:
            need = count - len(recs)
            nchunks = max(workers, int(need * 2.1 / chunk) + 1)
            jobs = [(n, base + (next_chunk + i) * chunk, base + (next_chunk + i + 1) * chunk)
                    for i in range(nchunks)]
            next_chunk += nchunks
            for part in pool.imap(_work, jobs):
                recs.extend(part)
    recs = recs[:count]
    wall = time.time() - t0
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    header = {"header": True, "n": n, "count": count, "seed_base": base,
              "seed_last": recs[-1]["seed"], "menu": MENU_NAMES, "wall_s": round(wall, 1),
              "cpu_s": round(sum(r["t"] for r in recs), 1), "workers": workers}
    with gzip.open(out, "wt") as fh:
        fh.write(json.dumps(header) + "\n")
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    return header


CURVE_T_MAX = 16384


def _work_curve(args) -> List[Dict]:
    """Supplementary positive-control frequency curve: CF signals only, T up to CURVE_T_MAX."""
    n, s0, s1 = args
    out = []
    for seed in range(s0, s1):
        g = instance(n, seed)
        if g is None:
            continue
        f = cf_features(g["N"], T_max=CURVE_T_MAX, T_list=(64,))      # operational: N only
        lab = labels(g["p"], g["q"])
        out.append({"seed": seed, "N": g["N"], "eps": lab["eps"], "qr2": lab["qr2"],
                    "first": {k[6:]: v for k, v in f.items() if k.startswith("first_")},
                    "cf_nonunit": f["cf_nonunit"], "cf_fb_ok": f["cf_fb_ok"]})
    return out


def run_curve(n: int, count: int, base: int, out: str, workers: int, chunk: int = 200) -> Dict:
    t0 = time.time()
    recs: List[Dict] = []
    next_chunk = 0
    with Pool(workers) as pool:
        while len(recs) < count:
            nchunks = max(workers, int((count - len(recs)) * 2.1 / chunk) + 1)
            jobs = [(n, base + (next_chunk + i) * chunk, base + (next_chunk + i + 1) * chunk)
                    for i in range(nchunks)]
            next_chunk += nchunks
            for part in pool.imap(_work_curve, jobs):
                recs.extend(part)
    recs = recs[:count]
    header = {"header": True, "kind": "pc_curve", "n": n, "count": count, "seed_base": base,
              "seed_last": recs[-1]["seed"], "T_max": CURVE_T_MAX, "wall_s": round(time.time() - t0, 1)}
    with gzip.open(out, "wt") as fh:
        fh.write(json.dumps(header) + "\n")
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    return header


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", required=True)
    ap.add_argument("--count", type=int, default=100000)
    ap.add_argument("--range", choices=sorted(SEED_BASE), required=True)
    ap.add_argument("--offset", type=int, default=0, help="added to the range base")
    ap.add_argument("--workers", type=int, default=30)
    ap.add_argument("--curve", action="store_true", help="supplementary PC frequency curve (CF only)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    for n in a.n:
        if a.curve:
            out = (a.out or "experiments/raw/eps_pccurve_n{n}.jsonl.gz").format(range=a.range, n=n)
            h = run_curve(n, a.count, SEED_BASE[a.range] + a.offset + n * 10 ** 9, out, a.workers)
        else:
            out = (a.out or "experiments/raw/eps_{range}_n{n}.jsonl.gz").format(range=a.range, n=n)
            h = run(n, a.count, SEED_BASE[a.range] + a.offset, out, a.workers)
        print(json.dumps({k: v for k, v in h.items() if k != "menu"}), flush=True)


if __name__ == "__main__":
    main()
