"""Experiment drivers for pm1bsgs (HARNESS side: uses factors for labels and
exact outcome prediction; the predictors are validated against operational runs
by `validate`).

Usage (from project root):
  .venv/bin/python -m src.pm1bsgs.experiments smooth   --count 4000
  .venv/bin/python -m src.pm1bsgs.experiments phase1   --count 600
  .venv/bin/python -m src.pm1bsgs.experiments validate --count 40
  .venv/bin/python -m src.pm1bsgs.experiments rho      --count 300
  .venv/bin/python -m src.pm1bsgs.experiments walltime
  .venv/bin/python -m src.pm1bsgs.analysis            (budget curves, restarts, fits)
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import random
import time
from multiprocessing import Pool
from typing import Dict, List

from src.common.baselines import pollard_rho_brent
from src.common.instances import balanced_semiprime
from src.common.metrics import Metrics

from . import costmodel as cm
from . import harness as H
from . import operational as op

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
CFG = os.path.join(ROOT, "experiments", "configs")

NS = (32, 40, 48, 56, 64, 72, 80)
NBASES = 32
B1_GRID = [1] + [1 << k for k in range(3, 31)]
BDAY_B1 = [1] + [1 << k for k in range(4, 29, 2)]
BDAY_KMAX = 1 << 21


def _w(path: str, recs: List[Dict]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    op_ = gzip.open if path.endswith(".gz") else open
    with op_(path, "wt") as f:
        for r in recs:
            f.write(json.dumps(r, sort_keys=True) + "\n")


def _cfg(name: str, d: Dict) -> None:
    os.makedirs(CFG, exist_ok=True)
    with open(os.path.join(CFG, f"pm1bsgs_{name}.json"), "w") as f:
        json.dump(d, f, indent=1, sort_keys=True)


# ------------------------------------------------------------------ E1 smooth parts
def _smooth_worker(args):
    n, seed = args
    inst = balanced_semiprime(n, seed)
    out = []
    for which, r in zip(("p", "q"), inst.labels["factors"]):
        P = H.PrimeInfo.of(r)
        out.append({"n": n, "seed": seed, "which": which, "bits": r.bit_length(),
                    "fac": sorted([[int(l), int(v)] for l, v in P.fac.items()])})
    return out


def run_smooth(count: int) -> None:
    jobs = [(n, s) for n in NS for s in range(count)]
    with Pool(32) as pool:
        res = pool.map(_smooth_worker, jobs, chunksize=64)
    _w(os.path.join(RAW, "pm1bsgs_smooth.jsonl.gz"), [r for rr in res for r in rr])
    _cfg("smooth", {"family": "balanced_semiprime", "ns": NS, "seeds": [0, count], "note":
                    "factorisation of p-1 and q-1 (harness), for smooth-part statistics"})


# ------------------------------------------------------------------ phase 1
def _phase1_worker(args):
    n, seed = args
    t0 = time.time()
    inst = balanced_semiprime(n, seed)
    N = inst.N
    p, q = inst.labels["factors"]
    P, Q = H.PrimeInfo.of(p), H.PrimeInfo.of(q)
    bases = []
    for bseed in range(NBASES):
        bs = H.base_state(N, P, Q, bseed)
        bases.append([bs.op, bs.oq])
    bs0 = H.base_state(N, P, Q, 0)
    bday = {}
    for B1 in BDAY_B1:
        st = H.stage1_outcome(bs0, B1)
        if st[0] != "ok":
            bday[str(B1)] = st[0]
        else:
            fp = H.birthday_first_pure(N, 0, st[1], st[2], BDAY_KMAX)
            bday[str(B1)] = fp if fp is not None else "none"
    return {"n": n, "seed": seed, "N": str(N), "p": str(p), "q": str(q),
            "pfac": sorted([[int(l), int(v)] for l, v in P.fac.items()]),
            "qfac": sorted([[int(l), int(v)] for l, v in Q.fac.items()]),
            "bases": [[str(a), str(b)] for a, b in bases], "bday_first_pure_base0": bday,
            "wall_s": round(time.time() - t0, 3)}


def run_phase1(count: int) -> None:
    jobs = [(n, s) for n in NS for s in range(count)]
    random.Random(0).shuffle(jobs)
    with Pool(32) as pool:
        res = pool.map(_phase1_worker, jobs, chunksize=4)
    res.sort(key=lambda r: (r["n"], r["seed"]))
    _w(os.path.join(RAW, "pm1bsgs_phase1.jsonl.gz"), res)
    _cfg("phase1", {"family": "balanced_semiprime", "ns": NS, "seeds": [0, count],
                    "train_seeds": [0, count // 2], "holdout_seeds": [count // 2, count],
                    "nbases": NBASES, "B1_grid": B1_GRID, "bday_B1": BDAY_B1, "bday_Kmax": BDAY_KMAX,
                    "base_choice": "operational.choose_base(N, base_seed)"})


# ------------------------------------------------------------------ validation (operational vs predictor)
def _validate_worker(args):
    n, seed, rs = args
    rng = random.Random(f"validate|{n}|{seed}|{rs}")
    inst = balanced_semiprime(n, seed)
    N = inst.N
    p, q = inst.labels["factors"]
    P, Q = H.PrimeInfo.of(p), H.PrimeInfo.of(q)
    bseed = rng.randrange(NBASES)
    bs = H.base_state(N, P, Q, bseed)
    out = []
    smax = min(17, n // 4 + 2)
    # (a)/(d)
    for _ in range(3):
        B1 = rng.choice([1] + [1 << k for k in range(3, 17)])
        s = rng.randrange(1, smax)
        mb, J = 1 << s, (1 << s) * rng.choice([1, 2])
        mode = rng.choice(["N", "order"])
        met = Metrics()
        r = op.pm1_bsgs(N, bseed, B1, mb, J, mode, met)
        lab = H.outcome_from_orders(bs.op, bs.opf, bs.oq, bs.oqf, N, "bsgs",
                                    {"B1": B1, "m": mb, "J": J, "mode": mode})
        out.append(_vrec("bsgs", n, seed, bseed, {"B1": B1, "m": mb, "J": J, "mode": mode}, r, lab, met))
    # (b)
    B1 = rng.choice([1 << k for k in range(3, 17)])
    B2 = B1 << rng.randrange(1, 8)
    met = Metrics()
    r = op.pm1_prime(N, bseed, B1, B2, met)
    lab = H.outcome_from_orders(bs.op, bs.opf, bs.oq, bs.oqf, N, "prime", {"B1": B1, "B2": B2})
    out.append(_vrec("prime", n, seed, bseed, {"B1": B1, "B2": B2}, r, lab, met))
    # (c)
    B1 = rng.choice(BDAY_B1[:8])
    K = 1 << rng.randrange(2, smax + 1)
    met = Metrics()
    r = op.pm1_birthday(N, bseed, B1, K, met)
    st = H.stage1_outcome(bs, B1)
    if st[0] == "ok":
        fp = H.birthday_first_pure(N, bseed, st[1], st[2], K)
        lab = "success2" if fp is not None else "fail2"
    else:
        lab = st[0]
    out.append(_vrec("bday", n, seed, bseed, {"B1": B1, "K": K}, r, lab, met))
    return out


def _vrec(method, n, seed, bseed, params, r, lab, met):
    ok = (r["factor"] is not None) == (lab in ("stage1", "success2"))
    rec = met.to_record(method=method, n=n, seed=seed, base_seed=bseed, params=params,
                        found=r["factor"] is not None, stage=r["stage"], cls=r["cls"],
                        predicted=lab, agree=ok)
    return rec


def run_validate(count: int) -> None:
    jobs = [(n, s, k) for n in NS for s in range(count) for k in range(2)]
    with Pool(32) as pool:
        res = pool.map(_validate_worker, jobs, chunksize=1)
    recs = [r for rr in res for r in rr]
    _w(os.path.join(RAW, "pm1bsgs_validate.jsonl"), recs)
    bad = [r for r in recs if not r["agree"]]
    print(f"validate: {len(recs)} runs, {len(bad)} disagreements")
    _cfg("validate", {"ns": NS, "seeds": [0, count], "runs_per_seed": 2})


# ------------------------------------------------------------------ rho baseline
def _rho_worker(args):
    n, seed = args
    inst = balanced_semiprime(n, seed)
    met = Metrics()
    d = pollard_rho_brent(inst.N, seed=seed, m=met)
    c = met.counts
    return {"n": n, "seed": seed, "found": d is not None and 1 < d < inst.N,
            "rn_mul": c.get("rn_mul", 0), "gcd": c.get("gcd", 0), "rn_add": c.get("rn_add", 0),
            "eq_ops": c.get("rn_mul", 0) + c.get("gcd", 0), "wall_s": round(met.elapsed(), 4),
            "p_bits": inst.labels["factors"][0].bit_length()}


def run_rho(count: int, start: int) -> None:
    jobs = [(n, s) for n in NS for s in range(start, start + count)]
    with Pool(32) as pool:
        res = pool.map(_rho_worker, jobs, chunksize=2)
    _w(os.path.join(RAW, "pm1bsgs_rho.jsonl"), res)
    _cfg("rho", {"ns": NS, "seeds": [start, start + count], "impl": "src.common.baselines.pollard_rho_brent"})


# ------------------------------------------------------------------ wall-time calibration
def run_walltime() -> None:
    recs = []
    inst = balanced_semiprime(80, 12345)
    N = inst.N
    for s in range(8, 19, 2):
        mb = 1 << s
        met = Metrics()
        t0 = time.perf_counter()
        b = op.choose_base(N, 1)
        op.bsgs_stage2(N, pow(b, 1 << 20, N), mb, mb, "order", met)
        dt = time.perf_counter() - t0
        eq = met.counts.get("rn_mul", 0) + met.counts.get("gcd", 0) + met.counts.get("poly_eq", 0)
        recs.append({"what": "bsgs_stage2", "m": mb, "J": mb, "eq_ops": eq, "wall_s": dt,
                     "ns_per_eq": 1e9 * dt / eq, "counts": met.counts})
        print(recs[-1]["what"], mb, eq, round(dt, 3), round(1e9 * dt / eq, 1))
    for K in (1 << 10, 1 << 14, 1 << 18):
        met = Metrics()
        t0 = time.perf_counter()
        x = 3
        for _ in range(K):
            x = x * x % N
        dt = time.perf_counter() - t0
        recs.append({"what": "scalar_modmul_loop", "K": K, "eq_ops": K, "wall_s": dt, "ns_per_eq": 1e9 * dt / K})
        print("scalar", K, round(1e9 * dt / K, 1))
    _w(os.path.join(RAW, "pm1bsgs_walltime.jsonl"), recs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    ap.add_argument("--count", type=int, default=100)
    ap.add_argument("--start", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "smooth":
        run_smooth(a.count)
    elif a.cmd == "phase1":
        run_phase1(a.count)
    elif a.cmd == "validate":
        run_validate(a.count)
    elif a.cmd == "rho":
        run_rho(a.count, a.start)
    elif a.cmd == "walltime":
        run_walltime()
    else:
        raise SystemExit(a.cmd)


if __name__ == "__main__":
    main()
