"""HARNESS runner for E2 (unit-twist census).  Usage: run_e2 cfg.json out.jsonl [nproc]"""
from __future__ import annotations
import json
import sys
from multiprocessing import Pool

from ..nfsfruit.instances import composite
from .unit_twist import cf_period, classify_unit


def _job(a):
    cls, bits, seed, cap = a
    N, fs = composite(cls, bits, seed)
    r = cf_period(N, cap)
    rec = dict(cls=cls, bits=bits, seed=seed, N=str(N))
    if r is None:
        rec["capped"] = True
        return rec
    rec.update(P=r["P"], sign=r["sign"], ok_Q=r["ok_Q"], ok_A=r["ok_A"], ok_unit=r["ok_unit"],
               unit=classify_unit(N, r["a"], r["sign"], fs))
    return rec


def main():
    cfg = json.load(open(sys.argv[1]))
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    jobs = [(c, b, s, cfg["cap"]) for b in cfg["bits"] for c in cfg["classes"]
            for s in range(cfg["seeds"])]
    with Pool(nproc) as pool, open(sys.argv[2], "w") as fh:
        for rec in pool.imap_unordered(_job, jobs):
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
