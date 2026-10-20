"""E2 harness (notes/almostall13_prereg.md): one-step descendability, random N versus fixed N.

HARNESS/DIAGNOSTIC code.  It reuses src/nested/toy_e1.py by import (unchanged).  It uses only inst.N from each
instance and never reads factor labels.

Usage: .venv/bin/python -m src.almostall13.harness_e2 experiments/configs/almostall13_e2.json
"""
from __future__ import annotations

import json
import random
import sys
import time

from src.common.instances import balanced_semiprime
from src.thmA import arith
from src.thmA import numfield as nf
from src.nested import toy_e1 as T


def run_field(N: int, field_seed: str, cfg: dict, ctx: T.Ctx) -> dict:
    t0 = time.time()
    S = cfg["S"]
    rngA = random.Random(f"{field_seed}|field")
    rngU = random.Random(f"{field_seed}|uniform")
    rngB = random.Random(f"{field_seed}|matched")
    A = T.sample_field(N, rngA, M_bits=cfg["M_bits"])
    U = T.Uniform(A, ctx)
    cache: dict = {}
    nu = cfg["n_uniform"]
    ua = [T.d1(A, rngU.choice(U.alg), S, ctx, cache) for _ in range(nu)]
    ur = [T.d1(A, T.rational_target(A, rngU.choice(U.rat)), S, ctx, cache) for _ in range(nu)]
    parents = []
    tried = 0
    lo = cfg["ell_lo"]
    for ell in range(lo, lo + cfg["ell_span"]):
        if len(parents) >= cfg["max_parents"]:
            break
        if not arith.is_prime(ell):
            continue
        tried += 1
        node = T.parent_node(A, ell, S, ctx)
        if node is None:
            continue
        recs = []
        for ch in node["children"]:
            kind, p, _ = ch
            bt = U.draw(kind, p, rngB)
            recs.append({"kind": kind, "p": p, "d1": T.d1(A, ch, S, ctx, cache), "base_d1": T.d1(A, bt, S, ctx, cache)})
        parents.append({"ell": ell, "children": recs})
    roots50 = sum(len(nf.roots_mod(A["f"], r)) for r in arith.primes_upto(50) if A["f"][A["d"]] % r)
    return {"field_seed": field_seed, "N_bits": N.bit_length(),
            "field": {k: A[k] for k in ("g0", "g1", "w", "c", "f")},
            "u_alg": sum(ua), "u_rat": sum(ur), "n_uniform": nu, "pool_alg": len(U.alg), "pool_rat": len(U.rat),
            "roots50": roots50, "ell_tried": tried, "parents": parents, "time": time.time() - t0}


def main():
    cfg = json.load(open(sys.argv[1]))
    ctx = T.Ctx(cfg["B"], cfg["Y"])
    with open(cfg["out"], "w") as fh:
        jobs = []
        for s in cfg["fix_seeds"]:
            for j in range(cfg["fields_per_fixed"]):
                jobs.append(("FIX", s, f"almostall13|FIX|{s}|{j}"))
        for s in cfg["rnd_seeds"]:
            jobs.append(("RND", s, f"almostall13|RND|{s}"))
        for regime, s, fs in jobs:
            N = balanced_semiprime(cfg["bits"], s).N
            r = run_field(N, fs, cfg, ctx)
            r["regime"], r["N_seed"] = regime, s
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(f"{regime} N{s} {fs}: u_alg {r['u_alg']}/{r['n_uniform']} parents {len(r['parents'])} "
                  f"time {r['time']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
