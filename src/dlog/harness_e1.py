"""HARNESS / DIAGNOSTIC for E1 (notes/dlog_prereg.md).  Generates instances (p, g primitive root, h = g^X) with
sympy (factor-aware: primitive_root, totient) and runs the OPERATIONAL src.dlog.pipeline.dlog, which never sees
X or the factorisation of p - 1.  Writes JSONL."""
from __future__ import annotations

import json
import random
import sys
import time

import sympy

from src.dlog.pipeline import dlog


def make_instance(bits: int, seed: int):
    rng = random.Random(f"dlogE1|{bits}|{seed}")
    while True:
        p = rng.randrange(2 ** (bits - 1) + 1, 2 ** bits)
        if sympy.isprime(p):
            break
    g = sympy.primitive_root(p)
    X = rng.randrange(p - 1)
    return p, int(g), pow(int(g), X, p), X


def main(cfg_path: str, out_path: str):
    cfg = json.load(open(cfg_path))
    with open(out_path, "w") as fo:
        for blk in cfg["blocks"]:
            for inst in blk["instances"]:
                p, g, h, X = make_instance(inst["bits"], inst["seed"])
                fac = {int(k): int(v) for k, v in sympy.factorint(p - 1).items()}   # diagnostic only
                t0 = time.time()
                pl = dict(cfg["pipeline"]); pl.update(blk.get("pipeline_override", {}))
                try:
                    res = dlog(p, g, h, pl, inst["seed"])
                    err = None
                except Exception as e:  # recorded, never hidden
                    res, err = {"X": None, "ctr": {}, "deps": []}, repr(e)
                rec = {"block": blk["name"], "inst": inst, "p": p, "g": g, "h": h, "X_true": X,
                       "p1_factorisation": fac, "phi_ratio_num": int(sympy.totient(p - 1)), "X_out": res["X"],
                       "correct": res["X"] is not None and pow(g, res["X"], p) == h,
                       "X_equal": res["X"] is not None and (res["X"] - X) % (p - 1) == 0,
                       "error": err, "ctr": res["ctr"], "deps": res["deps"], "wall": round(time.time() - t0, 2)}
                fo.write(json.dumps(rec, default=str) + "\n")
                fo.flush()
                print(blk["name"], inst, rec["correct"], rec["ctr"].get("rows_last"), rec["ctr"].get("cols"),
                      rec["ctr"].get("kernel_dim"), rec["ctr"].get("unit_relations"), rec["ctr"].get("r1r2"),
                      rec["wall"], err, flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
