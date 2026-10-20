"""E1 harness (notes/nested_prereg.md).  HARNESS: uses only N from each instance (labels are never read).

Usage: .venv/bin/python -m src.nested.harness_e1 experiments/configs/nested_e1.json [--calib]
"""
from __future__ import annotations

import json
import random
import sys
import time

from src.common.instances import balanced_semiprime
from src.thmA import arith
from src.nested import toy_e1 as T


def run_instance(N: int, seed: int, cfg: dict) -> dict:
    t0 = time.time()
    B, Y, S = cfg["B"], cfg["Y"], cfg["S"]
    S1 = cfg.get("S1", S)
    ctx = T.Ctx(B, Y)
    rngA = random.Random(f"nested_e1|fieldA|{seed}")
    rngA2 = random.Random(f"nested_e1|fieldA2|{seed}")
    rngB = random.Random(f"nested_e1|baseline|{seed}")
    A = T.sample_field(N, rngA, M_bits=cfg.get("M_bits"))
    A2 = T.sample_field(N, rngA2, M_bits=cfg.get("M_bits"))
    U, U2 = T.Uniform(A, ctx), T.Uniform(A2, ctx)
    cA, cA2 = {}, {}
    parents = []
    lo = cfg["ell_lo"]
    tried = 0
    for ell in range(lo, lo + cfg["ell_span"]):
        if len(parents) >= cfg["max_parents"]:
            break
        if not arith.is_prime(ell):
            continue
        tried += 1
        node = T.parent_node(A, ell, S1, ctx)
        if node is None:
            continue
        recs = []
        for ch in node["children"]:
            kind, p, _ = ch
            rec = {"kind": kind, "p": p, "d1": T.d1(A, ch, S, ctx, cA)}
            bt = U.draw(kind, p, rngB)
            rec["base_p"] = bt[1]
            rec["base_d1"] = T.d1(A, bt, S, ctx, cA)
            if kind == "r":
                t2 = T.rational_target(A2, p)
                rec["d1_A2"] = None if t2 is None else T.d1(A2, t2, S, ctx, cA2)
                bt2 = U2.draw("r", p, rngB)
                rec["base_d1_A2"] = T.d1(A2, bt2, S, ctx, cA2)
            recs.append(rec)
        parents.append({"ell": ell, "a": node["a"], "b": node["b"], "degenerate": node["degenerate"],
                        "children": recs})
    out = {"seed": seed, "N_bits": N.bit_length(), "field": {k: A[k] for k in ("g0", "g1", "w", "c", "f")},
           "field2": {k: A2[k] for k in ("g0", "g1", "w", "c", "f")}, "ell_tried": tried, "parents": parents}
    if seed in cfg.get("dead_seeds", []):
        rngD = random.Random(f"nested_e1|dead|{seed}")
        dead = {}
        for kind, pool in (("r", [T.rational_target(A, q) for q in U.rat]), ("a", U.alg)):
            sample = [rngD.choice(pool) for _ in range(cfg["dead_sample"])]
            for Sx in (S // 2, S, 2 * S):
                dead[f"{kind}_S{Sx}"] = sum(not T.d1(A, t, Sx, ctx, cA) for t in sample) / len(sample)
        out["dead"] = dead
    out["time"] = time.time() - t0
    return out


def main():
    cfg = json.load(open(sys.argv[1]))
    calib = "--calib" in sys.argv
    seeds = cfg["calib_seeds"] if calib else cfg["seeds"]
    if calib and len(sys.argv) > 3:
        cfg["S"] = int(sys.argv[3])
    path = cfg["calib_out"] if calib else cfg["out"]
    with open(path, "a" if calib else "w") as fh:
        for s in seeds:
            inst = balanced_semiprime(cfg["bits"], s)
            r = run_instance(inst.N, s, cfg)
            r["S"] = cfg["S"]
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            nch = sum(len(p["children"]) for p in r["parents"])
            print(f"seed {s} S {cfg['S']}: parents {len(r['parents'])}/{r['ell_tried']} children {nch} "
                  f"time {r['time']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
