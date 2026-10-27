"""HARNESS. E1 identity checks and E2 trichotomy decoder (notes/nativeinv_prereg.md)."""
from __future__ import annotations

import json
import os
import sys
import time

from src.common.instances import balanced_semiprime
from . import harness as H
from .trichotomy import decode, global_tier, pair_tier

RAW = "experiments/raw/nativeinv_e2.jsonl"


def e1(small=((3, 5), (3, 7), (3, 11), (5, 7))):
    rows = []
    for p, q in small:
        N = p * q
        r = {"N": N}
        r["kSL2"] = (H.k_SL2(N), (p + 4) * (q + 4), H.gcdsum_SL2(N))
        if N <= 21:
            r["kGL2"] = (H.k_GL2(N), (p * p - 1) * (q * q - 1), H.gcdsum_GL2(N))
        r["kAff"] = (H.k_Aff(N), N)
        degs = H.aff_degree_multiset(p, q)
        r["aff_deg_check"] = (sum(m * d * d for d, m in degs), N * (p - 1) * (q - 1), sum(m for _, m in degs), N)
        r["kHeis"] = (H.k_Heis(N), (p * p + p - 1) * (q * q + q - 1))
        r["pillai"] = (H.pillai(N), H.zero_product_count(N), (2 * p - 1) * (2 * q - 1))
        if N <= 15:
            r["subgroups"] = (H.subgroups_ZN2(N), (p + 3) * (q + 3))
        rows.append(r)
        print(r, flush=True)
    return rows


def e2(bits=(32, 64, 128), seeds=range(40)):
    out = []
    for n in bits:
        for s in seeds:
            inst = balanced_semiprime(n, s)
            N = inst.N
            p, q = inst.labels["p"], inst.labels["q"]
            for key, (F, tier) in H.REGISTRY.items():
                V = F.value(p) * F.value(q)          # harness computes the invariant from the factors
                t0 = time.perf_counter()
                res = decode(N, V, F)                # operational: sees N, V, F only
                dt = time.perf_counter() - t0
                true_pair = pair_tier(F, p % F.m if F.m > 1 else 0, q % F.m if F.m > 1 else 0)
                split_ok = res.split == (p, q)
                rec = {"n": n, "seed": s, "F": key, "declared": tier, "global_tier": global_tier(F),
                       "true_pair_tier": true_pair, "split": split_ok, "split_any": res.split is not None,
                       "candidates": res.candidates_tested, "sec": round(dt, 4)}
                if tier == "AGG":
                    import math
                    chi = 1 if N % 4 == 1 else -1
                    rec["agg_formula_ok"] = (V == chi * N * N)
                out.append(rec)
    os.makedirs(os.path.dirname(RAW), exist_ok=True)
    with open(RAW, "w") as fh:
        for r in out:
            fh.write(json.dumps(r) + "\n")
    # summary
    summ = {}
    for r in out:
        k = r["F"]
        d = summ.setdefault(k, {"declared": r["declared"], "global": r["global_tier"], "N": 0, "split": 0,
                                "fact_pairs": 0, "split_on_fact_pairs": 0, "agg_ok": 0})
        d["N"] += 1
        d["split"] += r["split"]
        d["fact_pairs"] += r["true_pair_tier"] == "FACT"
        d["split_on_fact_pairs"] += r["split"] and r["true_pair_tier"] == "FACT"
        d["agg_ok"] += r.get("agg_formula_ok", False)
    for k, d in summ.items():
        print(k, d)
    return summ


if __name__ == "__main__":
    which = sys.argv[1:] or ["e1", "e2"]
    if "e1" in which:
        e1()
    if "e2" in which:
        e2()
