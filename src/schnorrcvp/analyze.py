"""Analysis of E1/E2 raw data: the pre-registered statistics, plus clearly labelled POST-HOC models.

Usage: .venv/bin/python -m src.schnorrcvp.analyze e1|e2|posthoc_e1|posthoc_e2
"""
from __future__ import annotations

import json
import math
import statistics as st
import sys
from collections import defaultdict


def e1():
    rows = [json.loads(l) for l in open("experiments/raw/schnorrcvp_E1.jsonl")]
    cells = defaultdict(lambda: [0, 0.0, 0.0, 0, 0])
    for r in rows:
        c = cells[(r["bits"], r["n"])]
        c[0] += r["obs_primitive"]; c[1] += r["pred_primitive"]; c[2] += r["pred_naive"]; c[3] += r["obs_all"]; c[4] += 1
    tot_o = sum(c[0] for c in cells.values()); tot_p = sum(c[1] for c in cells.values())
    print("cell bits n | seeds | obs_all | obs_prim | pred_H_S' | ratio | pred_naive | naive>=obs")
    p1_ok = True; p2_ok = True
    for k in sorted(cells):
        o, p, nv, a, s = cells[k]
        ratio = o / p if p else float("nan")
        if p >= 10 and not (1 / 3 <= ratio <= 3):
            p1_ok = False
        if nv < o:
            p2_ok = False
        print(k, s, a, o, round(p, 2), round(ratio, 3), round(nv, 1), nv >= o)
    R = tot_o / tot_p
    print("pooled obs/pred H_S' =", round(R, 3), " E1-P1 pooled in [1/2,2]:", 0.5 <= R <= 2, " per-cell:", p1_ok)
    print("E1-P2 naive >= obs_primitive in every cell:", p2_ok)


def e2(path="experiments/raw/schnorrcvp_E2.jsonl"):
    rows = [json.loads(l) for l in open(path)]
    agg = defaultdict(lambda: [0, 0.0, 0.0, 0])  # rel, sum p_ref, sum p_naive, outputs
    by_cell = defaultdict(lambda: [0, 0])
    s21 = defaultdict(list)
    list0 = defaultdict(lambda: [0, set()])
    for r in rows:
        for o in r["outputs"]:
            if o["p_ref"] is None:
                continue
            a = agg[r["variant"]]
            a[0] += o["rel"]; a[1] += o["p_ref"]; a[2] += o["p_naive"]; a[3] += 1
            if o["kind"].startswith("list"):
                bc = by_cell[(r["n"], r["variant"], r["bits"])]
                bc[0] += o["rel"]; bc[1] += 1
                if r["variant"].startswith("S21") and r["bits"] >= 32 and o["log2w"] is not None:
                    s21[r["variant"]].append(o["log2w"] / r["bits"])
            if o["kind"] == "list0" and o["rel"]:
                l0 = list0[(r["bits"], r["seed"], r["n"], r["variant"])]
                l0[0] += 1; l0[1].add((o["u"], o["w"]))
    tr = sum(a[0] for a in agg.values()); tp = sum(a[1] for a in agg.values()); tn = sum(a[2] for a in agg.values())
    print("variant | outputs | rel | sum p_ref | ratio | sum p_naive | ratio_naive")
    p1 = True
    for v, a in sorted(agg.items()):
        rr = a[0] / a[1] if a[1] else float("nan")
        if a[1] >= 10 and not (1 / 3 <= rr <= 3):
            p1 = False
        print(v, a[3], a[0], round(a[1], 2), round(rr, 3), round(a[2], 2), round(a[0] / a[2], 3) if a[2] else None)
    R = tr / tp
    print("E2-P1 pooled ratio", round(R, 3), "in [1/2,2]:", 0.5 <= R <= 2, "per-variant:", p1)
    print("E2-P2 refined better:", abs(math.log(tr / tp)) < abs(math.log(tr / tn)), round(tr / tn, 3))
    p3 = True; viol = []
    for (n, var, bits), (rel, cnt) in by_cell.items():
        if bits == 20 and rel >= 1:
            hi = by_cell.get((n, var, 40))
            if hi and hi[0] / hi[1] >= rel / cnt:
                p3 = False; viol.append((n, var, rel / cnt, hi[0] / hi[1]))
    print("E2-P3 yield(40) < yield(20) in every eligible cell:", p3, viol)
    for v, xs in s21.items():
        print("E2-P4", v, "median log2|w|/bits (bits>=32):", round(st.median(xs), 3), ">= 0.5:", st.median(xs) >= 0.5)
    # yield table per (n, variant, bits)
    print("yield per list output: n variant | bits:yield")
    keys = sorted(set((n, v) for (n, v, b) in by_cell))
    for n, v in keys:
        print(n, v, " ".join("%d:%.3f" % (b, by_cell[(n, v, b)][0] / by_cell[(n, v, b)][1])
                            for b in sorted(set(b for (nn, vv, b) in by_cell if nn == n and vv == v))))
    print("E2-P5 duplication (list0 relations: count vs distinct), cells with >=3:")
    dup = [(k, c, len(s)) for k, (c, s) in list0.items() if c >= 3]
    for k, c, d in sorted(dup)[:40]:
        print(" ", k, c, d)
    if dup:
        print("  total list0 relations", sum(c for _, c, _ in dup), "distinct", sum(d for _, _, d in dup))
    try:
        for l in open(path.replace("E2.jsonl", "E2_factor.jsonl")):
            f = json.loads(l)
            print("E2-P6", f["bits"], f["seed"], f["distinct_relations"], f["fb_size"], f["factor"], f["divides"])
    except FileNotFoundError:
        pass


def posthoc_e1():
    """POST-HOC: H_S'' (full CRT local factors, parity) on the E1 boxes."""
    from .core import first_primes
    from .psi import expected_primitive_local
    rows = [json.loads(l) for l in open("experiments/raw/schnorrcvp_E1.jsonl")]
    cells = defaultdict(lambda: [0, 0.0])
    out = open("experiments/raw/schnorrcvp_E1_posthoc.jsonl", "w")
    for r in rows:
        e = expected_primitive_local(int(r["N"]), first_primes(r["n"]), r["V"], r["W"], K=2000, h=0.04)["E"]
        out.write(json.dumps({"bits": r["bits"], "n": r["n"], "seed": r["seed"], "pred_local": e,
                              "obs_primitive": r["obs_primitive"]}) + "\n")
        c = cells[(r["bits"], r["n"])]
        c[0] += r["obs_primitive"]; c[1] += e
    for k in sorted(cells):
        o, p = cells[k]
        print(k, o, round(p, 2), round(o / p, 3))
    to = sum(c[0] for c in cells.values()); tp = sum(c[1] for c in cells.values())
    print("POST-HOC pooled obs/pred H_S'' =", round(to / tp, 3))


def posthoc_e2():
    """POST-HOC: H_S'' local factors (psi.local_factor) applied to E2 outputs (w the unknown member)."""
    import math as _m
    from .core import first_primes
    from .psi import local_factor
    from .run_e2 import Res, supp
    rows = [json.loads(l) for l in open("experiments/raw/schnorrcvp_E2.jsonl")]
    resc = {}
    agg = defaultdict(lambda: [0, 0.0])
    for r in rows:
        key = (r["N"], r["n"])
        if key not in resc:
            resc[key] = Res(first_primes(r["n"]), 6 * _m.log(int(r["N"])) + 40)
        res = resc[key]
        P = res.P
        for o in r["outputs"]:
            if o["p_ref"] is None or o["u"] is None:
                continue
            u, v, w = int(o["u"]), int(o["v"]), int(o["w"])
            lx = _m.log(abs(w))
            p = res.d(lx) * local_factor(res.d, lx, P, supp(u * v, P))
            a = agg[r["variant"]]
            a[0] += o["rel"]; a[1] += p
    for k, (o, p) in sorted(agg.items()):
        print("POST-HOC", k, o, round(p, 1), round(o / p, 3))
    to = sum(a[0] for a in agg.values()); tp = sum(a[1] for a in agg.values())
    print("POST-HOC pooled", round(to / tp, 3))


def corrected_e2(W0: int = 200000, h: float = 0.1):
    """POST-HOC correction of the pre-registered H_res' implementation (codex schnorrcvp_2 #18: 2066 predictions
    exceeded 1 and 66 relations had p = 0).  For |w| <= W0: exact discrete probability
        #{m in [|w|e^{-h/2}, |w|e^{h/2}] : gcd(m, Q) = 1, m smooth} / #{m in window : gcd(m, Q) = 1},  Q = supp(uv);
    for |w| <= p_n it is exactly 1 (every integer <= B is B-smooth).  For |w| > W0: the pre-registered value clipped to
    [0, 1].  Also reports the plain clipping diagnostic."""
    import numpy as np
    from .core import first_primes
    from .enum_rel import smooth_mask_window
    from .run_e2 import supp
    rows = [json.loads(l) for l in open("experiments/raw/schnorrcvp_E2.jsonl")]
    masks = {}
    agg = defaultdict(lambda: [0, 0.0, 0.0, 0])
    over1 = zero_rel = 0
    for r in rows:
        n = r["n"]
        P = first_primes(n)
        if n not in masks:
            masks[n] = np.concatenate([[False], smooth_mask_window(1, int(W0 * 1.2) + 2, P)])
        mask = masks[n]
        for o in r["outputs"]:
            if o["p_ref"] is None:
                continue
            over1 += o["p_ref"] > 1
            zero_rel += (o["p_ref"] == 0 and o["rel"])
            a = agg[r["variant"]]
            a[0] += o["rel"]; a[2] += min(1.0, o["p_ref"]); a[3] += 1
            if o["w"] is None:
                a[1] += min(1.0, o["p_ref"]); continue
            w = abs(int(o["w"]))
            if w <= P[-1]:
                p = 1.0
            elif w <= W0:
                lo = max(1, int(w * math.exp(-h / 2))); hi = int(w * math.exp(h / 2)) + 1
                m = np.arange(lo, hi + 1)
                cop = np.ones(len(m), dtype=bool)
                for q in supp(int(o["u"]) * int(o["v"]), P):
                    cop &= (m % q != 0)
                denom = int(cop.sum())
                p = float((cop & mask[lo:hi + 1]).sum()) / denom if denom else 0.0
            else:
                p = min(1.0, o["p_ref"])
            a[1] += p
    print("pre-registered p_ref > 1:", over1, " relations with p_ref = 0:", zero_rel)
    for k, (o, p, pc, cnt) in sorted(agg.items()):
        print("CORRECTED", k, o, round(p, 1), round(o / p, 3), " clipped-only", round(o / pc, 3))
    to = sum(a[0] for a in agg.values()); tp = sum(a[1] for a in agg.values()); tc = sum(a[2] for a in agg.values())
    print("CORRECTED pooled", round(to / tp, 3), " clipped-only pooled", round(to / tc, 3))


if __name__ == "__main__":
    {"e1": e1, "e2": e2, "posthoc_e1": posthoc_e1, "posthoc_e2": posthoc_e2, "corrected_e2": corrected_e2}[sys.argv[1]]()
