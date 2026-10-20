"""Cost-model validation experiments V1-V3 (pre-registered in notes/cadobench_prereg.md).

V1: single-side smoothness of the actual norms at random points of real special-q lattices vs Dickman
    rho at the actual u (with / without Murphy alpha).
V2: exact joint (both sides) smoothness on exhaustive windows vs the product of the marginals.
V3: las recall on the same windows (las run on that special-q with the production parameters).

All smoothness decisions are exact (smoothcount). Floats only for rho / alpha / summaries.
Diagnostic harness: it reads N only through CADO's polynomial files; nothing here touches p, q.

Usage:
  python -m cadobench.validate v1 --digits 90 --out experiments/raw/cadobench_v1_c90.json
  python -m cadobench.validate v23 --digits 70 --out experiments/raw/cadobench_v23_c70.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import random
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

from .lasbench import Testbed, run_las
from .relstats import qlattices
from .smoothcount import (QLattice, SideParams, Smoother, murphy_alpha, rho, sample_points,
                          window_points)

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")
WORKDIR_RAW = os.path.expanduser("~/cadobench-work/lasraw")
os.makedirs(WORKDIR_RAW, exist_ok=True)


def fast_hom(coeffs: Sequence[int], a: int, b: int) -> int:
    d = len(coeffs) - 1
    acc = coeffs[d]
    bp = 1
    for k in range(d - 1, -1, -1):
        bp *= b
        acc = acc * a + coeffs[k] * bp
    return acc


def setup(digits: int, tag: str = "run1"):
    from .run_cado import WORK
    if not os.path.isdir(os.path.join(WORK, f"c{digits}_{tag}")) and digits == 60:
        tag = "smoke"  # the c60 baseline is the smoke-test factorization
    tb = Testbed.from_workdir(f"c{digits}_{tag}")
    P = tb.read_poly()
    p = tb.params
    side = [SideParams(int(p["lim0"]), int(p["lpb0"]), int(p["mfb0"])),
            SideParams(int(p["lim1"]), int(p["lpb1"]), int(p["mfb1"]))]
    return tb, P, side


def classify_points(points, f0, f1, q, sm: Smoother):
    n0 = [abs(fast_hom(f0, a, b)) for (_, _, a, b) in points]
    n1 = []
    for (_, _, a, b) in points:
        v = abs(fast_hom(f1, a, b))
        assert v % q == 0
        n1.append(v // q)
    c0 = sm.classify(n0, 0)
    c1 = sm.classify(n1, 1)
    return n0, n1, c0, c1


def _u(n: int, lpb: int, alpha: float = 0.0) -> float:
    if n <= 1:
        return 0.0
    return (math.log(n) + alpha) / (lpb * math.log(2))


def v1(digits: int, nq: int = 20, per_q: int = 10000, seed: int = 0) -> dict:
    tb, P, side = setup(digits)
    f0, f1 = P["f0"], P["f1"]
    alpha = [murphy_alpha(f0), murphy_alpha(f1)]
    sm = Smoother(side)
    lats = qlattices(tb.workdir)
    pick = [lats[k * (len(lats) - 1) // max(nq - 1, 1)] for k in range(nq)]
    rows = []
    agg = {k: 0.0 for k in ("pts", "s0", "s1", "c0", "c1", "sB", "sC", "r0", "r1", "r0a", "r1a",
                              "rj", "rja", "u0", "u1")}
    bins: Dict[str, Dict[str, float]] = {}
    for lat, nrel in pick:
        rng = random.Random(f"cadobench-v1|{digits}|{lat.q}|{lat.rho}|{seed}")
        pts = sample_points(lat, per_q, rng)
        n0, n1, c0, c1 = classify_points(pts, f0, f1, lat.q, sm)
        r = {"q": lat.q, "rho": lat.rho, "I": lat.I, "J": lat.J, "las_rels": nrel, "pts": len(pts)}
        r["s0"] = sum(c0["smooth"]); r["s1"] = sum(c1["smooth"])
        r["c0"] = sum(c0["cado"]); r["c1"] = sum(c1["cado"])
        r["sB"] = sum(x and y for x, y in zip(c0["smooth"], c1["smooth"]))
        r["sC"] = sum(x and y for x, y in zip(c0["cado"], c1["cado"]))
        u0 = [_u(x, side[0].lpb) for x in n0]
        u1 = [_u(x, side[1].lpb) for x in n1]
        u0a = [_u(x, side[0].lpb, alpha[0]) for x in n0]
        u1a = [_u(x, side[1].lpb, alpha[1]) for x in n1]
        r["r0"] = sum(map(rho, u0)); r["r1"] = sum(map(rho, u1))
        r["r0a"] = sum(map(rho, u0a)); r["r1a"] = sum(map(rho, u1a))
        r["rj"] = sum(rho(x) * rho(y) for x, y in zip(u0, u1))
        r["rja"] = sum(rho(x) * rho(y) for x, y in zip(u0a, u1a))
        r["u0"] = sum(u0); r["u1"] = sum(u1)
        for k in agg:
            agg[k] += r[k]
        # binned (alpha-corrected u, side) measured vs predicted
        for s, ua, flags in ((0, u0a, c0["smooth"]), (1, u1a, c1["smooth"])):
            for x, f in zip(ua, flags):
                key = f"{s}:{math.floor(x * 4) / 4:.2f}"
                b = bins.setdefault(key, {"n": 0, "smooth": 0, "rho_sum": 0.0})
                b["n"] += 1; b["smooth"] += int(f); b["rho_sum"] += rho(x)
        # las joint fraction over the whole region: coprime points ~ (6/pi^2) * I * J
        r["las_frac"] = nrel / (6 / math.pi ** 2 * lat.I * lat.J)
        rows.append(r)
        print(f"[v1 c{digits}] q={lat.q} s0={r['s0']} s1={r['s1']} sB={r['sB']}", file=sys.stderr, flush=True)
    n = agg["pts"]
    summ = {
        "points": int(n),
        "alpha": alpha,
        "mean_u": [agg["u0"] / n, agg["u1"] / n],
        "frac_smooth": [agg["s0"] / n, agg["s1"] / n],
        "frac_cado_admissible": [agg["c0"] / n, agg["c1"] / n],
        "mean_rho": [agg["r0"] / n, agg["r1"] / n],
        "mean_rho_alpha": [agg["r0a"] / n, agg["r1a"] / n],
        "ratio_meas_over_rho": [agg["s0"] / agg["r0"], agg["s1"] / agg["r1"]],
        "ratio_meas_over_rho_alpha": [agg["s0"] / agg["r0a"], agg["s1"] / agg["r1a"]],
        "joint_SB_sample": int(agg["sB"]), "joint_SCADO_sample": int(agg["sC"]),
        "pred_joint_rho": agg["rj"] / n, "pred_joint_rho_alpha": agg["rja"] / n,
        "las_joint_frac_mean": sum(r["las_frac"] for r in rows) / len(rows),
    }
    return {"digits": digits, "params": tb.params, "f0": f0, "f1": f1, "nq": nq, "per_q": per_q, "seed": seed,
            "summary": summ, "per_q_rows": rows, "bins": dict(sorted(bins.items()))}


WINDOWS = {60: None, 70: (-512, 512, 0, 1024), 80: (-512, 512, 0, 1024)}


def las_relations_for_q(tb: Testbed, lat: QLattice, threads: int = 2) -> Dict[Tuple[int, int], str]:
    """Relations las (production params) outputs for exactly this (q, rho)."""
    res = run_las(tb, lat.q, lat.q + 1, threads=threads, timings=False, keep_relations=True)
    out = {}
    for line in res["relation_lines"]:
        a, b = map(int, line.split(":", 1)[0].split(","))
        if (a - lat.rho * b) % lat.q == 0:
            out[(a, b)] = line
    return out


def v23(digits: int, nq: int = 5, window=None, threads: int = 2) -> dict:
    tb, P, side = setup(digits)
    f0, f1 = P["f0"], P["f1"]
    sm = Smoother(side)
    lats = qlattices(tb.workdir)
    pick = [lats[k * (len(lats) - 1) // max(nq - 1, 1)] for k in range(nq)]
    window = window if window is not None else WINDOWS.get(digits)
    rows, hits = [], []
    for lat, nrel in pick:
        t0 = time.time()
        w = window or (-lat.I // 2, lat.I // 2, 0, lat.J)
        pts = window_points(lat, *w)
        n0, n1, c0, c1 = classify_points(pts, f0, f1, lat.q, sm)
        s0, s1 = c0["smooth"], c1["smooth"]
        k0, k1 = c0["cado"], c1["cado"]
        npts = len(pts)
        sB = [x and y for x, y in zip(s0, s1)]
        sC = [x and y for x, y in zip(k0, k1)]
        las = las_relations_for_q(tb, lat, threads)
        las_in_w = {}
        for (a, b), line in las.items():
            ij = lat.ij(a, b)
            if ij and w[0] <= ij[0] < w[1] and w[2] <= ij[1] < w[3]:
                las_in_w[(a, b)] = ij
        sC_set = {(pts[k][2], pts[k][3]) for k in range(npts) if sC[k]}
        sB_set = {(pts[k][2], pts[k][3]) for k in range(npts) if sB[k]}
        found = set(las_in_w)
        r = {"q": lat.q, "rho": lat.rho, "I": lat.I, "J": lat.J, "window": list(w), "points": npts,
             "m0": sum(s0), "m1": sum(s1), "joint_SB": len(sB_set), "joint_SCADO": len(sC_set),
             "mc0": sum(k0), "mc1": sum(k1),
             "las_in_window": len(found), "las_total_q": len(las), "las_found_SCADO": len(found & sC_set),
             "las_found_SB": len(found & sB_set), "las_not_in_SB": len(found - sB_set), "secs": None}
        r["indep_ratio_SB"] = (r["joint_SB"] / npts) / ((r["m0"] / npts) * (r["m1"] / npts)) if r["m0"] and r["m1"] else None
        r["indep_ratio_SCADO"] = (r["joint_SCADO"] / npts) / ((r["mc0"] / npts) * (r["mc1"] / npts)) if r["mc0"] and r["mc1"] else None
        r["recall_SCADO"] = r["las_found_SCADO"] / r["joint_SCADO"] if r["joint_SCADO"] else None
        r["secs"] = time.time() - t0
        rows.append(r)
        for k in range(npts):
            if sB[k]:
                i, j, a, b = pts[k]
                hits.append({"q": lat.q, "rho": lat.rho, "i": i, "j": j, "a": a, "b": b, "SB": True,
                             "SCADO": bool(sC[k]), "las": (a, b) in found})
        print(f"[v23 c{digits}] {json.dumps(r)}", file=sys.stderr, flush=True)
    tot = {k: sum(r[k] for r in rows) for k in ("points", "m0", "m1", "mc0", "mc1", "joint_SB", "joint_SCADO",
                                                 "las_in_window", "las_found_SCADO", "las_found_SB", "las_not_in_SB")}
    # pooled independence ratio: sum joint / sum over q of (m0 m1 / points)
    exp_B = sum(r["m0"] * r["m1"] / r["points"] for r in rows)
    exp_C = sum(r["mc0"] * r["mc1"] / r["points"] for r in rows)
    tot["indep_ratio_SB_pooled"] = tot["joint_SB"] / exp_B if exp_B else None
    tot["indep_ratio_SCADO_pooled"] = tot["joint_SCADO"] / exp_C if exp_C else None
    tot["recall_SCADO_pooled"] = tot["las_found_SCADO"] / tot["joint_SCADO"] if tot["joint_SCADO"] else None
    return {"digits": digits, "params": tb.params, "rows": rows, "totals": tot, "hits": hits}


PER_POINT_CATS = ("2 sieving on side 0", "3 sieving (not differentiated)", "4 norms on side 0",
                  "258 sieving on side 1", "260 norms on side 1")
COFACT_CATS = ("5 cofactoring on side 0", "6 cofactoring (not differentiated)", "261 cofactoring on side 1")


TREE_GROUPS = {
    # per (factor-base prime x special-q) and per bucket update: lattice setup + bucket fill
    "bucket_fill": ("fill_in_buckets_toplevel", "fill_in_buckets_lowlevel", "downsort", "allocate_buckets"),
    # per sieve-region point: norm init, small sieve, applying buckets, survivor scan
    "per_point": ("init norms", "small sieve", "prepare small sieve", "apply buckets", "resieve",
                  "search_survivors"),
    # per survivor: trial division (bucket purge, S minus S) and cofactorization (P-1/P+1/ECM chain)
    "trial_division": ("purge buckets", "S minus S (2)"),
    "cofactorization": ("cofactoring_sync", "cofactoring"),
}


def tree_shares(tree: Dict[str, float]) -> Dict[str, float]:
    """Shares of las's -T self-times by our groups. NB: las files the fuzzy 'cofactoring_sync' timer
    under category -1 ('uncategorized'), so category-based cofactorization shares undercount it."""
    agg = {g: 0.0 for g in TREE_GROUPS}
    total = 0.0
    for k, v in tree.items():
        name = k.split(":", 1)[1].split(" ; ")[0].strip()
        total += v
        for g, names in TREE_GROUPS.items():
            if name in names:
                agg[g] += v
    agg["other"] = total - sum(agg.values())
    return {g: v / total for g, v in agg.items()} if total else {}


def v4(digits: int, nq: int = 200, repeats: int = 3, threads: int = 10, dI=(-1, 0, 1),
       keep_relations_at_default: bool = True) -> dict:
    """las CPU / yield vs sieve area on a fixed set of nq special-q starting at qmin."""
    tb, P, side = setup(digits)
    I0 = int(tb.params["I"])
    q0 = tb.qmin
    runs = []
    for d in dI:
        for rep in range(repeats):
            keep = keep_relations_at_default and d == 0 and rep == 0
            out_path = os.path.join(WORKDIR_RAW, f"c{digits}_v4_I{I0}_t{threads}.las.out") if keep else None
            r = run_las(tb, q0, 1 << 40, threads=threads, overrides={"I": str(I0 + d)}, timings=True,
                        out_path=out_path, extra=["-nq", str(nq)])
            cats = r.get("categories", {})
            tot = sum(cats.values()) or None
            r["I"] = I0 + d
            r["area_per_q"] = (1 << (I0 + d)) * r.get("avg_J", 0)
            r["per_point_share"] = sum(cats.get(c, 0.0) for c in PER_POINT_CATS) / tot if tot else None
            r["cofact_share"] = sum(cats.get(c, 0.0) for c in COFACT_CATS) / tot if tot else None
            r["survivors_cofact"] = r["funnel"].get("enter_cofactoring")
            r["tree_shares"] = tree_shares(r.get("tree", {}))
            r["special_q"] = len(r.get("special_q", []))
            r["raw_out"] = out_path
            runs.append(r)
            print(f"[v4 c{digits}] I={I0 + d} rep={rep} cpu={r.get('cpu_s')} rels={r['relations']} "
                  f"pp={r['per_point_share']} cof={r['cofact_share']}", file=sys.stderr, flush=True)
    # fits
    import numpy as np
    byI = {}
    for r in runs:
        byI.setdefault(r["I"], []).append(r)
    pts = []
    for I, rs in sorted(byI.items()):
        cpu = [x["cpu_s"] for x in rs]
        pts.append({"I": I, "area_per_q": rs[0]["area_per_q"], "cpu_mean": float(np.mean(cpu)),
                    "cpu_sd": float(np.std(cpu, ddof=1)) if len(cpu) > 1 else 0.0,
                    "relations": rs[0]["relations"], "unique_ab": rs[0]["unique_ab"],
                    "survivors_cofact": rs[0]["survivors_cofact"],
                    "per_point_share": float(np.mean([x["per_point_share"] for x in rs])),
                    "cofact_share": float(np.mean([x["cofact_share"] for x in rs])),
                    "tree_shares": {g: float(np.mean([x["tree_shares"][g] for x in rs])) for g in rs[0]["tree_shares"]},
                    "n_special_q": rs[0].get("n_special_q")})
    la = np.log([p["area_per_q"] for p in pts])
    exp_cpu = float(np.polyfit(la, np.log([p["cpu_mean"] for p in pts]), 1)[0])
    exp_yield = float(np.polyfit(la, np.log([p["relations"] for p in pts]), 1)[0])
    X = np.array([[p["area_per_q"] * p["n_special_q"], p["survivors_cofact"]] for p in pts], dtype=float)
    y = np.array([p["cpu_mean"] for p in pts])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    # post-hoc (disclosed): CPU = c0 + c1 * area_per_q, i.e. with a per-special-q fixed cost
    A = np.array([p["area_per_q"] for p in pts], dtype=float)
    c1, c0 = np.polyfit(A, y, 1)
    k0 = [p["I"] for p in pts].index(int(tb.params["I"]))
    return {"digits": digits, "params": tb.params, "q0": q0, "nq": nq, "threads": threads, "runs": runs, "by_I": pts,
            "fit_intercept_posthoc": {"c0_fixed_cpu_s": float(c0), "c1_per_area_per_q": float(c1),
                                      "fixed_share_at_default": float(c0 / y[k0]),
                                      "cpu_ms_per_relation_at_default": 1000 * float(y[k0]) / pts[k0]["relations"],
                                      "points_per_relation_at_default": float(A[k0] * nq / pts[k0]["relations"])},
            "exp_cpu_vs_area": exp_cpu, "exp_yield_vs_area": exp_yield,
            "fit_cpu_area_out": {"c_area_per_point": float(coef[0]), "c_out_per_survivor": float(coef[1]),
                                 "area_term_share_at_default": float(coef[0] * X[len(pts) // 2][0] / y[len(pts) // 2])}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("exp", choices=["v1", "v23", "v4"])
    ap.add_argument("--digits", type=int, required=True)
    ap.add_argument("--nq", type=int, default=None)
    ap.add_argument("--per-q", type=int, default=10000)
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.exp == "v4":
        res = v4(a.digits, nq=a.nq or 200, threads=a.threads)
        print(json.dumps({k: res[k] for k in ("by_I", "exp_cpu_vs_area", "exp_yield_vs_area", "fit_cpu_area_out")}, indent=1))
    elif a.exp == "v1":
        res = v1(a.digits, nq=a.nq or 20, per_q=a.per_q)
        print(json.dumps(res["summary"], indent=1))
    else:
        res = v23(a.digits, nq=a.nq or 5)
        hits = res.pop("hits")
        with gzip.open(a.out.replace(".json", "_hits.jsonl.gz"), "wt") as fh:
            for h in hits:
                fh.write(json.dumps(h) + "\n")
        print(json.dumps(res["totals"], indent=1))
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
