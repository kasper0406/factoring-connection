"""Pre-registered experiments for `unbalanced` (notes/unbalanced_prereg.md).

  PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments config
  PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments lemma
  PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments oper1
  PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments oper2
  PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments exc
  PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments summarize

HARNESS code: the operational algorithms (algorithm.dclw, dcltheorem.algorithm.dcl_det)
receive only N; p, q are used afterwards for classification / certificates.
"""
from __future__ import annotations

import json
import os
import sys
import time
from math import exp
from multiprocessing import Pool
from typing import Dict, List

from src.common.metrics import Metrics, append_jsonl
from src.dcltheorem import algorithm as dalg
from src.dcltheorem.harness import classify_hit as d_classify
from src.dcltheorem.params import schedule_design

from .algorithm import build_targets, dclw, metric_E
from .harness import (classify_hit, covering, in_range, lemma_checks, ord2, predict,
                      sample_pair_ratio)
from .params import WideDesign, regime1_Jmax, schedule_regime2, wide_from_C1

RAW = "experiments/raw"
CFG = "experiments/configs/unbalanced_experiments.json"

CONFIG = {
    "lemma": {"n": [32, 40], "seeds": list(range(0, 10)), "js": "0..J1(n)", "design": "W1"},
    "oper1": {"n": [48, 64, 80], "seeds": list(range(100, 124)), "js": [0, 2, 4, "J1"],
              "designs": ["W1", "C1"]},
    "oper2": {"nj": [[64, 7], [64, 9], [80, 7], [80, 9], [80, 11], [96, 7], [96, 9],
                     [96, 11], [96, 13]], "seeds": list(range(200, 216)), "design": "W2"},
    "exc": {"n": [64, 80], "js": [0, 3, "J1"], "seeds": list(range(5000, 5300)), "design": "W1"},
    "designs": {
        "C1": "dcltheorem.params.schedule_design(n, 3.0, K=4) (Theorem D' finite schedule)",
        "W1": "wide_from_C1(n, js=0..J1): H, Y, K=4 from C1; ratio bands j=0..J1, J1 = max j with 2^j <= H^2",
        "W2": "schedule_regime2(n, j, 3.0): weights (1,B), B in [2^j, 2^(j+1)), ratio window (B-Delta, B+Delta], per-band Y_j",
    },
    "metric_E": "dcltheorem.algorithm.metric_E (unchanged)",
}

_DC: Dict = {}


def J1(n: int) -> int:
    return regime1_Jmax(schedule_design(n, 3.0).H) if ("H", n) not in _DC else regime1_Jmax(_DC[("H", n)])


def get_design(name: str, n: int, j: int = 0):
    key = (name, n, j if name == "W2" else 0)
    if key not in _DC:
        if name == "C1":
            _DC[key] = schedule_design(n, 3.0, K=4, mod4=True, label="C1")
        elif name == "W1":
            H = get_design("C1", n).H
            _DC[key] = wide_from_C1(n, list(range(regime1_Jmax(H) + 1)))
        elif name == "W2":
            _DC[key] = schedule_regime2(n, j, 3.0)
        else:
            raise ValueError(name)
    return _DC[key]


def _strip(blk: Dict) -> Dict:
    return {k: v for k, v in blk.items()
            if k not in ("targets", "giant_meta", "giant_gcds", "gamma", "u")}


def _struct_cover(N: int, p: int, q: int, blk: Dict) -> Dict:
    """Design-agnostic structural coverage of a recorded block: targets with
    p == r (mod l) and g = (Bp+Aq-c)/l^2 in [g_lo, g_lo+W)  (p-side), and the
    q-side analogue (q == r, (Bq+Ap-c)/l^2 in window)."""
    cp = cq = 0
    for t in blk.get("targets", []):
        A, B, l, r, c, g_lo, W = t["A"], t["B"], t["l"], t["r"], t["c"], t["g_lo"], t["W"]
        l2 = l * l
        if p % l == r:
            x = B * p + A * q - c
            if x % l2 == 0 and g_lo <= x // l2 < g_lo + W:
                cp += 1
        if q % l == r:
            x = B * q + A * p - c
            if x % l2 == 0 and g_lo <= x // l2 < g_lo + W:
                cq += 1
    return {"cover_p": cp, "cover_q": cq}


def run_wide(args) -> Dict:
    mode, dname, n, j, seed = args
    N, p, q = sample_pair_ratio(n, j, seed)
    D = get_design(dname, n, j)
    d = ord2(q)
    met = Metrics()
    t0 = time.perf_counter()
    r = dclw(N, D, met, early_exit=False, record=True)
    wall = time.perf_counter() - t0
    ok = r["factor"] is None or (1 < r["factor"] < N and N % r["factor"] == 0)
    blocks, Es = [], []
    for b in r["blocks"]:
        E = metric_E(b["met"]["counts"], n)
        Es.append(E)
        rec = _strip(b)
        rec["E"] = E
        if "targets" in b:
            rec["lemma"] = lemma_checks(N, p, q, D, b, d)
            rec["class"] = classify_hit(N, p, q, b) if b.get("factor") else None
            rec.update(_struct_cover(N, p, q, b))
        else:
            rec["class"] = "GCD_" + b.get("status", "") if b.get("factor") else None
        blocks.append(rec)
    ms = [b.get("m", 0) for b in r["blocks"]]
    pr = predict(N, p, q, D, ms, d)
    first = next((i for i, b in enumerate(r["blocks"]) if b.get("factor")), None)
    return {"mode": mode, "design": dname, "n": n, "j": j, "seed": seed, "N": N,
            "ratio_bits": round((q / p).bit_length() - 1 + 0.0, 1),
            "design_summary": D.summary(), "success": r["factor"] is not None, "correct": ok,
            "first_block": first, "E_first": sum(Es[:first + 1]) if first is not None else None,
            "E_full": sum(Es), "E_block0": Es[0] if Es else None,
            "list_block0": blocks[0].get("m", 0) + blocks[0].get("n_giants", 0) if blocks else None,
            "certified": pr["certified"], "cert_blocks": pr["blocks"], "ord_q2_bits": d.bit_length(),
            "q_bits": q.bit_length(), "p_bits": p.bit_length(), "blocks": blocks,
            "wall_s": round(wall, 3)}


def run_c1(args) -> Dict:
    """Theorem D' (C1) on a ratio-band instance: same N as run_wide with the same (n, j, seed)."""
    mode, _, n, j, seed = args
    N, p, q = sample_pair_ratio(n, j, seed)
    D = get_design("C1", n)
    met = Metrics()
    t0 = time.perf_counter()
    r = dalg.dcl_det(N, D, met, early_exit=False, record=True)
    wall = time.perf_counter() - t0
    ok = r["factor"] is None or (1 < r["factor"] < N and N % r["factor"] == 0)
    blocks, Es = [], []
    for b in r["blocks"]:
        E = metric_E(b["met"]["counts"], n)
        Es.append(E)
        rec = _strip(b)
        rec["E"] = E
        if "targets" in b:
            rec["class"] = d_classify(N, p, q, b) if b.get("factor") else None
            rec.update(_struct_cover(N, p, q, b))
        else:
            rec["class"] = "GCD_" + b.get("status", "") if b.get("factor") else None
        blocks.append(rec)
    first = next((i for i, b in enumerate(r["blocks"]) if b.get("factor")), None)
    return {"mode": mode, "design": "C1", "n": n, "j": j, "seed": seed, "N": N,
            "success": r["factor"] is not None, "correct": ok, "first_block": first,
            "E_first": sum(Es[:first + 1]) if first is not None else None,
            "E_full": sum(Es), "E_block0": Es[0] if Es else None, "blocks": blocks,
            "wall_s": round(wall, 3)}


def run_exc(args) -> Dict:
    """Certificate components only (no stage 3), W1."""
    n, j, seed = args
    N, p, q = sample_pair_ratio(n, j, seed)
    D = get_design("W1", n)
    ms = []
    for k in range(D.K):
        t = build_targets(N, D, k, Metrics())
        from src.dcltheorem.algorithm import choose_m
        ms.append(choose_m(t, 0))
    pr = predict(N, p, q, D, ms)
    lam = []
    for k in range(D.K):
        s = sum(1.0 / l for l in D.bands[k])
        lam.append(sum(1 for t in D.insts if in_range(N, p, t)) * s)
    return {"mode": "exc", "n": n, "j": j, "seed": seed, "certified": pr["certified"],
            "blocks": pr["blocks"], "lam": lam, "ord_q2_bits": pr["ord_q2"].bit_length(),
            "q_bits": q.bit_length()}


def _pool(fn, jobs, path, procs=30):
    if os.path.exists(path):
        os.remove(path)
    with Pool(procs, maxtasksperchild=4) as pool:
        for rec in pool.imap_unordered(fn, jobs):
            append_jsonl(path, rec)


def predictions() -> Dict:
    """n-only predictions frozen in the config (computed before any run)."""
    out = {}
    for n in [48, 64, 80]:
        D = get_design("C1", n)
        W = get_design("W1", n)
        s0 = sum(1.0 / l for l in W.bands[0])
        J = regime1_Jmax(D.H)
        lamW = {}
        for j in range(J + 1):
            nv = sum(1 for t in W.insts if t.band == j)
            lamW[j] = round(nv * s0, 3)
        # C1: weight (A,B) has window [Lo, max(F(a0), F(b0))]; for B >= sqrt2 A the
        # window contains F(p) iff N/p^2 <= (B/A)^2 (and N/p^2 >= 1): count weights
        # covering the band-centre ratio 1.5 * 2^j.
        lamC = {}
        for j in range(J + 1):
            Rc = 1.5 * (1 << j)
            cnt = 0
            for (A, B) in D.V:
                # exact window test at a reference N = 3*2^(n-2), p = sqrt(N/Rc)
                Nr = 3 << (n - 2)
                pr = (Nr / Rc) ** 0.5
                a0 = (Nr / 2) ** 0.5
                b0 = Nr ** 0.5
                F = lambda x: B * x + A * Nr / x
                if F(pr) <= max(F(a0), F(b0)):
                    cnt += 1
            lamC[j] = round(cnt * sum(1.0 / l for l in D.bands[0]), 3)
        out[n] = {"H": D.H, "Y": D.Y, "J1": J, "lam_W1_block0": lamW,
                  "lam_C1_block0_at_band_centre": lamC,
                  "uncovered_W1_block0_poisson": {j: round(exp(-v), 4) for j, v in lamW.items()}}
    r2 = {}
    from .params import _band_primes, _expected_cost_r2
    for n, j in CONFIG["oper2"]["nj"]:
        W = get_design("W2", n, j)
        r2[f"{n},{j}"] = {"Y": W.Y, "Delta": W.Delta, "n_insts": len(W.insts),
                          "band_sizes": W.summary()["band_sizes"],
                          "expected_block0_list": _expected_cost_r2(n, j, W.Delta, W.bands[0])}
    out["regime2"] = r2
    return out


def main(mode: str):
    os.makedirs(RAW, exist_ok=True)
    if mode == "config":
        os.makedirs(os.path.dirname(CFG), exist_ok=True)
        cfg = dict(CONFIG)
        cfg["predictions"] = predictions()
        cfg["design_summaries"] = {}
        for n in [32, 40, 48, 64, 80]:
            cfg["design_summaries"][f"W1@{n}"] = get_design("W1", n).summary()
            cfg["design_summaries"][f"C1@{n}"] = get_design("C1", n).summary()
        with open(CFG, "w") as fh:
            json.dump(cfg, fh, indent=1, default=str)
        print(json.dumps(cfg["predictions"], indent=1, default=str))
    elif mode == "lemma":
        c = CONFIG["lemma"]
        jobs = [("lemma", "W1", n, j, s) for n in c["n"] for j in range(J1(n) + 1) for s in c["seeds"]]
        _pool(run_wide, jobs, f"{RAW}/unbalanced_lemma.jsonl")
    elif mode == "oper1":
        c = CONFIG["oper1"]
        jobs = []
        for n in c["n"]:
            js = sorted({J1(n) if j == "J1" else j for j in c["js"]})
            for j in js:
                for s in c["seeds"]:
                    jobs.append(("oper1", "W1", n, j, s))
        jobs.sort(key=lambda t: -t[2])
        _pool(run_wide, jobs, f"{RAW}/unbalanced_oper1_W1.jsonl")
        _pool(run_c1, jobs, f"{RAW}/unbalanced_oper1_C1.jsonl")
    elif mode == "oper2":
        c = CONFIG["oper2"]
        jobs = [("oper2", "W2", n, j, s) for (n, j) in c["nj"] for s in c["seeds"]]
        jobs.sort(key=lambda t: (-t[2], t[3]))
        _pool(run_wide, jobs, f"{RAW}/unbalanced_oper2.jsonl", procs=24)
    elif mode == "exc":
        c = CONFIG["exc"]
        jobs = []
        for n in c["n"]:
            js = sorted({J1(n) if j == "J1" else j for j in c["js"]})
            jobs += [(n, j, s) for j in js for s in c["seeds"]]
        _pool(run_exc, jobs, f"{RAW}/unbalanced_exc.jsonl")
    elif mode == "summarize":
        from .analyze import summarize
        summarize()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main(sys.argv[1])
