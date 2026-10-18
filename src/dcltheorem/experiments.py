"""Pre-registered experiments for dcltheorem (notes/dcltheorem_prereg.md).

  PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments config
  PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments lemma
  PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments oper
  PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments hh
  PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments exc
  PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments summarize

HARNESS code: the operational algorithm (algorithm.dcl_det, hh_baseline.hh_run)
receives only N; p, q are used afterwards for classification / certificates.
"""
from __future__ import annotations

import json
import os
import sys
import time
from math import gcd, log2
from multiprocessing import Pool
from typing import Dict, List

from src.common.metrics import Metrics, append_jsonl, read_jsonl

from .algorithm import build_targets, choose_m, dcl_det, metric_E
from .harness import classify_hit, lemma_checks, ord2, predict, sample_pair
from .params import Design, lam_full, schedule_design

RAW = "experiments/raw"
CFG = "experiments/configs/dcltheorem_experiments.json"
ANA = "experiments/analysis/dcltheorem_summary.json"

CONFIG = {
    "lemma": {"n": [28, 32, 36, 40], "seeds": list(range(0, 60)), "design": "C1"},
    "oper": {"n": [40, 48, 56, 64, 72, 80], "seeds": list(range(100, 132)),
             "designs": ["C1", "C3"], "C2_n": [40, 48, 56, 64], "C2_seeds": list(range(100, 108))},
    "hh": {"n": [40, 48, 56, 64, 72, 80], "seeds": list(range(100, 116))},
    "exc": {"n": [48, 64, 80, 96], "seeds": list(range(5000, 6000)), "design": "C1"},
    "designs": {
        "C1": "schedule_design(n, 3.0, K=4, mod4=True)",
        "C2": "schedule_design(n, lam_full(n), K=4, mod4=True)",
        "C3": "schedule_design(n, 3.0, K=4, mod4=False)",
    },
    "metric_E": "mulN + mul_small + poly_cost + ceil(log2 n)*(gcdN+invN) + ceil(bigmul_bits/n)",
    "fit_range": [48, 80],
}

_DCACHE: Dict = {}


def get_design(name: str, n: int) -> Design:
    key = (name, n)
    if key not in _DCACHE:
        if name == "C1":
            _DCACHE[key] = schedule_design(n, 3.0, K=4, mod4=True, label="C1")
        elif name == "C2":
            _DCACHE[key] = schedule_design(n, lam_full(n), K=4, mod4=True, label="C2")
        elif name == "C3":
            _DCACHE[key] = schedule_design(n, 3.0, K=4, mod4=False, label="C3")
        else:
            raise ValueError(name)
    return _DCACHE[key]


def _accidental_null(N: int, p: int, q: int, blk: Dict) -> float:
    """Expected number of random giant/baby coincidences mod p or mod q over
    the block: n_giants * m * (1/ord_p(gamma) + 1/ord_q(gamma))."""
    from sympy import n_order
    g = blk["gamma"]
    tot = 0.0
    for r in (p, q):
        gr = g % r
        o = 1 if gr == 1 else int(n_order(gr, r))
        tot += blk["n_giants"] * blk["m"] / o
    return tot


def _strip(blk: Dict) -> Dict:
    return {k: v for k, v in blk.items()
            if k not in ("targets", "giant_meta", "giant_gcds", "gamma", "u")}


def run_one(args) -> Dict:
    mode, dname, n, seed = args
    N, p, q = sample_pair(n, seed)
    D = get_design(dname, n)
    met = Metrics()
    t0 = time.perf_counter()
    r = dcl_det(N, D, met, early_exit=False, record=True)
    wall = time.perf_counter() - t0
    ok = r["factor"] is None or (1 < r["factor"] < N and N % r["factor"] == 0)
    blocks = []
    Es = []
    for b in r["blocks"]:
        E = metric_E(b["met"]["counts"], n)
        Es.append(E)
        rec = _strip(b)
        rec["E"] = E
        if "targets" in b:
            rec["lemma"] = lemma_checks(N, p, q, D, b)
            rec["class"] = classify_hit(N, p, q, b) if b.get("factor") else None
            rec["acc_null"] = _accidental_null(N, p, q, b)
        else:
            rec["class"] = "GCD_" + b.get("status", "") if b.get("factor") else None
        blocks.append(rec)
    ms = [b.get("m", 0) for b in r["blocks"]]
    pr = predict(N, p, q, D, ms)
    first = next((i for i, b in enumerate(r["blocks"]) if b.get("factor")), None)
    out = {"mode": mode, "design": dname, "n": n, "seed": seed, "N": N,
           "design_summary": D.summary(), "success": r["factor"] is not None,
           "correct": ok, "first_block": first,
           "E_first": sum(Es[:first + 1]) if first is not None else None,
           "E_full": sum(Es), "E_block0": Es[0],
           "list_block0": blocks[0].get("m", 0) + blocks[0].get("n_giants", 0),
           "list_full": sum(b.get("m", 0) + b.get("n_giants", 0) for b in blocks),
           "certified": pr["certified"], "cert_blocks": pr["blocks"], "ord_q2": pr["ord_q2"],
           "q_bits": q.bit_length(), "blocks": blocks, "wall_s": round(wall, 4)}
    return out


def run_hh(args) -> Dict:
    from .hh_baseline import hh_E, hh_params, hh_run
    n, seed = args
    N, p, q = sample_pair(n, seed)
    t0 = time.perf_counter()
    prm = hh_params(N)
    r = hh_run(N, prm, seed=seed)
    wall = time.perf_counter() - t0
    c = r["met"]["counts"]
    f = r.get("factor")
    return {"mode": "hh", "n": n, "seed": seed, "N": N, "params": prm,
            "success": bool(f) and 1 < f < N and N % f == 0, "status": r["status"],
            "E": hh_E(c, n), "list": prm["work"], "counts": c, "wall_s": round(wall, 4)}


def run_exc(args) -> Dict:
    """Certificate components only (no stage 3): clean blocks, coverage,
    FULL counts, per block; m from the public choose_m."""
    n, seed = args
    N, p, q = sample_pair(n, seed)
    D = get_design("C1", n)
    ms = []
    for k in range(D.K):
        t = build_targets(N, D, k, Metrics())
        ms.append(choose_m(t, D.mbar[k]))
    pr = predict(N, p, q, D, ms)
    return {"mode": "exc", "n": n, "seed": seed, "certified": pr["certified"],
            "blocks": pr["blocks"], "ord_q2_bits": pr["ord_q2"].bit_length(),
            "q_bits": q.bit_length()}


def _pool_run(fn, jobs, path, procs=24):
    if os.path.exists(path):
        os.remove(path)
    with Pool(procs) as pool:
        for rec in pool.imap_unordered(fn, jobs):
            append_jsonl(path, rec)


def main(mode: str):
    os.makedirs(RAW, exist_ok=True)
    if mode == "config":
        os.makedirs(os.path.dirname(CFG), exist_ok=True)
        cfg = dict(CONFIG)
        cfg["design_summaries"] = {f"{d}@{n}": get_design(d, n).summary()
                                   for d in ("C1", "C3") for n in CONFIG["oper"]["n"]}
        cfg["design_summaries"].update({f"C2@{n}": get_design("C2", n).summary()
                                        for n in CONFIG["oper"]["C2_n"]})
        with open(CFG, "w") as fh:
            json.dump(cfg, fh, indent=1)
        print(json.dumps(cfg["design_summaries"], indent=0)[:3000])
    elif mode == "lemma":
        c = CONFIG["lemma"]
        jobs = [("lemma", c["design"], n, s) for n in c["n"] for s in c["seeds"]]
        _pool_run(run_one, jobs, f"{RAW}/dcltheorem_lemma.jsonl")
    elif mode == "oper":
        c = CONFIG["oper"]
        jobs = [("oper", d, n, s) for d in c["designs"] for n in c["n"] for s in c["seeds"]]
        jobs += [("oper", "C2", n, s) for n in c["C2_n"] for s in c["C2_seeds"]]
        jobs.sort(key=lambda j: -j[2])            # big first
        _pool_run(run_one, jobs, f"{RAW}/dcltheorem_oper.jsonl")
    elif mode == "hh":
        c = CONFIG["hh"]
        jobs = [(n, s) for n in c["n"] for s in c["seeds"]]
        _pool_run(run_hh, jobs, f"{RAW}/dcltheorem_hh.jsonl")
    elif mode == "exc":
        c = CONFIG["exc"]
        jobs = [(n, s) for n in c["n"] for s in c["seeds"]]
        _pool_run(run_exc, jobs, f"{RAW}/dcltheorem_exc.jsonl")
    elif mode == "summarize":
        from .analyze import summarize
        summarize()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main(sys.argv[1])
