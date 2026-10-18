"""Pre-registered experiments for multimod (notes/multimod_prereg.md).

  .venv/bin/python -m src.multimod.experiments model      # E-M model counts
  .venv/bin/python -m src.multimod.experiments coverage   # E-C exhaustive coverage
  .venv/bin/python -m src.multimod.experiments mc         # E-S sampled success
  .venv/bin/python -m src.multimod.experiments oper       # E-O operational runs
  .venv/bin/python -m src.multimod.experiments crt        # CRT index check (Thm C)
  .venv/bin/python -m src.multimod.experiments summarize

HARNESS: instance labels (p) are used only for classification / coverage
checks, never passed to the operational search.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from math import gcd

import numpy as np

from src.common.instances import balanced_semiprime
from src.common.metrics import append_jsonl, read_jsonl
from src.cycloharvey.calibration import best_harvey, best_hh
from src.cycloharvey.operational import harvey_search
from .dcl import best_design, candidate_range, uncovered_candidates, verify_windows
from .operational import dcl_search

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
ANA = os.path.join(ROOT, "experiments", "analysis")
CFG = os.path.join(ROOT, "experiments", "configs", "multimod_experiments.json")

CONFIG = {
    "Hmax": 12,
    "regimes": {"R1": "lam = 3.0", "R2": "lam = ln(sqrt N) + 2"},
    "EM": {"bits": [40, 48, 56, 64, 72, 80, 88, 96], "fit_bits": [56, 64, 72, 80, 88, 96],
           "N": "2^n + 2^(n/2) + 1"},
    "EC": {"bits": [32, 36, 40, 44, 48], "seeds": [0, 1, 2, 3], "special": True,
           "lam_star_grid": [x / 2 for x in range(2, 61)]},
    "ES": {"bits": {"56": 2000, "64": 2000, "80": 1000, "96": 300},
           "m": "baby length m*(n, regime) fixed from the E-M special N of the same size"},
    "EO": {"bits": [40, 48, 56, 64], "seeds": list(range(16))},
}


def lam_of(regime: str, N: int) -> float:
    if regime == "R1":
        return 3.0
    if regime == "R2":
        return math.log(math.isqrt(N)) + 2.0
    raise ValueError(regime)


def _write_cfg():
    os.makedirs(os.path.dirname(CFG), exist_ok=True)
    with open(CFG, "w") as f:
        json.dump(CONFIG, f, indent=1)


def _fresh(name):
    out = os.path.join(RAW, name)
    if os.path.exists(out):
        os.remove(out)
    return out


# ---------------------------------------------------------------------------
def em():
    out = _fresh("multimod_model.jsonl")
    H = CONFIG["Hmax"]
    for n in CONFIG["EM"]["bits"]:
        N = 2 ** n + 2 ** (n // 2) + 1
        rec = {"bits": n}
        for reg in ("R1", "R2"):
            t = time.time()
            d = best_design(N, lam_of(reg, N), H)
            rec[reg] = {k: d[k] for k in ("lam", "m", "giants", "cost", "n_targets",
                                          "n_moduli", "Mmin", "Mmax", "Hmax_used")}
            rec[reg]["design_s"] = round(time.time() - t, 2)
        rec["HH"] = best_hh(N)["cost"]
        rec["Harvey"] = best_harvey(N)["cost"]
        append_jsonl(out, rec)
        print(n, rec["R1"]["cost"], rec["R2"]["cost"], rec["HH"], rec["Harvey"], flush=True)


def _instances(bits, seeds, special):
    out = [("special", 2 ** bits + 2 ** (bits // 2) + 1, None)] if special else []
    for s in seeds:
        inst = balanced_semiprime(bits, s)
        out.append((s, inst.N, inst.labels["p"]))
    return out


def _ec_one(args):
    n, seed, N, p = args
    H = CONFIG["Hmax"]
    rec = {"bits": n, "seed": seed}
    for reg in ("R1", "R2"):
        lam = lam_of(reg, N)
        d = best_design(N, lam, H)
        lo, hi, unc = uncovered_candidates(N, d["targets"])
        rec[reg] = {"lam": d["lam"], "cost": d["cost"], "candidates": int(len(unc)),
                    "uncovered": int(unc.sum()),
                    "poisson_pred": float(len(unc) * math.exp(-d["lam"])),
                    "window_violations": verify_windows(N, d["targets"]),
                    "true_p_covered": (None if p is None else bool(not unc[p - lo]))}
    # diagnostic: smallest lam on the grid with full coverage
    lam_star = None
    for lam in CONFIG["EC"]["lam_star_grid"]:
        d = best_design(N, lam, H)
        _, _, unc = uncovered_candidates(N, d["targets"])
        if not unc.any():
            lam_star = {"lam_target": lam, "lam": d["lam"], "cost": d["cost"]}
            break
    rec["lam_star"] = lam_star
    return rec


def ec(procs: int = 25):
    from multiprocessing import Pool
    out = _fresh("multimod_coverage.jsonl")
    jobs = []
    for n in CONFIG["EC"]["bits"]:
        for seed, N, p in _instances(n, CONFIG["EC"]["seeds"], CONFIG["EC"]["special"]):
            jobs.append((n, seed, N, p))
    with Pool(procs) as pool:
        recs = pool.map(_ec_one, jobs, chunksize=1)
    for rec in recs:
        append_jsonl(out, rec)
        print(rec["bits"], rec["seed"], rec["R1"]["uncovered"], rec["R2"]["uncovered"],
              rec["lam_star"], flush=True)


def _es_one(args):
    n, s, ms = args
    from .dcl import design
    inst = balanced_semiprime(n, s)
    N, p = inst.N, inst.labels["p"]
    rec = {"bits": n, "seed": s}
    for reg in ("R1", "R2"):
        d = design(N, lam_of(reg, N), ms[reg], CONFIG["Hmax"])
        cov = any(p % t["M"] == t["r"] for t in d["targets"])
        rec[reg] = {"covered": cov, "cost": d["cost"], "lam": d["lam"], "m": d["m"]}
    return rec


def es(procs: int = 24):
    from multiprocessing import Pool
    out = _fresh("multimod_mc.jsonl")
    em_rows = {r["bits"]: r for r in read_jsonl(os.path.join(RAW, "multimod_model.jsonl"))}
    for ns, cnt in CONFIG["ES"]["bits"].items():
        n = int(ns)
        ms = {reg: em_rows[n][reg]["m"] for reg in ("R1", "R2")}
        with Pool(procs) as pool:
            recs = pool.map(_es_one, [(n, s, ms) for s in range(cnt)], chunksize=4)
        for rec in recs:
            append_jsonl(out, rec)
        print(n, len(recs), {reg: sum(r[reg]["covered"] for r in recs) for reg in ("R1", "R2")},
              flush=True)


def _classify(N, res, p, d_targets, m):
    """STRUCTURAL iff the colliding target's equation holds for the factor x
    actually found (x = p or q; the q-side is the same mechanism with the roles
    of the factors swapped): x == r (mod M) and
    g_lo + i + k m == (B x + A (N/x) - c)/M^2.  PM1: factor from a gcd on
    gamma's order (baby_order / gamma_order).  Else ACCIDENTAL.
    (Corrected after the first E-O run, which checked only x = p; see notes.)"""
    if res["factor"] not in (p, N // p):
        return None
    st = res["status"]
    if st in ("gcd_baby_order", "gcd_gamma_order"):
        return "PM1"
    if st.startswith("gcd_"):
        return "ACCIDENTAL"
    hit = res.get("hit")
    if hit is None:
        return "STRUCTURAL" if st == "exact_match" else "ACCIDENTAL"
    t = d_targets[hit["target"]]
    x = res["factor"]
    y = N // x
    A, B, M, r, c, g_lo = t["A"], t["B"], t["M"], t["r"], t["c"], t["g_lo"]
    if x % M != r or (B * x + A * y - c) % (M * M):
        return "ACCIDENTAL"
    g = (B * x + A * y - c) // (M * M)
    return "STRUCTURAL" if g == g_lo + hit["baby"] + hit["k"] * m else "ACCIDENTAL"


def _eo_one(args):
    n, s = args
    H = CONFIG["Hmax"]
    inst = balanced_semiprime(n, s)
    N, p = inst.N, inst.labels["p"]
    rec = {"bits": n, "seed": s}
    for reg in ("R1", "R2"):
        lam = lam_of(reg, N)
        d = best_design(N, lam, H)
        cov = any(p % t["M"] == t["r"] for t in d["targets"])
        t0 = time.time()
        res = dcl_search(N, lam, s, H, m=d["m"], early_exit=False)
        wall = time.time() - t0
        rec[reg] = {"success": res["factor"] in (p, N // p), "status": res["status"],
                    "class": _classify(N, res, p, d["targets"], d["m"]),
                    "true_p_covered": cov, "model_cost": d["cost"],
                    "model_giant": d["giants"], "impl_giant": res["met"].get("giant"),
                    "baby": d["m"], "counts": res["met"]["counts"], "wall_s": round(wall, 3)}
    b = best_harvey(N)
    t0 = time.time()
    hv = harvey_search(N, b["r"], b["m"], s, early_exit=False)
    rec["Harvey"] = {"success": hv["factor"] in (p, N // p), "status": hv["status"],
                     "model_cost": b["cost"], "impl_giant": hv["met"].get("giant"),
                     "baby": b["m"], "counts": hv["met"]["counts"],
                     "wall_s": round(time.time() - t0, 3)}
    return rec


def eo(procs: int = 16):
    from multiprocessing import Pool
    out = _fresh("multimod_oper.jsonl")
    for n in CONFIG["EO"]["bits"]:
        with Pool(procs) as pool:
            recs = pool.map(_eo_one, [(n, s) for s in CONFIG["EO"]["seeds"]], chunksize=1)
        for rec in recs:
            append_jsonl(out, rec)
            print(n, rec["seed"], rec["R1"]["success"], rec["R1"]["class"], rec["R1"]["model_cost"],
                  rec["R2"]["success"], rec["R2"]["class"], rec["R2"]["model_cost"],
                  rec["Harvey"]["success"], rec["Harvey"]["model_cost"], flush=True)


def crt():
    """Theorem C check: [L_{r,M1 M2}:Lambda] = [L_{r1,M1}:Lambda][L_{r2,M2}:Lambda]
    computed directly by HNF (cycloharvey.lattice), small cases."""
    from src.cycloharvey.lattice import class_lattice_index
    out = _fresh("multimod_crt.jsonl")
    N = 1000003 * 1000033
    for I, emax in (((-1, 1), 3), ((-2, -1, 1), 4)):
        for M1, M2 in ((5, 7), (3, 7), (5, 11)):
            if len(I) == 3 and M1 * M2 > 35:
                continue
            M = M1 * M2
            for r in [x for x in range(1, M) if gcd(x, M) == 1][:4]:
                a = class_lattice_index(I, N, r, M, emax=emax)["index"]
                b1 = class_lattice_index(I, N, r % M1, M1, emax=emax)["index"]
                b2 = class_lattice_index(I, N, r % M2, M2, emax=emax)["index"]
                rec = {"I": list(I), "M1": M1, "M2": M2, "r": r, "index": a,
                       "product": b1 * b2, "equal": a == b1 * b2}
                append_jsonl(out, rec)
                print(rec, flush=True)


def _fit(xs, ys):
    xs, ys = np.array(xs, float), np.log2(np.array(ys, float))
    return float(np.polyfit(xs, ys, 1)[0])


def summarize():
    res = {}
    f = os.path.join(RAW, "multimod_model.jsonl")
    if os.path.exists(f):
        rows = read_jsonl(f)
        fb = CONFIG["EM"]["fit_bits"]
        sel = [r for r in rows if r["bits"] in fb]
        res["EM"] = {k: _fit([r["bits"] for r in sel],
                             [(r[k]["cost"] if isinstance(r[k], dict) else r[k]) for r in sel])
                     for k in ("R1", "R2", "HH", "Harvey")}
        res["EM"]["rows"] = [{"bits": r["bits"], "R1": r["R1"]["cost"], "R2": r["R2"]["cost"],
                              "HH": r["HH"], "Harvey": r["Harvey"]} for r in rows]
    f = os.path.join(RAW, "multimod_coverage.jsonl")
    if os.path.exists(f):
        rows = read_jsonl(f)
        res["EC"] = [{"bits": r["bits"], "seed": r["seed"],
                      "R1_unc": r["R1"]["uncovered"], "R1_pred": round(r["R1"]["poisson_pred"], 1),
                      "R2_unc": r["R2"]["uncovered"], "cand": r["R1"]["candidates"],
                      "winviol": r["R1"]["window_violations"] + r["R2"]["window_violations"],
                      "lam_star": r["lam_star"]} for r in rows]
    f = os.path.join(RAW, "multimod_mc.jsonl")
    if os.path.exists(f):
        rows = read_jsonl(f)
        out = {}
        for n in sorted({r["bits"] for r in rows}):
            rr = [r for r in rows if r["bits"] == n]
            out[n] = {reg: {"success": sum(r[reg]["covered"] for r in rr) / len(rr),
                            "count": len(rr),
                            "median_cost": float(np.median([r[reg]["cost"] for r in rr]))}
                      for reg in ("R1", "R2")}
        res["ES"] = out
    f = os.path.join(RAW, "multimod_oper.jsonl")
    if os.path.exists(f):
        rows = read_jsonl(f)
        out = {}
        for n in sorted({r["bits"] for r in rows}):
            rr = [r for r in rows if r["bits"] == n]
            d = {}
            for reg in ("R1", "R2"):
                cls = {}
                for r in rr:
                    cls[str(r[reg]["class"])] = cls.get(str(r[reg]["class"]), 0) + 1
                d[reg] = {"success": sum(r[reg]["success"] for r in rr), "n": len(rr),
                          "classes": cls,
                          "covered": sum(r[reg]["true_p_covered"] for r in rr),
                          "list_eq_model": sum((r[reg]["impl_giant"] == r[reg]["model_giant"])
                                               for r in rr if r[reg]["impl_giant"] is not None),
                          "list_runs": sum(r[reg]["impl_giant"] is not None for r in rr),
                          "mean_model_cost": float(np.mean([r[reg]["model_cost"] for r in rr])),
                          "mean_mults": float(np.mean([r[reg]["counts"].get("rn_mul", 0)
                                                       + r[reg]["counts"].get("rn_mul_pow", 0)
                                                       for r in rr])),
                          "mean_wall": float(np.mean([r[reg]["wall_s"] for r in rr]))}
            d["Harvey"] = {"success": sum(r["Harvey"]["success"] for r in rr),
                           "mean_model_cost": float(np.mean([r["Harvey"]["model_cost"] for r in rr])),
                           "mean_mults": float(np.mean([r["Harvey"]["counts"].get("rn_mul", 0)
                                                        + r["Harvey"]["counts"].get("rn_mul_pow", 0)
                                                        for r in rr])),
                           "mean_wall": float(np.mean([r["Harvey"]["wall_s"] for r in rr]))}
            out[n] = d
        res["EO"] = out
    os.makedirs(ANA, exist_ok=True)
    with open(os.path.join(ANA, "multimod_summary.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    _write_cfg()
    {"model": em, "coverage": ec, "mc": es, "oper": eo, "crt": crt,
     "summarize": summarize}.get(sys.argv[1], lambda: None)()


# ---------------------------------------------------------------------------
# POST-HOC diagnostic (added after the E-O run; not pre-registered).  HARNESS:
# uses p, q and sympy's n_order on the PRIME p (never on N).
def diag():
    import random as _r
    from sympy import n_order
    out = _fresh("multimod_oper_diag.jsonl")
    rows = read_jsonl(os.path.join(RAW, "multimod_oper.jsonl"))
    H = CONFIG["Hmax"]
    for rec in rows:
        n, s = rec["bits"], rec["seed"]
        inst = balanced_semiprime(n, s)
        N, p = inst.N, inst.labels["p"]
        q = N // p
        o = {"bits": n, "seed": s}
        for reg in ("R1", "R2"):
            d = best_design(N, lam_of(reg, N), H)
            L = 1
            for M in sorted({t["M"] for t in d["targets"]}):
                L *= M * M
            beta = _r.Random(f"multimod-dcl|{N}|{s}").randrange(2, N - 1)
            op = n_order(pow(beta, L, p), p) if pow(beta, L, p) != 1 else 1
            oq = n_order(pow(beta, L, q), q) if pow(beta, L, q) != 1 else 1
            pairs = d["m"] * d["giants"]
            o[reg] = {"class": rec[reg]["class"], "covered": rec[reg]["true_p_covered"],
                      "ord_p_gamma": op, "ord_q_gamma": oq,
                      "log2_gcd_L_pm1": round(math.log2(gcd(L, p - 1)), 2),
                      "null_accidental": pairs / op + pairs / oq,
                      "m": d["m"]}
        append_jsonl(out, o)
        print(n, s, {k: (o[k]["class"], o[k]["ord_p_gamma"], round(o[k]["null_accidental"], 3))
                     for k in ("R1", "R2")}, flush=True)


if __name__ == "__main__" and sys.argv[1] == "diag":
    diag()
