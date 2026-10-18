"""HARNESS: pre-registered sign-splitting experiments (notes/signsplit_prereg.md).

The operational searches (hh.py) receive only N, (M, m0, anti, tR, g) and a
seed.  This harness knows the factors and uses them ONLY to (i) classify each
factor discovery as STRUCTURAL / ACCIDENTAL and (ii) verify that the true branch
is covered (coverage check), after the operational run has finished.

Run:
  .venv/bin/python -m src.signsplit.run main [--procs 24]
  .venv/bin/python -m src.signsplit.run summarize
Raw: experiments/raw/signsplit_runs.jsonl ; summary: experiments/analysis/signsplit_summary.json
"""
from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool
from statistics import median
from typing import Dict, List

from src.common.instances import balanced_semiprime
from src.common.metrics import append_jsonl, read_jsonl

from . import accounting as acc
from .hh import (buckets, crt, factor_squarefree, hh_search, lam, prop33, tau,
                 tau0_floor, v2_counts, v2_search, weight_lattice)
from .identity import carry_range, sign_decomposition, sqrt_mod_prime

CONFIG = {
    "bits": [32, 40, 48, 56, 64],
    "seeds": list(range(20)),
    "Ms": [6, 30, 210, 2310],
    "c_ratio": 2,
    "v2_operational_cap": 300_000,
    "base_seed_retries": 3,
}
RAW = "experiments/raw/signsplit_runs.jsonl"
SUMMARY = "experiments/analysis/signsplit_summary.json"


# --------------------------------------------------------------------------
# harness-only checks (use p, q)
# --------------------------------------------------------------------------
def true_branch(N: int, p: int, M: int, m0: int, anti=()) -> Dict:
    q = N // p
    S0 = buckets(N, m0)
    s0 = max(s for s in S0 if s <= p)
    sigma = p % M
    if anti:  # canonical sign on anti primes, as in hh_search
        res, mods = [], []
        for m in factor_squarefree(M):
            r = sigma % m
            if m in anti and r > m // 2:
                r = m - r
            res.append(r)
            mods.append(m)
        sigma = crt(res, mods)
    return {"sigma": sigma, "s0": s0, "q": q}


def coverage_shared(N: int, p: int, M: int, m0: int, anti=()) -> Dict:
    """Recompute the true branch's (a, b, j) and check u = a q + b p is in the
    giant window: u - j in M^2 [0, 2 lambda]."""
    tb = true_branch(N, p, M, m0, anti)
    primes = factor_squarefree(M)
    MA = 1
    for m in anti:
        MA *= m
    I = M * MA * MA
    lmb = lam(N, M, m0, I)
    lat = weight_lattice(N, primes, tb["sigma"], anti)
    from src.common.metrics import Metrics
    a, b = prop33(N, m0, tb["s0"], lat, Metrics())
    M2 = M * M
    j = M2 * (tau0_floor(a, b, N, tb["s0"], M2) - lmb) + tau(a, b, tb["sigma"], N, M)
    u = a * tb["q"] + b * p
    k, r = divmod(u - j, M2)
    return {"covered": r == 0 and 0 <= k <= 2 * lmb, "i_true": k, "sigma": tb["sigma"],
            "s0": tb["s0"], "a": a, "b": b,
            "abs_c": abs(N * a), "abs_d": abs(-N * m0 * a + m0 * tb["s0"] ** 2 * b)}


def classify(res: Dict, cov: Dict) -> str:
    st = res["status"]
    if res.get("factor") is None:
        return "NONE"
    if st in ("collision", "exact_match") and res.get("hit"):
        tgt = res["hit"].get("target")
        if tgt and tgt[0] == cov["sigma"] and tgt[1] == cov["s0"]:
            if res["hit"].get("baby") in (None, cov["i_true"]):
                return "STRUCTURAL"
        return "ACCIDENTAL"
    return "ACCIDENTAL"   # gcd_base, gcd_baby_order, gcd_inverse, gcd_M


def coverage_v2(N: int, p: int, M: int, m0: int, tR: int, g: int) -> Dict:
    q = N // p
    S0 = buckets(N, m0)
    s0 = max(s for s in S0 if s <= p)
    primes = factor_squarefree(M)
    M2 = M * M
    lmb = lam(N, M, m0, M)
    gm = (N * pow(p * p, -1, M)) % M
    rr = []
    for m in primes:
        roots = sqrt_mod_prime((N * pow(gm, -1, m)) % m, m)
        rr.append(min(r for r in roots if r % m) if m > 2 else 1)
    from src.common.metrics import Metrics
    a, b = prop33(N, m0, s0, weight_lattice(N, primes, crt(rr, primes), ()), Metrics())
    dec = sign_decomposition(a, b, rr, primes, N)
    odd = [m for m in primes if m % 2]
    eps = [1 if (p % m) == r else -1 for m, r in zip(primes, rr) if m % 2]
    S = dec["Delta0"] + sum(e * d for e, d in zip(eps, dec["Delta"]))
    clo, chi = carry_range(dec)
    K = 2 * lmb + 1 + chi - clo
    u = a * q + b * p
    k, r = divmod(u - S - M2 * (tau0_floor(a, b, N, s0, M2) - lmb - chi), M2)
    return {"covered": r == 0 and 0 <= k < K, "i_true": k, "gamma": gm, "s0": s0, "K": K}


def classify_v2(res: Dict, cov: Dict, g: int) -> str:
    if res.get("factor") is None:
        return "NONE"
    hit = res.get("hit") or {}
    rt = hit.get("right")
    if res["status"] in ("collision", "exact_match") and rt:
        gm, s0, SR, i = rt
        WR = -(-cov["K"] // g)
        if gm == cov["gamma"] and s0 == cov["s0"] and i == cov["i_true"] % WR:
            # right index i = i_true mod W_R (the left index h absorbs the rest)
            return "STRUCTURAL"
    return "ACCIDENTAL"


# --------------------------------------------------------------------------
# one instance
# --------------------------------------------------------------------------
def run_instance(args) -> List[Dict]:
    bits, seed = args
    inst = balanced_semiprime(bits, seed)
    N, p = inst.N, inst.labels["p"]
    out = []
    for M in CONFIG["Ms"]:
        primes = factor_squarefree(M)
        odd = tuple(m for m in primes if m % 2)
        t = len(odd)
        base = acc.best_shared(N, M, "V0")
        m0 = base["m0"]
        jobs = [("E1", "V0", m0, ()), ("E1", "V1", m0, ()), ("E1", "V3", m0, (3,)),
                ("E1", "V3", m0, odd)]
        for A in ((3,), odd):
            b3 = acc.best_shared(N, M, "V3", A)
            jobs.append(("E2", "V3", b3["m0"], A))
        jobs.append(("E2", "V0", m0, ()))
        for exp, var, mm0, A in jobs:
            if exp == "E2" and var == "V0":
                continue  # identical to E1 V0 (recorded once, used for both)
            for rs in range(CONFIG["base_seed_retries"]):
                t0 = time.time()
                res = hh_search(N, M, mm0, var, anti=A, seed=rs)
                if not str(res["status"]).startswith("fail_"):
                    break
            cov = coverage_shared(N, p, M, mm0, A if var == "V3" else ())
            met = res["met"]
            out.append({"exp": exp, "variant": var, "bits": bits, "seed": seed, "N": N,
                        "M": M, "t": t, "m0": mm0, "anti": list(A), "base_seed": rs,
                        "status": res["status"], "success": res["factor"] in (p, N // p),
                        "class": classify(res, cov), "covered": cov["covered"],
                        "giant": met.get("giant"), "baby": met.get("baby"),
                        "work": met.get("giant") + met.get("baby"),
                        "lambda": met.get("lambda"), "index": met.get("index"),
                        "abs_c_bits": cov["abs_c"].bit_length(), "abs_d_bits": cov["abs_d"].bit_length(),
                        "counts": met["counts"], "wall_s": round(time.time() - t0, 3)})
        # V2: matched (E1) and own optimum (E2)
        for exp, prm in (("E1", acc.v2_matched(N, M, m0)), ("E2", acc.best_v2(N, M))):
            cnt = v2_counts(N, M, prm["m0"], prm["tR"], prm["g"])
            rec = {"exp": exp, "variant": "V2", "bits": bits, "seed": seed, "N": N, "M": M,
                   "t": t, "m0": prm["m0"], "tR": prm["tR"], "g": prm["g"],
                   "work": cnt["work"], "left": cnt["left"], "right": cnt["right"],
                   "formula": prm["formula"], "lambda": cnt["lambda"]}
            cov = coverage_v2(N, p, M, prm["m0"], prm["tR"], prm["g"])
            rec["covered"] = cov["covered"]
            if cnt["work"] <= CONFIG["v2_operational_cap"]:
                t0 = time.time()
                res = v2_search(N, M, prm["m0"], prm["tR"], prm["g"], seed=0)
                rec.update({"operational": True, "status": res["status"],
                            "success": res["factor"] in (p, N // p),
                            "class": classify_v2(res, cov, prm["g"]),
                            "counts": res["met"]["counts"], "wall_s": round(time.time() - t0, 3)})
            else:
                rec.update({"operational": False})
            out.append(rec)
    return out


def main(procs: int = 24):
    jobs = [(b, s) for b in CONFIG["bits"] for s in CONFIG["seeds"]]
    append_jsonl("experiments/configs/signsplit_config.jsonl", {"config": CONFIG, "time": time.time()})
    with Pool(procs) as pool:
        for recs in pool.imap_unordered(run_instance, jobs):
            for r in recs:
                append_jsonl(RAW, r)
            print(recs[0]["bits"], recs[0]["seed"], "done", flush=True)


def summarize():
    rows = list(read_jsonl(RAW))
    out = {}
    for bits in CONFIG["bits"]:
        for M in CONFIG["Ms"]:
            sel = [r for r in rows if r["bits"] == bits and r["M"] == M]
            if not sel:
                continue
            base = {r["seed"]: r["work"] for r in sel if r["variant"] == "V0"}
            ent = {}
            keys = sorted({(r["exp"], r["variant"], tuple(r.get("anti", []))) for r in sel})
            for k in keys:
                rr = [r for r in sel if (r["exp"], r["variant"], tuple(r.get("anti", []))) == k]
                ratios = [r["work"] / base[r["seed"]] for r in rr if r["seed"] in base]
                op = [r for r in rr if r.get("success") is not None]
                ent["%s/%s/%s" % (k[0], k[1], ",".join(map(str, k[2])))] = {
                    "n": len(rr), "median_work": median(r["work"] for r in rr),
                    "median_ratio_vs_V0": median(ratios) if ratios else None,
                    "min_ratio_vs_V0": min(ratios) if ratios else None,
                    "max_ratio_vs_V0": max(ratios) if ratios else None,
                    "covered": sum(bool(r.get("covered")) for r in rr),
                    "operational": len(op),
                    "success": sum(bool(r.get("success")) for r in op),
                    "structural": sum(r.get("class") == "STRUCTURAL" for r in op),
                    "accidental": sum(r.get("class") == "ACCIDENTAL" for r in op),
                    "median_rn_mul": median(r["counts"].get("rn_mul", 0) for r in op) if op else None,
                    "median_gcd": median(r["counts"].get("gcd", 0) for r in op) if op else None,
                    "median_lattice_red": median(r["counts"].get("lattice_red", 0) for r in op) if op else None,
                }
            out[f"{bits}/{M}"] = ent
    # E3: best over the menu (M, anti set, own parameters) per instance
    e3 = {}
    for bits in CONFIG["bits"]:
        per = {}
        for r in rows:
            if r["bits"] != bits:
                continue
            d = per.setdefault(r["seed"], {"V0": None, "V3": None, "V2": None})
            v = r["variant"]
            if v in d and (d[v] is None or r["work"] < d[v]):
                d[v] = r["work"]
        r3 = [d["V3"] / d["V0"] for d in per.values() if d["V0"] and d["V3"]]
        r2 = [d["V2"] / d["V0"] for d in per.values() if d["V0"] and d["V2"]]
        rs = [min(d["V3"], d["V2"]) / d["V0"] for d in per.values() if d["V0"] and d["V3"] and d["V2"]]
        if r3:
            e3[bits] = {"n": len(rs), "median_bestV3_over_bestV0": median(r3),
                        "min_bestV3_over_bestV0": min(r3),
                        "median_bestV2_over_bestV0": median(r2),
                        "median_bestSignSplit_over_bestV0": median(rs),
                        "min_bestSignSplit_over_bestV0": min(rs),
                        "median_bestV0": median(d["V0"] for d in per.values())}
    out["E3_best_over_menu"] = e3
    sz = [e3[b]["median_bestSignSplit_over_bestV0"] for b in (56, 64) if b in e3]
    if len(sz) == 2:
        out["decision"] = ("SUPPORTED" if max(sz) <= 0.5 else
                           "REFUTED" if min(sz) >= 0.9 else "INCONCLUSIVE")
    with open(SUMMARY, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    return out


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "main"
    if cmd == "main":
        procs = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[2] == "--procs" else 24
        main(procs)
    elif cmd == "summarize":
        print(json.dumps(summarize(), indent=1))
