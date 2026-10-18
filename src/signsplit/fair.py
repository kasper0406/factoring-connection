"""POST-HOC measurements requested by the orchestrator after the pre-registered run
(not part of the pre-registered decision).

(1) Fair baseline: V0 and V1 with naive square-and-multiply giants vs fixed-base
    windowed exponentiation (V0) / global Gray-code step tables (V1), at the E1
    parameters, on the same 100 instances.  Stage-3 counters (product tree) are
    reported alongside.
(2) 2-adic pairing V4 (count-only): M8 = 8 * M_odd.  The 4 roots {+-r, +-r + 4} of
    one 2-adic square class share the HH lattice (r^2 == 1 mod 8).  If a is odd
    the offset shift is +-16 E_2 (verify_two_adic), so one giant plus baby offsets
    {0, 16 E_2} covers the pair r, r+4; if a is even and the shift is 0 the pair
    collapses for free; otherwise 2 giants.  Baby list doubled (+1 for wrap).

Run: .venv/bin/python -m src.signsplit.fair [--procs 24]
Raw: experiments/raw/signsplit_fair.jsonl ; summary printed and written to
experiments/analysis/signsplit_fair.json
"""
from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from statistics import median

from src.common.instances import balanced_semiprime
from src.common.metrics import Metrics, append_jsonl

from . import accounting as acc
from .hh import (buckets, crt, factor_squarefree, hh_search, lam, prop33, square_classes,
                 tau, weight_lattice)
from .identity import sqrt_mod_prime

BITS = [32, 40, 48, 56, 64]
SEEDS = list(range(20))
MS = [6, 30, 210, 2310]
RAW = "experiments/raw/signsplit_fair.jsonl"
OUT = "experiments/analysis/signsplit_fair.json"


def factor_any(M):
    out, x, d = [], M, 2
    while d * d <= x:
        while x % d == 0:
            out.append(d)
            x //= d
        d += 1
    if x > 1:
        out.append(x)
    return out


def v4_counts(N: int, M8: int, m0: int) -> dict:
    """Exact V4 list counts (needs one Prop. 3.3 reduction per (gamma, sigma0))."""
    met = Metrics()
    Modd = M8 // 8
    odd = factor_squarefree(Modd) if Modd > 1 else []
    t = len(odd)
    S0 = buckets(N, m0)
    giants = 0
    stats = {"odd_a": 0, "even_shift0": 0, "even_shift32": 0}
    # square classes mod M8: gamma = N / s^2
    cls = {}
    for s in range(1, M8, 2):
        if all(s % m for m in odd):
            cls.setdefault((N * pow(s * s, -1, M8)) % M8, s)
    for gm, s in cls.items():
        lat = weight_lattice_pp(N, M8, s)
        for s0 in S0:
            a, b = prop33(N, m0, s0, lat, met)
            if a % 2:
                giants += 2 * (1 << t)
                stats["odd_a"] += 1
            else:
                d = (tau(a, b, (s + 4) % 8, N, 8) - tau(a, b, s % 8, N, 8)) % 64
                if d == 0:
                    giants += 2 * (1 << t)
                    stats["even_shift0"] += 1
                else:
                    giants += 4 * (1 << t)
                    stats["even_shift32"] += 1
    lmb = lam(N, M8, m0, M8)
    baby = 2 * (2 * lmb + 2)
    return {"giant": giants, "baby": baby, "work": giants + baby, **stats}


def weight_lattice_pp(N: int, M: int, s: int):
    """HH lock lattice b == (N/s^2) a (mod M) for any (not nec. squarefree) M."""
    return 1, (N * pow(s * s, -1, M)) % M, M


def v0_counts_pp(N: int, M: int, m0: int) -> dict:
    phiM = sum(1 for s in range(1, M) if __import__("math").gcd(s, M) == 1)
    lmb = lam(N, M, m0, M)
    g = phiM * len(buckets(N, m0))
    return {"giant": g, "baby": 2 * lmb + 1, "work": g + 2 * lmb + 1}


def best_over_m0(fn, N, M):
    best = None
    for m0 in acc.m0_grid(N, M):
        c = fn(N, M, m0)
        if best is None or c["work"] < best[1]["work"]:
            best = (m0, c)
    m00 = best[0]
    for m0 in range(max(1, m00 - 3), m00 + 4):
        c = fn(N, M, m0)
        if c["work"] < best[1]["work"]:
            best = (m0, c)
    return {"m0": best[0], **best[1]}


def run_instance(args):
    bits, seed = args
    inst = balanced_semiprime(bits, seed)
    N, p = inst.N, inst.labels["p"]
    out = []
    for M in MS:
        m0 = acc.best_shared(N, M, "V0")["m0"]
        for var in ("V0", "V1"):
            for expo in ("naive", "fixed"):
                r = hh_search(N, M, m0, var, seed=0, expo=expo)
                c = r["met"]["counts"]
                out.append({"kind": "fair", "bits": bits, "seed": seed, "M": M, "m0": m0,
                            "variant": var, "expo": expo, "status": r["status"],
                            "success": r["factor"] in (p, N // p),
                            "rn_mul": c.get("rn_mul", 0), "rn_pow": c.get("rn_pow", 0),
                            "lattice_red": c.get("lattice_red", 0),
                            "poly_mul": c.get("poly_mul", 0),
                            "poly_mul_dimsum": c.get("poly_mul_dimsum", 0),
                            "giant": r["met"].get("giant"), "baby": r["met"].get("baby")})
    # 2-adic pairing, count-only, own optimum vs HH (M = 2 M_odd and M = 8 M_odd)
    for Modd in (1, 3, 15, 105, 1155):
        if bits <= 40 and Modd == 1155:
            continue
        v0_2 = best_over_m0(v0_counts_pp, N, 2 * Modd)
        v0_8 = best_over_m0(v0_counts_pp, N, 8 * Modd)
        # V4: search m0 near V0(8 M_odd)'s optimum (lattice reductions are costly)
        best = None
        for f in (0.5, 0.71, 1.0, 1.41, 2.0):
            mm = max(1, round(v0_8["m0"] * f))
            c = v4_counts(N, 8 * Modd, mm)
            if best is None or c["work"] < best["work"]:
                best = {"m0": mm, **c}
        out.append({"kind": "two_adic", "bits": bits, "seed": seed, "M_odd": Modd,
                    "V0_2": v0_2, "V0_8": v0_8, "V4_8": best})
    return out


def main(procs=24):
    jobs = [(b, s) for b in BITS for s in SEEDS]
    with Pool(procs) as pool:
        for recs in pool.imap_unordered(run_instance, jobs):
            for r in recs:
                append_jsonl(RAW, r)
    summarize()


def summarize():
    rows = [json.loads(l) for l in open(RAW)]
    out = {"fair": {}, "two_adic": {}}
    for b in BITS:
        for M in MS:
            sel = [r for r in rows if r["kind"] == "fair" and r["bits"] == b and r["M"] == M]
            if not sel:
                continue
            d = {}
            for s in SEEDS:
                x = {(r["variant"], r["expo"]): r for r in sel if r["seed"] == s}
                if len(x) < 4:
                    continue
                d.setdefault("V1fixed_over_V0fixed", []).append(x[("V1", "fixed")]["rn_mul"] / x[("V0", "fixed")]["rn_mul"])
                d.setdefault("V1fixed_over_V0naive", []).append(x[("V1", "fixed")]["rn_mul"] / x[("V0", "naive")]["rn_mul"])
                d.setdefault("V0fixed_over_V0naive", []).append(x[("V0", "fixed")]["rn_mul"] / x[("V0", "naive")]["rn_mul"])
                d.setdefault("stage3_dimsum_equal", []).append(
                    len({x[k]["poly_mul_dimsum"] for k in x}) == 1)
                d.setdefault("success_all", []).append(all(x[k]["success"] for k in x))
                d.setdefault("stage3_dimsum_over_V1fixed_rnmul", []).append(
                    x[("V1", "fixed")]["poly_mul_dimsum"] / x[("V1", "fixed")]["rn_mul"])
            out["fair"][f"{b}/{M}"] = {k: (median(v) if not isinstance(v[0], bool) else all(v))
                                       for k, v in d.items()}
        for Modd in (1, 3, 15, 105, 1155):
            sel = [r for r in rows if r["kind"] == "two_adic" and r["bits"] == b and r["M_odd"] == Modd]
            if not sel:
                continue
            out["two_adic"][f"{b}/{Modd}"] = {
                "V4_over_V0_8": median(r["V4_8"]["work"] / r["V0_8"]["work"] for r in sel),
                "V4_over_V0_2": median(r["V4_8"]["work"] / r["V0_2"]["work"] for r in sel),
                "V0_8_over_V0_2": median(r["V0_8"]["work"] / r["V0_2"]["work"] for r in sel),
                "frac_odd_a": median(r["V4_8"]["odd_a"] / (r["V4_8"]["odd_a"] + r["V4_8"]["even_shift0"] + r["V4_8"]["even_shift32"]) for r in sel),
                "frac_even_shift0": median(r["V4_8"]["even_shift0"] / (r["V4_8"]["odd_a"] + r["V4_8"]["even_shift0"] + r["V4_8"]["even_shift32"]) for r in sel),
            }
        # best over menu
        per = {}
        for r in rows:
            if r["kind"] != "two_adic" or r["bits"] != b:
                continue
            d = per.setdefault(r["seed"], [None, None])
            v0 = min(r["V0_2"]["work"], r["V0_8"]["work"])
            d[0] = v0 if d[0] is None else min(d[0], v0)
            d[1] = r["V4_8"]["work"] if d[1] is None else min(d[1], r["V4_8"]["work"])
        if per:
            out["two_adic"][f"{b}/best_over_menu"] = {
                "median_bestV4_over_bestV0": median(v[1] / v[0] for v in per.values())}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "summarize":
        print(json.dumps(summarize(), indent=1))
    else:
        main()
