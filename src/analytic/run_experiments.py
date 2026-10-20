"""Experiments E1 (Hiary blocks / hull vertices) and E2 (cubic Weyl-sum reduction).

HARNESS. Pre-registration: notes/analytic_prereg.md.  Usage:
    .venv/bin/python -m src.analytic.run_experiments e1
    .venv/bin/python -m src.analytic.run_experiments e2
    .venv/bin/python -m src.analytic.run_experiments summary
Raw JSONL: experiments/raw/analytic_e{1,2}.jsonl; summary: experiments/analysis/analytic_summary.json.
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from fractions import Fraction
from math import isqrt

import gmpy2

from src.common.instances import balanced_semiprime, prime_instance, random_prime, _rng
from src.analytic.hiary import hiary_factor, predicted_blocks, icbrt_ceil
from src.analytic.hull import hull_vertices
from src.analytic.cubic import weyl_oracle, cubic_decode

RAW = "experiments/raw"
ANA = "experiments/analysis"


def _append(path, rec):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(rec) + "\n")


# ------------------------------------------------------------------ E1
E1_HIARY_BITS = (24, 28, 32, 36, 40, 44, 48)
E1_HULL_BITS = (24, 28, 32, 36, 40)
E1_SEEDS = 10
E1_HULL_SEEDS = 5


def run_e1():
    path = f"{RAW}/analytic_e1.jsonl"
    if os.path.exists(path):
        os.remove(path)
    for bits in E1_HIARY_BITS:
        for seed in range(E1_SEEDS):
            for fam in ("prime", "balanced"):
                inst = prime_instance(bits, seed) if fam == "prime" else balanced_semiprime(bits, seed)
                N = inst.N
                st = {}
                t = time.perf_counter()
                f = hiary_factor(N, stats=st)
                wall = time.perf_counter() - t
                ok = (f is None) if fam == "prime" else (f is not None and N % f == 0 and 1 < f < N)
                rec = {"exp": "e1_hiary", "bits": bits, "seed": seed, "family": fam, "N": N,
                       "factor": f, "correct": ok, "stats": st, "pred_blocks": predicted_blocks(N),
                       "wall_s": wall}
                _append(path, rec)
        print("e1 hiary", bits, flush=True)
    for bits in E1_HULL_BITS:
        for seed in range(E1_HULL_SEEDS):
            inst = balanced_semiprime(bits, seed)
            N = inst.N
            lo = icbrt_ceil(N)
            hi = isqrt(N)
            t = time.perf_counter()
            V = hull_vertices(N, lo, hi)
            xs = [v[0] for v in V]
            # per-dyadic-window counts, windows [lo*2^k, lo*2^{k+1})
            wins = []
            X = lo
            while X <= hi:
                Y = min(2 * X - 1, hi)
                wins.append({"X1": X, "X2": Y, "count": sum(1 for x in xs if X <= x <= Y),
                             "full": Y == 2 * X - 1})
                X = 2 * X
            # HARNESS: divisor check uses the labels only to evaluate P1.5
            p, q = inst.labels["factors"]
            divs_in = [d for d in (p,) if lo <= d <= hi]
            divs_are_vertices = all((d, N // d) in set(V) for d in divs_in)
            st = {}
            hiary_factor(N, stats=st, trial_division=False, lo=lo, hi=hi)  # blocks over same range
            rec = {"exp": "e1_hull", "bits": bits, "seed": seed, "N": N, "lo": lo, "hi": hi,
                   "vertices": len(V), "windows": wins, "divisors_checked": len(divs_in),
                   "divisors_are_vertices": divs_are_vertices, "hiary_blocks_same_range": st["blocks"],
                   "wall_s": time.perf_counter() - t}
            _append(path, rec)
        print("e1 hull", bits, flush=True)


# ------------------------------------------------------------------ E2
E2_BITS = (14, 16, 18, 20, 22)
E2_INST = 20
E2_NOISE_INST = 10
E2_DELTAS = (Fraction(1, 100), Fraction(3, 100), Fraction(1, 10), Fraction(3, 10), Fraction(1),
             Fraction(3), Fraction(10))
E2_NOISE_SEEDS = 3


def semiprime_mod3(bits: int, seed: int, residue: int):
    """HARNESS generator: balanced N = pq, p = residue mod 3 and q = 1 mod 3 (residue=1) or any q."""
    params = {"bits": bits, "res": residue}
    rng = _rng("analytic_mod3", params, seed)
    while True:
        p = random_prime(rng, (bits + 1) // 2, top2=True)
        q = random_prime(rng, bits // 2, top2=True)
        if p == q or (p * q).bit_length() != bits:
            continue
        if residue == 1 and (p % 3 != 1 or q % 3 != 1):
            continue
        if residue == 2 and p % 3 != 2:
            continue
        return p * q, sorted((p, q))


def _oracle_table(N):
    cache = {}

    def orc(a):
        a %= N
        if a not in cache:
            cache[a] = Fraction(weyl_oracle(a, N))
        return cache[a]
    return orc, cache


def run_e2():
    path = f"{RAW}/analytic_e2.jsonl"
    if os.path.exists(path):
        os.remove(path)
    for bits in E2_BITS:
        for idx in range(E2_INST):
            for res in (1, 2):
                N, (p, q) = semiprime_mod3(bits, idx, res)
                orc, cache = _oracle_table(N)
                st = {}
                t = time.perf_counter()
                f = cubic_decode(N, orc, random.Random(f"e2|{bits}|{idx}|{res}"), stats=st)
                wall = time.perf_counter() - t
                maxabs = max(abs(float(v)) for v in cache.values()) if (res == 2 and cache) else None
                rec = {"exp": "e2_clean", "bits": bits, "idx": idx, "res": res, "N": N,
                       "factor": f, "success": f is not None and 1 < f < N and N % f == 0,
                       "found_in": st.get("found_in"), "stats": {k: v for k, v in st.items() if k != "found_in"},
                       "control_max_abs_S": maxabs, "wall_s": wall}
                _append(path, rec)
        print("e2 clean", bits, flush=True)
        for idx in range(E2_NOISE_INST):
            N, (p, q) = semiprime_mod3(bits, idx, 1)
            orc, _ = _oracle_table(N)
            for delta in E2_DELTAS:
                for ns in range(E2_NOISE_SEEDS):
                    nrng = random.Random(f"noise|{bits}|{idx}|{delta}|{ns}")
                    noisy = {}

                    def orc_noisy(a, orc=orc, nrng=nrng, noisy=noisy, delta=delta):
                        a %= N
                        if a not in noisy:
                            # uniform rational noise in [-delta, delta] with 2^20 grid points
                            k = nrng.randrange(-(1 << 20), (1 << 20) + 1)
                            noisy[a] = orc(a) + delta * Fraction(k, 1 << 20)
                        return noisy[a]
                    st = {}
                    f = cubic_decode(N, orc_noisy, random.Random(f"e2n|{bits}|{idx}|{ns}"), stats=st)
                    ok = f is not None and 1 < f < N and N % f == 0
                    rec = {"exp": "e2_noise", "bits": bits, "idx": idx, "N": N, "delta": str(delta),
                           "noise_seed": ns, "success": ok, "found_in": st.get("found_in")}
                    _append(path, rec)
        print("e2 noise", bits, flush=True)


# ------------------------------------------------------------------ summary
def _fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / sxx


def summary():
    out = {}
    e1 = [json.loads(l) for l in open(f"{RAW}/analytic_e1.jsonl")]
    H = [r for r in e1 if r["exp"] == "e1_hiary"]
    out["e1_correct"] = f"{sum(r['correct'] for r in H)}/{len(H)}"
    pr = [r for r in H if r["family"] == "prime"]
    ratios = [r["stats"]["blocks"] / r["pred_blocks"] for r in pr]
    out["e1_prime_blocks_over_pred_min_max"] = [min(ratios), max(ratios)]
    by = {}
    for r in pr:
        by.setdefault(r["bits"], []).append(r)
    xs, ys, rows = [], [], []
    for b in sorted(by):
        lnN = sum(math.log(r["N"]) for r in by[b]) / len(by[b])
        blk = sum(r["stats"]["blocks"] for r in by[b]) / len(by[b])
        cf = sum(r["stats"]["cf_steps"] for r in by[b]) / len(by[b])
        xs.append(lnN)
        ys.append(math.log(blk / lnN))
        rows.append({"bits": b, "mean_blocks": blk, "mean_cf_steps": cf,
                     "blocks_over_N13lnN": blk / (math.exp(lnN / 3) * lnN)})
    out["e1_prime_rows"] = rows
    out["e1_slope_blocks_over_lnN"] = _fit(xs, ys)
    U = [r for r in e1 if r["exp"] == "e1_hull"]
    out["e1_divisors_vertices"] = f"{sum(r['divisors_are_vertices'] for r in U)}/{len(U)}"
    byh = {}
    for r in U:
        byh.setdefault(r["bits"], []).append(r)
    xs, ys, rows = [], [], []
    for b in sorted(byh):
        lnN = sum(math.log(r["N"]) for r in byh[b]) / len(byh[b])
        V = sum(r["vertices"] for r in byh[b]) / len(byh[b])
        B = sum(r["hiary_blocks_same_range"] for r in byh[b]) / len(byh[b])
        full = [w["count"] for r in byh[b] for w in r["windows"] if w["full"] and w["X1"] >= 2 * r["lo"]]
        xs.append(lnN)
        ys.append(math.log(V / lnN))
        rows.append({"bits": b, "mean_vertices": V, "V_over_N13lnN": V / (math.exp(lnN / 3) * lnN),
                     "V_over_hiary_blocks": V / B,
                     "window_count_max_over_min": (max(full) / min(full)) if full else None,
                     "window_count_over_N13": (sum(full) / len(full) / math.exp(lnN / 3)) if full else None})
    out["e1_hull_rows"] = rows
    out["e1_slope_vertices_over_lnN"] = _fit(xs, ys)
    e2 = [json.loads(l) for l in open(f"{RAW}/analytic_e2.jsonl")]
    C = [r for r in e2 if r["exp"] == "e2_clean"]
    rows = []
    for b in E2_BITS:
        c1 = [r for r in C if r["bits"] == b and r["res"] == 1]
        c2 = [r for r in C if r["bits"] == b and r["res"] == 2]
        rows.append({"bits": b,
                     "res1_success": sum(r["success"] for r in c1), "res1_structural": sum(r["found_in"] == "structural_zw_gcd" for r in c1),
                     "res1_accidental": sum(r["found_in"] == "gcd_accidental" for r in c1), "res1_n": len(c1),
                     "res2_success": sum(r["success"] for r in c2), "res2_accidental": sum(r["found_in"] == "gcd_accidental" for r in c2),
                     "res2_max_abs_S": max((r["control_max_abs_S"] or 0) for r in c2), "res2_n": len(c2)})
    out["e2_clean_rows"] = rows
    Nn = [r for r in e2 if r["exp"] == "e2_noise"]
    rows, xs, ys = [], [], []
    for b in E2_BITS:
        rb = [r for r in Nn if r["bits"] == b]
        rates = {}
        for d in E2_DELTAS:
            rr = [r for r in rb if r["delta"] == str(d)]
            rates[str(d)] = sum(r["success"] for r in rr) / len(rr)
        ok = [d for d in E2_DELTAS if rates[str(d)] >= 0.5]
        dstar = max(ok) if ok else None
        lnN = sum(math.log(r["N"]) for r in rb) / len(rb)
        rows.append({"bits": b, "rates": rates, "delta_star": str(dstar) if dstar else None})
        if dstar:
            xs.append(lnN)
            ys.append(math.log(float(dstar)))
    out["e2_noise_rows"] = rows
    out["e2_slope_log_deltastar_vs_lnN"] = _fit(xs, ys) if len(xs) >= 2 else None
    os.makedirs(ANA, exist_ok=True)
    with open(f"{ANA}/analytic_summary.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    {"e1": run_e1, "e2": run_e2, "summary": summary}[sys.argv[1]]()
