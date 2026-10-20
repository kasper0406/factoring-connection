"""floorsum experiments E1-E5 (HARNESS).  Pre-registration: notes/floorsum_prereg.md.

    PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments e1|e2|e3|e4|e5|summary
Raw JSONL: experiments/raw/floorsum_e{1..5}.jsonl; summary: experiments/analysis/floorsum_summary.json
Floats are used only in the summary fits.
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

from src.common.instances import balanced_semiprime, prime_instance
from src.analytic.hull import _lower_hull, hull_vertices
from src.analytic.ladder import _Block, iroot_floor
from src.floorsum.farey import Counters, FloorSumConvex, floor_sum_brute, floor_sum_poly
from src.floorsum.ladder_ops import _ladder_coef, ladder_factor_op
from src.floorsum.lattice_block import best_params, lattice_factor

RAW = "experiments/raw"
ANA = "experiments/analysis"


def _append(path, rec):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(rec) + "\n")


def _fresh(path):
    if os.path.exists(path):
        os.remove(path)


def _ceil_div(a, b):
    return -((-a) // b)


def _blocks(N, j, c=Fraction(1, 2)):
    sN = isqrt(N)
    B = iroot_floor(N, j + 2)
    x0 = min(sN + 1, -((-B * c.denominator) // c.numerator))
    x = x0
    while x <= sN:
        H = max(1, (c.numerator * x) // (c.denominator * B))
        H = min(H, sN - x + 1)
        yield x, H
        x += H


# ------------------------------------------------------------------ E1
E1_BITS = (28, 32, 36, 40, 44)
E1_SEEDS = 5
E1_J = (2, 3, 4)


def run_e1():
    path = f"{RAW}/floorsum_e1.jsonl"
    _fresh(path)
    for bits in E1_BITS:
        for seed in range(E1_SEEDS):
            N = balanced_semiprime(bits, seed).N
            for j in E1_J:
                t = time.perf_counter()
                Vtot = nblk = kc = 0
                x0 = None
                xs_last = None
                for x, H in _blocks(N, j):
                    if x0 is None:
                        x0 = x
                    blk = _Block(N, x, H, j)
                    den = blk.den
                    pts = []
                    rho_num, rho_den = blk.rho.numerator, blk.rho.denominator
                    for h in range(H):
                        v = blk.num(h)
                        pts.append((h, _ceil_div(v, den)))
                        r = v % den  # {P} = r/den
                        if j % 2 == 0:
                            if r * rho_den <= rho_num * den:
                                kc += 1
                        else:
                            if (den - r) * rho_den <= rho_num * den and r != 0:
                                kc += 1
                    Vtot += len(_lower_hull(pts))
                    nblk += 1
                    xs_last = x + H - 1
                if x0 is None:
                    continue
                Vg = len(hull_vertices(N, x0, xs_last))
                B = iroot_floor(N, j + 2)
                rec = {"exp": "e1", "bits": bits, "seed": seed, "j": j, "N": N, "x0": x0,
                       "x_end": xs_last, "B": B, "blocks": nblk, "V_tot": Vtot, "V_glob": Vg,
                       "K_cand": kc, "wall_s": time.perf_counter() - t}
                _append(path, rec)
            print("e1", bits, seed, flush=True)


# ------------------------------------------------------------------ E2
E2_SCALE = {2: (48, 64, 80, 96, 112, 128), 3: (48, 64, 80, 96, 112), 4: (48, 64, 80, 96)}


def run_e2():
    path = f"{RAW}/floorsum_e2.jsonl"
    _fresh(path)
    rng = random.Random("floorsum-e2")
    # exactness A: random convex polynomials
    bad = 0
    for t in range(3000):
        deg = rng.choice([2, 3, 4])
        den = rng.randint(1, 10 ** 15)
        H = rng.randint(1, 3000)
        coef = [rng.randint(-10 ** 15, 10 ** 15), rng.randint(-10 ** 15, 10 ** 15),
                rng.randint(1, 10 ** rng.randint(1, 12))] + [0] * (deg - 2)
        if deg >= 3:
            coef[3] = rng.randint(0, 10 ** 3)
        if deg == 4:
            coef[4] = rng.randint(0, 10)
        a = rng.randint(-100, 100)
        if floor_sum_poly(coef, den, a, a + H) != floor_sum_brute(coef, den, a, a + H):
            bad += 1
    _append(path, {"exp": "e2_exact_random", "cases": 3000, "mismatches": bad})
    print("e2 random exactness", bad, flush=True)
    # exactness B: ladder blocks
    for j in (2, 3, 4):
        for bits in (40, 48, 56):
            bad = tot = 0
            for t in range(50):
                N = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
                B = iroot_floor(N, j + 2)
                x = rng.randint(2 * B, isqrt(N))
                H = max(1, x // (2 * B))
                blk = _Block(N, x, H, j)
                eta = Fraction(1, 4 * blk.den * blk.rho.denominator)
                for sh in (Fraction(0), -blk.rho - eta, blk.rho):
                    c = _ladder_coef(blk, sh.numerator, sh.denominator)
                    d = blk.den * sh.denominator
                    F = FloorSumConvex(c, d, 0, H)
                    for _ in range(4):
                        h1 = rng.randint(0, H)
                        h2 = rng.randint(h1, H)
                        tot += 1
                        if F.range_sum(h1, h2) != floor_sum_brute(c, d, h1, h2):
                            bad += 1
            _append(path, {"exp": "e2_exact_ladder", "j": j, "bits": bits, "queries": tot,
                           "mismatches": bad})
            print("e2 ladder exactness", j, bits, bad, tot, flush=True)
    # scaling
    for j, sizes in E2_SCALE.items():
        for bits in sizes:
            for t in range(5):
                N = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
                B = iroot_floor(N, j + 2)
                x = isqrt(N) // 2 + rng.randint(0, 1000)
                H = max(1, x // (2 * B))
                blk = _Block(N, x, H, j)
                cnt = Counters()
                t0 = time.perf_counter()
                F = FloorSumConvex(_ladder_coef(blk, 0, 1), blk.den, 0, H, cnt)
                wall = time.perf_counter() - t0
                rec = {"exp": "e2_scale", "j": j, "bits": bits, "N": N, "x": x, "H": H, "B": B,
                       "Q": getattr(F, "Q", None), "windows": cnt.get("windows", 0),
                       "runs": cnt.get("runs", 0), "evals": cnt.get("evals", 0),
                       "bsearch": cnt.get("bsearch", 0), "wall_s": wall}
                _append(path, rec)
            print("e2 scale", j, bits, flush=True)


# ------------------------------------------------------------------ E3
E3_BITS = (28, 32, 36, 40, 44)


def run_e3():
    path = f"{RAW}/floorsum_e3.jsonl"
    _fresh(path)
    for bits in E3_BITS:
        for j in (2, 3):
            for fam, seeds in (("prime", range(3)), ("balanced", range(5))):
                for seed in seeds:
                    inst = prime_instance(bits, seed) if fam == "prime" else balanced_semiprime(bits, seed)
                    N = inst.N
                    st = {}
                    t = time.perf_counter()
                    f = ladder_factor_op(N, j, stats=st)
                    wall = time.perf_counter() - t
                    ok = (f is None) if fam == "prime" else (f is not None and 1 < f < N and N % f == 0)
                    _append(path, {"exp": "e3", "bits": bits, "j": j, "family": fam, "seed": seed,
                                   "N": N, "factor": f, "correct": ok, "stats": st, "wall_s": wall})
            print("e3", bits, j, flush=True)


# ------------------------------------------------------------------ E4
E4_BITS = (48, 56, 64)


def run_e4():
    path = f"{RAW}/floorsum_e4.jsonl"
    _fresh(path)
    rng = random.Random("floorsum-e4")
    for bits in E4_BITS:
        for t in range(200):
            N = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            j = 2
            B = iroot_floor(N, j + 2)
            sN = isqrt(N)
            x = rng.randint(sN // 4, sN)
            H = max(1, x // (2 * B))
            blk = _Block(N, x, H, j)
            den = blk.den
            K0 = max(1, math.isqrt(H - 1) + 1)  # ceil(sqrt H)
            n_d = 0
            imax = max(0, H.bit_length())
            hist = [0] * (imax + 2)  # hist[i] = #h counted for K = 2^i (cumulative later)
            for h in range(H):
                r = blk.num(h) % den
                dist = min(r, den - r)
                if K0 * dist < den:
                    n_d += 1
                # largest i with 2^i * dist < den
                if dist == 0:
                    hist[imax + 1] += 1
                else:
                    lim = (den - 1) // dist  # K <= lim
                    i = min(imax + 1, lim.bit_length() - 1) if lim >= 1 else -1
                    if i >= 0:
                        hist[i] += 1
            # n_{1/K} for K = 2^i is the number of h with top index >= i
            cum = [0] * (imax + 2)
            s = 0
            for i in range(imax + 1, -1, -1):
                s += hist[i]
                cum[i] = s
            M = min((1 << i) + cum[i] for i in range(0, imax + 1))
            _append(path, {"exp": "e4", "bits": bits, "N": N, "x": x, "H": H, "K": K0,
                           "n_delta": n_d, "ratio": n_d / (2 * H / K0), "M": M,
                           "M_over_sqrtH": M / math.sqrt(H)})
        print("e4", bits, flush=True)


# ------------------------------------------------------------------ E5
E5_BITS = (36, 40, 44, 48)
E5_SEEDS = 3


def run_e5():
    path = f"{RAW}/floorsum_e5.jsonl"
    _fresh(path)
    for bits in E5_BITS:
        for seed in range(E5_SEEDS):
            inst = balanced_semiprime(bits, seed)
            N = inst.N
            st = {}
            t = time.perf_counter()
            f = lattice_factor(N, stats=st)
            wall = time.perf_counter() - t
            sN = isqrt(N)
            m, n, X = best_params(N, sN)
            Hlad = sN // (2 * iroot_floor(N, 4))
            _append(path, {"exp": "e5", "bits": bits, "seed": seed, "N": N, "factor": f,
                           "correct": f is not None and 1 < f < N and N % f == 0, "stats": st,
                           "X_at_sqrtN": X, "m": m, "n": n, "H_ladder_j2": Hlad,
                           "wall_s": wall})
        print("e5", bits, flush=True)


# ------------------------------------------------------------------ summary
def _load(name):
    p = f"{RAW}/floorsum_{name}.jsonl"
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


def _slope(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)


def summary():
    out = {}
    e1 = _load("e1")
    if e1:
        s = {}
        for j in E1_J:
            rows = [r for r in e1 if r["j"] == j]
            p11 = all(r["V_tot"] >= r["V_glob"] - r["K_cand"] for r in rows)
            ratios = [r["V_tot"] / r["V_glob"] for r in rows]
            xs, ys, per = [], [], {}
            for b in E1_BITS:
                rb = [r for r in rows if r["bits"] == b]
                if not rb:
                    continue
                lnN = sum(math.log(r["N"]) for r in rb) / len(rb)
                mV = sum(r["V_tot"] for r in rb) / len(rb)
                xs.append(lnN)
                ys.append(math.log(mV / lnN))
                per[b] = {"mean_V_tot": mV, "mean_V_glob": sum(r["V_glob"] for r in rb) / len(rb),
                          "mean_blocks": sum(r["blocks"] for r in rb) / len(rb),
                          "per_block_over_scale": sum(r["V_tot"] / r["blocks"] /
                                                      ((2 * r["N"]) ** (1 / 3) / (2 * r["B"]))
                                                      for r in rb) / len(rb),
                          "V_tot_over_N13lnN": mV / (math.exp(lnN / 3) * lnN)}
            s[j] = {"P1.1_holds": p11, "ratio_min": min(ratios), "ratio_max": max(ratios),
                    "slope": _slope(xs, ys) if len(xs) > 1 else None, "per_bits": per}
        out["e1"] = s
    e2 = _load("e2")
    if e2:
        s = {"exact": [r for r in e2 if r["exp"].startswith("e2_exact")]}
        for j in E2_SCALE:
            rows = [r for r in e2 if r["exp"] == "e2_scale" and r["j"] == j]
            if not rows:
                continue
            ratios = []
            xs, ys = [], []
            for b in E2_SCALE[j]:
                rb = [r for r in rows if r["bits"] == b]
                if not rb:
                    continue
                for r in rb:
                    beta = r["N"] / r["x"] ** 3
                    ratios.append((r["windows"] + r["runs"]) / (r["H"] * beta ** (1 / 3)))
                xs.append(sum(math.log(r["H"]) for r in rb) / len(rb))
                ys.append(math.log(sum(r["windows"] + r["runs"] for r in rb) / len(rb)))
            s[j] = {"ratio_min": min(ratios), "ratio_max": max(ratios),
                    "slope_vs_lnH": _slope(xs, ys), "pred": 2 * (j - 1) / (3 * j)}
        out["e2"] = s
    e3 = _load("e3")
    if e3:
        s = {"correct": sum(r["correct"] for r in e3), "total": len(e3)}
        for j in (2, 3):
            rows = [r for r in e3 if r["j"] == j and r["family"] == "prime"]
            xs, ys, per = [], [], {}
            for b in E3_BITS:
                rb = [r for r in rows if r["bits"] == b]
                if not rb:
                    continue
                lnN = sum(math.log(r["N"]) for r in rb) / len(rb)
                me = sum(r["stats"].get("op_evals", 0) for r in rb) / len(rb)
                xs.append(lnN)
                ys.append(math.log(me / lnN ** 2))
                per[b] = {"mean_evals": me,
                          "mean_blocks": sum(r["stats"]["blocks"] for r in rb) / len(rb),
                          "mean_calls": sum(r["stats"]["oracle_calls"] for r in rb) / len(rb),
                          "mean_windows": sum(r["stats"].get("op_windows", 0) for r in rb) / len(rb),
                          "mean_runs": sum(r["stats"].get("op_runs", 0) for r in rb) / len(rb),
                          "mean_wall": sum(r["wall_s"] for r in rb) / len(rb)}
            s[j] = {"slope": _slope(xs, ys) if len(xs) > 1 else None, "per_bits": per}
        out["e3"] = s
    e4 = _load("e4")
    if e4:
        s = {}
        for b in E4_BITS:
            rb = [r for r in e4 if r["bits"] == b]
            if rb:
                s[b] = {"mean_ratio": sum(r["ratio"] for r in rb) / len(rb),
                        "mean_M_over_sqrtH": sum(r["M_over_sqrtH"] for r in rb) / len(rb),
                        "mean_H": sum(r["H"] for r in rb) / len(rb)}
        out["e4"] = s
    e5 = _load("e5")
    if e5:
        out["e5"] = {"correct": sum(r["correct"] for r in e5), "total": len(e5),
                     "rows": [{k: r[k] for k in ("bits", "seed", "X_at_sqrtN", "H_ladder_j2", "m", "n",
                                                 "wall_s")} | {"blocks": r["stats"].get("blocks")}
                              for r in e5]}
    os.makedirs(ANA, exist_ok=True)
    with open(f"{ANA}/floorsum_summary.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    {"e1": run_e1, "e2": run_e2, "e3": run_e3, "e4": run_e4, "e5": run_e5,
     "summary": summary}[sys.argv[1]]()
