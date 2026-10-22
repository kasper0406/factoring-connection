"""Pre-registered experiments (notes/sinonobliv_prereg.md).

  .venv/bin/python -m src.sinonobliv.run_exp E1 | E2A | E2B | all

Raw output: experiments/raw/sinonobliv_E{1,2A,2B}.jsonl. Exact integer smoothness decisions throughout; floats
only in predictions and ratios.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

from .core import (DensityCurve, count_members_in_window, divisors_from_factored, family_pieces,
                   family_t_for_size, family_value, lcm_upto_factored, smooth_mask_values,
                   smooth_mask_window, value_of_factored, ch_log_threshold, johnson_log_threshold)
from .msample import (divisor_trick_yield, inheritance_scan, joint_tuples_indep, predicted_smooth_first,
                      predicted_tuples)
from .theory import consecutive_product_log_hits, toy_path

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, "sinonobliv_%s.jsonl" % name), "w")


# ------------------------------------------------------------------------------------------------ E1

E1_CELLS = {"cp2": (40, 48, 56), "cp3": (48, 56, 60), "aur": (40, 46, 52), "sq1": (40, 48, 56)}
E1_Y = (1 << 8, 1 << 10)


def _t_range(fam, bits):
    t0 = family_t_for_size(fam, 1 << bits)
    if fam in ("cp2", "sq1"):
        T = min(t0 // 16, 1 << 21)
    elif fam == "cp3":
        T = t0
    else:  # aur
        T = 3 * t0
    return t0, T


def _piece_masks(fam, t0, T, y):
    if fam in ("cp2", "cp3"):
        r = 2 if fam == "cp2" else 3
        m = smooth_mask_window(t0, T + r - 1, y)
        return [m[j:j + T] for j in range(r)]
    t = np.arange(t0, t0 + T, dtype=object)
    pieces = family_pieces(fam, t)
    return [smooth_mask_values(np.array([int(v) for v in p], dtype=np.uint64), y) for p in pieces]


def run_e1():
    f = _out("E1")
    for fam, bitss in E1_CELLS.items():
        for bits in bitss:
            t0, T = _t_range(fam, bits)
            tt = np.arange(t0, t0 + T, dtype=object)
            vals = family_value(fam, tt)
            vmin, vmax = int(vals[0]), int(vals[-1])
            for y in E1_Y:
                t_start = time.time()
                masks = _piece_masks(fam, t0, T, y)
                joint = np.ones(T, dtype=bool)
                for mk in masks:
                    joint &= mk
                nb = 16
                edges = np.linspace(0, T, nb + 1).astype(int)
                pred = 0.0
                for b in range(nb):
                    a, c = edges[b], edges[b + 1]
                    if c <= a:
                        continue
                    prod = 1.0
                    for mk in masks:
                        prod *= float(mk[a:c].mean())
                    pred += (c - a) * prod
                obs = int(joint.sum())
                # control density near member sizes
                b_lo = math.log2(vmin) - 0.25
                b_hi = min(math.log2(vmax) + 0.25, 62.4)
                curve = DensityCurve(y, b_lo, b_hi, width=1 << 21)
                sample = np.linspace(0, T - 1, min(T, 4096)).astype(int)
                ctrl = float(np.mean(curve(np.array([float(vals[i]) for i in sample])))) * T
                rec = {"family": fam, "bits": bits, "y": y, "t0": t0, "T": T, "obs": obs, "pred_indep": pred,
                       "R": obs / pred if pred > 0 else None, "control_expected": ctrl,
                       "gain_vs_random": obs / ctrl if ctrl > 0 else None,
                       "struct_gain": pred / ctrl if ctrl > 0 else None,
                       "piece_rates": [float(mk.mean()) for mk in masks],
                       "joint_rate": obs / T, "curve_zero_cells": curve.zero_cells,
                       "secs": round(time.time() - t_start, 2)}
                print(json.dumps(rec))
                f.write(json.dumps(rec) + "\n")
                f.flush()
    # P2: hit law in windows I_k (cp2, cp3)
    rng = np.random.default_rng(2022)
    for fam in ("cp2", "cp3"):
        for bits in (40, 48):
            y, h = toy_path(bits)
            obs = 0
            pred = 0.0
            for _ in range(20):
                N = int(rng.integers(1 << (bits - 1), 1 << bits)) | 1
                for k in range(1, 65):
                    obs += count_members_in_window(fam, k * N, k * N + h)
                    t = family_t_for_size(fam, k * N)
                    # F'(t): derivative of the product (exact integer), continuous prediction (h+1)/F'(t)
                    if fam == "cp2":
                        dF = 2 * t + 1
                    else:
                        dF = 3 * t * t + 6 * t + 2
                    pred += (h + 1) / dF
            rec = {"P2": True, "family": fam, "bits": bits, "h": h, "obs": obs, "pred": pred,
                   "ratio": obs / pred if pred else None}
            print(json.dumps(rec))
            f.write(json.dumps(rec) + "\n")
    # D1 model extrapolation
    for bits in (512, 1024, 2048, 4096):
        x = bits * math.log(2)
        rec = {"D1_model": True, "bits": bits,
               "ln_smooth_members_per_window": {r: consecutive_product_log_hits(x, r) for r in (2, 3, 4, 8)}}
        f.write(json.dumps(rec) + "\n")
    f.close()


# ------------------------------------------------------------------------------------------------ E2

def run_e2a():
    f = _out("E2A")
    rng = np.random.default_rng(3022)
    for bits in (32, 40):
        for i in range(10):
            N = int(rng.integers(1 << (bits - 1), 1 << bits)) | 1
            t = time.time()
            res = inheritance_scan(N, list(range(1, 9)), 1 << 9)
            rec = {"bits": bits, "N": N, "K": 8, "h": 1 << 9, **res, "secs": round(time.time() - t, 2)}
            print(json.dumps(rec))
            f.write(json.dumps(rec) + "\n")
            f.flush()
    f.close()


E2B_CONFIGS = [
    {"name": "C1", "Xb": 36, "hb": 10, "B": 32, "y": 1 << 8, "calls": 20000, "bare_calls": 2000},
    {"name": "C2", "Xb": 44, "hb": 12, "B": 41, "y": 1 << 10, "calls": 40000, "bare_calls": 2000},
]


def run_e2b():
    f = _out("E2B")
    for cfg in E2B_CONFIGS:
        X, h, B, y = 1 << cfg["Xb"], 1 << cfg["hb"], cfg["B"], cfg["y"]
        fac = lcm_upto_factored(B)
        P = value_of_factored(fac)
        divs = divisors_from_factored(fac)
        t, pi = math.log(h + 1), math.log(P)
        curve = DensityCurve(y, 2.0, cfg["Xb"] + 1.0, width=1 << 18)
        cells = [("J", 1, johnson_log_threshold(t, pi), cfg["calls"]),
                 ("CH", 2, ch_log_threshold(t, pi, 2), cfg["calls"]),
                 ("CH", 3, ch_log_threshold(t, pi, 3), cfg["calls"]),
                 ("bare", 2, math.log(2 * (h + 1)), cfg["bare_calls"]),
                 ("bare", 3, math.log(2 * (h + 1)), cfg["bare_calls"])]
        for idx, (kind, m, sig, calls) in enumerate(cells):
            rng = np.random.default_rng(4022 + 100 * idx + (0 if cfg["name"] == "C1" else 7))
            tt = time.time()
            res = joint_tuples_indep(X, h, P, B, m, sig, y, calls, rng)
            pt = calls * predicted_tuples(divs, h, m, sig)
            ps = calls * predicted_smooth_first(divs, h, m, sig, X, curve)
            rec = {"config": cfg["name"], "kind": kind, "m": m, "sigma": sig, "d_min": math.exp(sig), **res,
                   "pred_tuples": pt, "ratio_tuples": res["tuples"] / pt if pt else None,
                   "pred_z1_smooth": ps, "ratio_z1_smooth": res["tuples_z1_smooth"] / ps if ps else None,
                   "smooth_yield_per_window": res["distinct_smooth_z1"] / (calls * m),
                   "secs": round(time.time() - tt, 1)}
            print(json.dumps(rec))
            f.write(json.dumps(rec) + "\n")
            f.flush()
        rng = np.random.default_rng(5022 + (0 if cfg["name"] == "C1" else 7))
        res = divisor_trick_yield(X, h, y, 20000, rng)
        rec = {"config": cfg["name"], "kind": "divisor_trick", **res,
               "smooth_yield_per_window": res["successes"] / res["calls"]}
        print(json.dumps(rec))
        f.write(json.dumps(rec) + "\n")
    f.close()


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("E1", "all"):
        run_e1()
    if what in ("E2A", "all"):
        run_e2a()
    if what in ("E2B", "all"):
        run_e2b()
