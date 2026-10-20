"""topblock experiments (pre-registered in notes/topblock_prereg.md).  HARNESS.

  PYTHONPATH=. .venv/bin/python -m src.topblock.run_experiments {x1,x2,x3,x4,x5,x6,summary}

Raw JSONL: experiments/raw/topblock_<exp>.jsonl.  Summary: experiments/analysis/topblock_summary.json.
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from fractions import Fraction

from src.common.instances import balanced_semiprime
from src.topblock import core
from src.topblock import lattices as L

RAW = "experiments/raw"
ANA = "experiments/analysis"


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, f"topblock_{name}.jsonl"), "w")


def _rng(*key):
    return random.Random("|".join(str(k) for k in ("topblock",) + key))


# ------------------------------------------------------------------ X1: (W) characterisation

def x1():
    f = _out("x1")
    mism = tested = 0
    for bits in (20, 24, 28):
        for seed in range(10):
            inst = balanced_semiprime(bits, seed)
            N, p = inst.N, inst.factors[0]
            for j in (2, 3, 4):
                rng = _rng("x1", bits, seed, j)
                blocks = [core.top_block_containing(N, p, j)]
                B = core.ladder_B(N, j)
                for _ in range(3):
                    x = rng.randrange(math.isqrt(N) // 2, math.isqrt(N))
                    blocks.append((x, max(1, x // (2 * B))))
                for (x, H) in blocks:
                    rho = Fraction(N * H ** (j + 1), x ** (j + 2))
                    wins = [(Fraction(0), rho) if j % 2 == 0 else (1 - rho, Fraction(1) - Fraction(1, 10 ** 9))]
                    for _ in range(10):
                        a = Fraction(rng.randrange(0, 1000), 1000)
                        b = min(a + Fraction(rng.randrange(1, 200), 1000), Fraction(999, 1000))
                        wins.append((a, b))
                    for (a, b) in wins:
                        s1 = core.window_set_direct(N, x, j, H, a, b)
                        s2 = core.window_set_neardiv(N, x, j, H, a, b)
                        tested += 1
                        ok = s1 == s2
                        mism += not ok
                        f.write(json.dumps({"bits": bits, "seed": seed, "j": j, "x": x, "H": H,
                                            "a": str(a), "b": str(b), "n_window": len(s1), "match": ok}) + "\n")
    f.close()
    print(f"x1: {tested} (block, window) pairs, {mism} mismatches")


# ------------------------------------------------------------------ X2: reductions (iii), (iv)

def x2():
    f = _out("x2")
    for (bits, js) in ((40, (2, 3, 4)), (48, (2, 3)), (44, (4,))):
        for seed in range(10):
            inst = balanced_semiprime(bits, seed)
            N, p = inst.N, inst.factors[0]
            for j in js:
                x, H = core.top_block_containing(N, p, j)
                orc = core.WindowOracle(N, x, j, H)
                divs, calls = core.locate_divisor_via_window_oracle(N, x, j, H, orc.query)
                Cw = core.weighted_neardiv_bound(N, x, j, H)
                rec = {"exp": "iii", "bits": bits, "seed": seed, "j": j, "H": H, "found_p": p in divs,
                       "divs": divs, "calls": calls, "log2H": math.log2(H), "Cw_bound": Cw,
                       "calls_le_bound": calls <= Cw + 1e-9}
                f.write(json.dumps(rec) + "\n")
                print(rec["bits"], j, H, rec["found_p"], calls, round(Cw, 1), flush=True)
    # (iv) at 32 bits
    mism = tested = 0
    for seed in range(10):
        inst = balanced_semiprime(32, seed)
        N, p = inst.N, inst.factors[0]
        for j in (2, 3):
            x, H = core.top_block_containing(N, p, j)
            orc = core.WindowOracle(N, x, j, H)
            rng = _rng("x2iv", seed, j)
            for _ in range(20):
                h1 = rng.randrange(0, H - 1)
                h2 = rng.randrange(h1 + 1, H + 1)
                a = Fraction(rng.randrange(0, 1000), 1000)
                b = min(a + Fraction(rng.randrange(1, 30), 1000), Fraction(999, 1000))
                c1 = orc.query(h1, h2, a, b)
                c2, fc = core.window_count_via_divisor_finder(N, x, j, h1, h2, a, b)
                tested += 1
                mism += c1 != c2
                f.write(json.dumps({"exp": "iv", "seed": seed, "j": j, "h1": h1, "h2": h2, "a": str(a),
                                    "b": str(b), "count": c1, "count_iv": c2, "finder_calls": fc}) + "\n")
    f.close()
    print(f"x2 (iv): {tested} queries, {mism} mismatches")


# ------------------------------------------------------------------ X3/X4/X5/X6 lattice scans

U_LIST = [(3, 8), (4, 10), (5, 12), (6, 14), (7, 16), (9, 20), (11, 24), (13, 28), (15, 32)]


def _h0(bits, seed, lb2):
    """Planted offset for scan point log2 X = lb2/2 (identical for every lattice at that point)."""
    X = int(round(2 ** (lb2 / 2)))
    rng = _rng("h0", bits, seed, lb2)
    return X, rng.randrange(max(1, X // 2), X)


def scan(test, bits, seed, p, lb2_start=8, lb2_stop=None, extend_good=False):
    """Scan log2 X = lb2/2 upward.  Xmax = last success before two consecutive failures.
    If extend_good, keep scanning (1-bit steps) until lb2_stop recording HG-good counts."""
    recs = []
    lb2 = lb2_start
    fails = 0
    last = None
    while True:
        X, h0 = _h0(bits, seed, lb2)
        r = test(X, p - h0)
        r.update({"log2X": lb2 / 2})
        recs.append(r)
        if r["success"]:
            last = lb2 / 2
            fails = 0
        else:
            fails += 1
            if fails >= 2:
                break
        lb2 += 1
        if lb2_stop and lb2 > lb2_stop:
            break
    if extend_good and lb2_stop:
        lb2 += 2
        while lb2 <= lb2_stop:
            X, h0 = _h0(bits, seed, lb2)
            r = test(X, p - h0)
            r.update({"log2X": lb2 / 2, "extension": True})
            recs.append(r)
            lb2 += 2
    return last, recs


def _univ_table(N, p, bits, seed, f):
    tab = {}
    for (m, n) in U_LIST:
        xm, recs = scan(lambda X, x: L.univariate_attempt(N, x, m, n, X), bits, seed, p)
        tab[(m, n)] = xm
        f.write(json.dumps({"kind": "U", "bits": bits, "seed": seed, "m": m, "n": n, "Xmax": xm}) + "\n")
    return tab


def _ucomp(N, p, bits, seed, tab, M, nimg, f):
    key = (M, nimg)
    if key not in tab:
        xm, _ = scan(lambda X, x: L.univariate_attempt(N, x, M, nimg, X), bits, seed, p)
        tab[key] = xm
        f.write(json.dumps({"kind": "U", "bits": bits, "seed": seed, "m": M, "n": nimg, "Xmax": xm}) + "\n")
    return tab[key]


def _best_u_dim(tab, nmax):
    vals = [v for (m, n), v in tab.items() if n <= nmax and v is not None]
    return max(vals) if vals else None


X3_CONFIGS = {2: [(1, 1, 8), (2, 1, 10), (2, 2, 10)], 3: [(1, 1, 9), (2, 1, 11), (2, 2, 12)]}


def x3():
    f = _out("x3")
    for bits in (64, 96, 128):
        for seed in range(3):
            inst = balanced_semiprime(bits, seed)
            N, p = inst.N, inst.factors[0]
            tab = _univ_table(N, p, bits, seed, f)
            for j, cfgs in X3_CONFIGS.items():
                for (m1, m2, E) in cfgs:
                    win = [(0, j, m2)]
                    nmix = L.mixed_dim(win, E)
                    stop = int(2 * 0.45 * bits) if m1 == 1 else None
                    xm, recs = scan(lambda X, x: L.mixed_attempt(N, x, win, m1, E, X), bits, seed, p,
                                    lb2_stop=stop, extend_good=(m1 == 1))
                    M, nimg = m1 + m2, E + m2 + 1
                    ui = _ucomp(N, p, bits, seed, tab, M, nimg, f)
                    ud = _best_u_dim(tab, nmix)
                    goodX = [r["log2X"] for r in recs if r["good"] > 0]
                    nontriv_beyond = sum(r["good_nontrivial"] for r in recs
                                         if ui is not None and r["log2X"] > (ui + 3))
                    good_beyond = [r["log2X"] for r in recs if r["good"] > 0 and ui is not None and r["log2X"] >= ui + 4]
                    Fshort = [r["log2X"] for r in recs if L.F_relation_row_is_short(
                        N, p - _h0(bits, seed, int(2 * r["log2X"]))[1], j, int(round(2 ** r["log2X"])), m1, m2)]
                    rec = {"kind": "MIX", "bits": bits, "seed": seed, "j": j, "m1": m1, "m2": m2, "E": E,
                           "n_mix": nmix, "Xmax": xm, "U_img": [M, nimg], "Xmax_U_img": ui, "Xmax_U_dim": ud,
                           "max_log2X_good": max(goodX) if goodX else None,
                           "max_log2X_Fshort": max(Fshort) if Fshort else None,
                           "good_beyond_Uimg_plus4": good_beyond, "nontrivial_beyond_Uimg_plus3": nontriv_beyond,
                           "scan": recs}
                    f.write(json.dumps(rec) + "\n")
                    print(bits, seed, j, (m1, m2, E), nmix, "MIX", xm, "Uimg", ui, "Udim", ud,
                          "good<=", rec["max_log2X_good"], "F<=", rec["max_log2X_Fshort"], flush=True)
    f.close()


X4_CONFIGS = [  # (label, windows-as-function-of-X, m1, E)
    ("adjacent_j2", lambda X: [(0, 2, 1), (X, 2, 1)], 1, 7),
    ("adjacent_j2", lambda X: [(0, 2, 1), (X, 2, 1)], 2, 9),
    ("degrees_23", lambda X: [(0, 2, 1), (0, 3, 1)], 1, 8),
    ("degrees_23", lambda X: [(0, 2, 1), (0, 3, 1)], 2, 10),
]


def x4():
    f = _out("x4")
    for bits in (64, 96):
        for seed in range(3):
            inst = balanced_semiprime(bits, seed)
            N, p = inst.N, inst.factors[0]
            tab = _univ_table(N, p, bits, seed, f)
            for (label, wf, m1, E) in X4_CONFIGS:
                nmix = L.mixed_dim(wf(1), E)
                xm, recs = scan(lambda X, x: L.mixed_attempt(N, x, wf(X), m1, E, X), bits, seed, p)
                T = sum(w[2] for w in wf(1))
                ui = _ucomp(N, p, bits, seed, tab, m1 + T, E + T + 1, f)
                ud = _best_u_dim(tab, nmix)
                rec = {"kind": "MIX2", "label": label, "bits": bits, "seed": seed, "m1": m1, "E": E, "n_mix": nmix,
                       "Xmax": xm, "U_img": [m1 + T, E + T + 1], "Xmax_U_img": ui, "Xmax_U_dim": ud, "scan": recs}
                f.write(json.dumps(rec) + "\n")
                print(bits, seed, label, m1, E, nmix, "MIX2", xm, "Uimg", ui, "Udim", ud, flush=True)
    f.close()


X5_CONFIGS = [(2, 4), (3, 6)]


def x5():
    f = _out("x5")
    for bits in (64, 96):
        for seed in range(3):
            inst = balanced_semiprime(bits, seed)
            N, p = inst.N, inst.factors[0]
            tab = _univ_table(N, p, bits, seed, f)
            for (m, E) in X5_CONFIGS:
                def test(X, x):
                    lb = max(1, int(math.log2(X)))
                    Hs = 2 ** (lb // 2)
                    U = -(-X // Hs)
                    return L.hm_attempt(N, x, Hs, U, m, E)
                nhm = (E + 1) * (E + 2) // 2
                xm, recs = scan(test, bits, seed, p)
                ud = _best_u_dim(tab, nhm)
                good = sum(r["good"] for r in recs)
                line = sum(r["good_line"] for r in recs)
                rec = {"kind": "HM", "bits": bits, "seed": seed, "m": m, "E": E, "n": nhm, "Xmax": xm,
                       "Xmax_U_dim": ud, "good_total": good, "good_line_total": line, "scan": recs}
                f.write(json.dumps(rec) + "\n")
                print(bits, seed, (m, E), nhm, "HM", xm, "Udim", ud, "line", line, "/", good, flush=True)
    f.close()


X6_CONFIGS = [(3, 8), (4, 10), (5, 12), (6, 14), (7, 16)]


def x6():
    f = _out("x6")
    for bits in (64, 96):
        for seed in range(3):
            inst = balanced_semiprime(bits, seed)
            N, p = inst.N, inst.factors[0]
            for (m, n) in X6_CONFIGS:
                recs = []
                lb2 = 8
                fails = 0
                xl = xs = None
                while True:
                    X, h0 = _h0(bits, seed, lb2)
                    r = L.univariate_attempt(N, p - h0, m, n, X, use_svp=True)
                    r["log2X"] = lb2 / 2
                    recs.append(r)
                    if r["success"]:
                        xl = lb2 / 2
                    if r["svp_success_hg"]:
                        xs = lb2 / 2
                    if not r["success"] and not r["svp_success_hg"]:
                        fails += 1
                        if fails >= 2:
                            break
                    else:
                        fails = 0
                    lb2 += 1
                rec = {"kind": "SVP", "bits": bits, "seed": seed, "m": m, "n": n, "Xmax_LLL": xl, "Xmax_SVP": xs,
                       "scan": recs}
                f.write(json.dumps(rec) + "\n")
                print(bits, seed, (m, n), "LLL", xl, "SVP", xs, flush=True)
    f.close()


# ------------------------------------------------------------------ summary

def _load(name):
    path = os.path.join(RAW, f"topblock_{name}.jsonl")
    if not os.path.exists(path):
        return []
    return [json.loads(l) for l in open(path)]


def summary():
    out = {}
    x1r = _load("x1")
    out["x1"] = {"pairs": len(x1r), "mismatches": sum(not r["match"] for r in x1r)}
    x2r = _load("x2")
    iii = [r for r in x2r if r["exp"] == "iii"]
    iv = [r for r in x2r if r["exp"] == "iv"]
    if iii:
        ratios = sorted(r["calls"] / max(1, math.ceil(r["log2H"])) for r in iii)
        out["x2"] = {"blocks": len(iii), "found_p": sum(r["found_p"] for r in iii),
                     "calls_le_bound": sum(r["calls_le_bound"] for r in iii),
                     "median_calls_per_log2H": ratios[len(ratios) // 2], "max_calls_per_log2H": ratios[-1],
                     "iv_queries": len(iv), "iv_mismatches": sum(r["count"] != r["count_iv"] for r in iv)}
    for name, kind in (("x3", "MIX"), ("x4", "MIX2")):
        rs = [r for r in _load(name) if r.get("kind") == kind]
        if not rs:
            continue
        d_img = [r["Xmax"] - r["Xmax_U_img"] for r in rs if r["Xmax"] is not None and r["Xmax_U_img"] is not None]
        d_dim = [r["Xmax"] - r["Xmax_U_dim"] for r in rs if r["Xmax"] is not None and r["Xmax_U_dim"] is not None]
        o = {"cells": len(rs), "max_gap_vs_Uimg": max(d_img) if d_img else None,
             "max_gap_vs_Udim": max(d_dim) if d_dim else None,
             "cells_Xmax_None": sum(r["Xmax"] is None for r in rs)}
        if kind == "MIX":
            m1one = [r for r in rs if r["m1"] == 1]
            o["m1_1_cells"] = len(m1one)
            o["m1_1_good_beyond_Uimg_plus4"] = sum(bool(r["good_beyond_Uimg_plus4"]) for r in m1one)
            o["nontrivial_beyond_Uimg_plus3_total"] = sum(r["nontrivial_beyond_Uimg_plus3"] for r in rs)
            dF = [r["max_log2X_good"] - r["max_log2X_Fshort"] for r in m1one
                  if r["max_log2X_good"] is not None and r["max_log2X_Fshort"] is not None]
            o["m1_1_goodextent_minus_Fpred"] = [min(dF), max(dF)] if dF else None
        out[name] = o
    x5r = [r for r in _load("x5") if r.get("kind") == "HM"]
    if x5r:
        d = [r["Xmax"] - r["Xmax_U_dim"] for r in x5r if r["Xmax"] is not None and r["Xmax_U_dim"] is not None]
        out["x5"] = {"cells": len(x5r), "max_gap_vs_Udim": max(d) if d else None,
                     "line_fraction": sum(r["good_line_total"] for r in x5r) / max(1, sum(r["good_total"] for r in x5r)),
                     "cells_Xmax_None": sum(r["Xmax"] is None for r in x5r)}
    x6r = [r for r in _load("x6") if r.get("kind") == "SVP"]
    if x6r:
        d = [r["Xmax_SVP"] - r["Xmax_LLL"] for r in x6r if r["Xmax_SVP"] is not None and r["Xmax_LLL"] is not None]
        ratios = []
        for r in x6r:
            for s in r["scan"]:
                if r["Xmax_LLL"] is not None and s["log2X"] > r["Xmax_LLL"]:
                    ratios.append(s["log2_lambda1"] - s["log2_gh"])
        out["x6"] = {"cells": len(x6r), "gap_SVP_minus_LLL": [min(d), max(d)] if d else None,
                     "log2_lambda1_over_GH_beyond": [min(ratios), max(ratios)] if ratios else None}
    os.makedirs(ANA, exist_ok=True)
    json.dump(out, open(os.path.join(ANA, "topblock_summary.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    {"x1": x1, "x2": x2, "x3": x3, "x4": x4, "x5": x5, "x6": x6, "summary": summary}[sys.argv[1]]()
