"""E1 (pre-registered, notes/lv5unif_prereg.md): exact check of Lemma Kern-U.

Usage: .venv/bin/python -m src.lv5unif.e1_kern experiments/raw/lv5unif_e1.jsonl experiments/raw/lv5unif_e1_summary.json
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction as Fr
from math import gcd

from .kern import fibre_count_lattice, fibre_histogram_brute, triangular_bound, vvec, phi

SEED = 20261020
DS = [3, 4, 5, 6, 7, 9]
BOUND_CAP = 3 * 10 ** 6
VS = [3, 4, 6]
PAIRS_PER_D = 40
FIBRES_RANDOM = 8
BRUTE_LIMIT = 2 * 10 ** 6


def draw_pair(rnd: random.Random, V: int):
    while True:
        b, bp = rnd.randint(V // 2 + 1, V), rnd.randint(V // 2 + 1, V)
        a, ap = rnd.randint(-V, V), rnd.randint(-V, V)
        if gcd(a, b) == 1 and gcd(ap, bp) == 1 and a * bp - ap * b != 0:
            return (a, b), (ap, bp)


def run(out_path: str, summary_path: str) -> dict:
    rnd = random.Random(SEED)
    rows = []
    for d in DS:
        for j in range(PAIRS_PER_D):
            V = VS[j % len(VS)]
            Hs = [V * V, 2 * V * V, d * V * V]
            H = Hs[(j // len(VS)) % len(Hs)]
            redraws = 0
            while True:
                P, Pp = draw_pair(rnd, V)
                b, bp = P[1], Pp[1]
                bound = triangular_bound(d, b, bp, H)
                if bound <= BOUND_CAP:
                    break
                redraws += 1
                if redraws > 200:   # give up on this (d, V, H) cell; recorded
                    break
            if bound > BOUND_CAP:
                rows.append({"d": d, "V": V, "H": H, "skipped": True, "redraws": redraws})
                continue
            c0s = [[0] * d] + [[rnd.randint(-H, H) for _ in range(d)] for _ in range(FIBRES_RANDOM)]
            brute = (2 * H + 1) ** d <= BRUTE_LIMIT
            hist = fibre_histogram_brute(d, P, Pp, H) if brute else None
            counts = []
            cross_ok = True
            for c0 in c0s:
                cnt = fibre_count_lattice(d, P, Pp, H, c0)
                if hist is not None:
                    z = phi(c0, vvec(d, *P), vvec(d, *Pp))
                    if hist.get(z, 0) != cnt:
                        cross_ok = False
                counts.append(cnt)
            maxfib_all = max(hist.values()) if hist is not None else None
            maxcount = max(counts + ([maxfib_all] if maxfib_all is not None else []))
            # P3 statistic: (2H+1)^2 (bb')^(d-2) * maxfib / (2H+1)^d, as an exact fraction
            p3 = Fr((2 * H + 1) ** 2 * (b * bp) ** (d - 2) * maxcount, (2 * H + 1) ** d)
            rows.append({
                "d": d, "V": V, "H": H, "P": P, "Pp": Pp, "bound": bound,
                "central": counts[0], "counts": counts, "brute": brute,
                "max_fibre_all_z": maxfib_all, "cross_check_ok": cross_ok,
                "violation": maxcount > bound,
                "window": b > V * (1 - Fr(1, d)) and bp > V * (1 - Fr(1, d)),
                "H_ge_dV2": H >= d * V * V,
                "p3_num": p3.numerator, "p3_den": p3.denominator, "redraws": redraws, "skipped": False,
            })
    with open(out_path, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    # summary
    skipped = [r for r in rows if r.get("skipped")]
    rows = [r for r in rows if not r.get("skipped")]
    viol = sum(r["violation"] for r in rows)
    cross_bad = sum(not r["cross_check_ok"] for r in rows)
    per_d = {}
    p2_ok = True
    p3_max = Fr(0)
    for d in DS:
        rs = [r for r in rows if r["d"] == d]
        win = [Fr(r["central"], r["bound"]) for r in rs if r["window"] and r["H_ge_dV2"]]
        win.sort()
        med = win[len(win) // 2] if win else None
        if med is not None and med < Fr(1, 20):
            p2_ok = False
        p3d = max((Fr(r["p3_num"], r["p3_den"]) for r in rs if r["H_ge_dV2"]), default=Fr(0))
        p3_max = max(p3_max, p3d)
        per_d[d] = {
            "pairs": len(rs), "brute_pairs": sum(r["brute"] for r in rs),
            "violations": sum(r["violation"] for r in rs),
            "window_H_dV2_pairs": len(win),
            "median_central_over_bound_window": None if med is None else [med.numerator, med.denominator],
            "median_float_for_display": None if med is None else float(med),
            "max_p3_H_ge_dV2": [p3d.numerator, p3d.denominator], "max_p3_float_for_display": float(p3d),
        }
    summary = {
        "P1_violations": viol, "P1_holds": viol == 0,
        "cross_check_failures": cross_bad,
        "P2_holds": p2_ok, "P2_note": "median over window pairs with H = dV^2 of central/bound >= 0.05 for each d with such pairs",
        "P3_max": [p3_max.numerator, p3_max.denominator], "P3_max_float_for_display": float(p3_max),
        "P3_holds": p3_max <= Fr(66, 10),
        "per_d": per_d, "rows": len(rows), "skipped_cells": len(skipped),
        "total_redraws": sum(r["redraws"] for r in rows) + sum(r["redraws"] for r in skipped),
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=1)
    return summary


if __name__ == "__main__":
    s = run(sys.argv[1], sys.argv[2])
    print(json.dumps({k: v for k, v in s.items() if k != "per_d"}, indent=1))
