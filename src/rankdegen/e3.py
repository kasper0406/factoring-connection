"""Experiment E3 (pre-registered, notes/rankdegen_prereg.md): box points with small |F| in excess of the R1 volume
bound, classified by exceptional type (linear candidates L, products of linear candidates Rd, irreducible G).

Exact: integer norms via resultants; the main term is a certified arb ball (upper end used).
Run: PYTHONPATH=. .venv/bin/python -m src.rankdegen.e3   (writes experiments/raw/rankdegen_E3.jsonl)
"""
from __future__ import annotations

import bisect
import json
import math
import time

from src.rank3.forms import count_bound_terms, delta_F, poly_comb, resultant

OUT = "experiments/raw/rankdegen_E3.jsonl"
E = 40


def is_square(n: int) -> bool:
    if n < 0:
        return False
    s = math.isqrt(n)
    return s * s == n


def run_one(r_):
    f, g, P, Q = r_["f"], r_["g"], r_["P"], r_["Q"]
    info = delta_F(f, P, g, Q, 3)
    vals = {"L": [], "Rd": [], "G": []}
    for c in range(-E, E + 1):
        for b in range(-E, E + 1):
            for a in range(-E, E + 1):
                if a == 0 and b == 0 and c == 0:
                    continue
                A = poly_comb(P, (a, b, c))
                w = abs(resultant(f, A)) * abs(resultant(g, poly_comb(Q, (a, b, c))))
                if w == 0:
                    continue
                if c == 0:
                    vals["L"].append(w)
                elif is_square(b * b - 4 * a * c):
                    vals["Rd"].append(w)
                else:
                    vals["G"].append(w)
    for k in vals:
        vals[k].sort()
    allv = sorted(vals["L"] + vals["Rd"] + vals["G"])
    out = []
    p8 = False
    p9 = True
    last_excess = None
    for k in range(0, allv[-1].bit_length() + 1):
        X = 1 << k
        cnt = {t: bisect.bisect_right(vals[t], X) for t in vals}
        T = sum(cnt.values())
        main, _ = count_bound_terms(3, 6, E, X, info)
        mu = float(main.upper())
        if T > mu:
            p8 = True
            last_excess = k
        ok9 = cnt["G"] <= 2 * mu + 20
        p9 &= ok9
        out.append({"log2X": k, "T": T, **cnt, "main_upper": mu, "P9_ok": ok9})
    return {"kind": r_["kind"], "seed": r_["seed"], "N": r_["N"], "f": f, "g": g, "E": E,
            "Q2_bits": abs(int(info["Q2"])).bit_length(), "P8": p8, "P9": p9,
            "last_excess_log2X": last_excess, "grid": out}


def main():
    rows = [json.loads(l) for l in open("experiments/raw/rank3_E1.jsonl")]
    sel = [r_ for r_ in rows if r_["r"] == 3 and (r_["d"], r_["d2"]) == (3, 3) and r_["kind"] in ("A", "B")]
    n8 = n9 = 0
    with open(OUT, "w") as fhd:
        for r_ in sel:
            t0 = time.time()
            res = run_one(r_)
            res["secs"] = round(time.time() - t0, 1)
            fhd.write(json.dumps(res) + "\n")
            fhd.flush()
            n8 += res["P8"]
            n9 += res["P9"]
            worst = max(res["grid"], key=lambda z: z["G"] / max(z["main_upper"], 1e-300))
            print(res["kind"], res["seed"], "N", res["N"], "P8", res["P8"], "P9", res["P9"], "last excess k",
                  res["last_excess_log2X"], "worst G/main at k", worst["log2X"], worst["G"], round(worst["main_upper"], 3),
                  res["secs"], "s", flush=True)
    print(f"P8 held in {n8}/{len(sel)}; P9 held in {n9}/{len(sel)}")


if __name__ == "__main__":
    main()
