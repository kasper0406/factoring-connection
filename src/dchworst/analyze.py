"""Evaluate the frozen predictions of notes/dchworst_prereg.md on experiments/raw/dchworst_*.
HARNESS.  Run:  PYTHONPATH=. .venv/bin/python -m src.dchworst.analyze  (writes
experiments/analysis/dchworst_summary.json and prints the tables used in the notes)."""
from __future__ import annotations

import json
import math
import os
from math import isqrt

import numpy as np

RAW = "experiments/raw"


def _load(name):
    path = os.path.join(RAW, f"dchworst_{name}.jsonl")
    if not os.path.exists(path):
        return []
    return [json.loads(l) for l in open(path)]


def tail_cells(rec, Hmin=0):
    """(H, obs, exp) cells with expected >= 20 (P3)."""
    out = []
    for H, (o, e) in enumerate(zip(rec["tail_obs"], rec["tail_exp"])):
        if H >= Hmin and e >= 20 and H <= rec.get("Hcap", 64):
            out.append((H, o, e))
    return out


def p3(cells):
    if not cells:
        return {"cells": 0}
    inside2 = sum(1 for _, o, e in cells if 0.5 * e <= o <= 2 * e)
    inside14 = sum(1 for _, o, e in cells if 0.7 * e <= o <= 1.4 * e)
    worst = max(cells, key=lambda c: abs(math.log(max(c[1], 0.5) / c[2])))
    return {"cells": len(cells), "in_[0.5,2]": inside2, "in_[0.7,1.4]": inside14,
            "frac_in_[0.7,1.4]": round(inside14 / len(cells), 4),
            "worst_cell": [worst[0], worst[1], round(worst[2], 1)],
            "pass": inside2 == len(cells) and inside14 >= 0.95 * len(cells)}


def e1():
    rows = _load("e1")
    out = {"P1": [], "P2": [], "P3": [], "P4": []}
    groups = {}
    for r in rows:
        groups.setdefault((r["n"], r["Y"], r["signed"]), []).append(r)
    for (n, Y, s), rs in sorted(groups.items()):
        rs.sort(key=lambda r: r["J"])
        if s:
            r1 = rs[0]
            out["P1"].append({"n": n, "Y": Y, "l1": r1["l_max"], "hstar1": r1["hstar"],
                              "bound": isqrt(r1["l_max"]), "pass": r1["hstar"] <= isqrt(r1["l_max"]),
                              "pairs_above_bound": r1["tail_obs"][isqrt(r1["l_max"])],
                              "pairs_with_t1_expected": round(r1["npairs"] / (r1["l_max"] - 1))})
        for r in rs:
            if r["hpred"] <= r["Hcap"]:
                ok = r["hpred"] - 4 <= r["hstar"] <= r["hpred"] + 3
                out["P2"].append({"n": n, "Y": Y, "signed": s, "J": r["J"], "hstar": r["hstar"],
                                  "hpred": r["hpred"], "pass": ok, "worst_pair": r["worst_pair"]})
        cells = [c for r in rs for c in tail_cells(r)]
        out["P3"].append({"n": n, "Y": Y, "signed": s, **p3(cells)})
        if Y in (2 ** 11, 2 ** 13):
            pts = [(r["J"], r["hstar"], r["hpred"]) for r in rs if r["J"] >= 16]
            x = np.log([p[0] for p in pts])
            c_obs = -np.polyfit(x, np.log([p[1] for p in pts]), 1)[0]
            c_pred = -np.polyfit(x, np.log([p[2] for p in pts]), 1)[0]
            out["P4"].append({"n": n, "Y": Y, "signed": s, "c_obs": round(c_obs, 4),
                              "c_pred": round(c_pred, 4), "pass": abs(c_obs - c_pred) <= 0.08,
                              "series": pts})
    return out


def e2():
    rows = _load("e2")
    out = {"P2": [], "P3": []}
    for r in rows:
        r = dict(r)
        r["Hcap"] = 64
        ok = r["hpred"] - 4 <= r["hstar"] <= r["hpred"] + 3 if r["hpred"] <= 64 else None
        out["P2"].append({"n": r["n"], "block": r["block"], "signed": r["signed"], "J": r["J"],
                          "hstar": r["hstar"], "hpred": r["hpred"], "pass": ok})
        out["P3"].append({"n": r["n"], "block": r["block"], "signed": r["signed"], **p3(tail_cells(r))})
    return out


def e3():
    rows = _load("e3")
    out = []
    base = {r["signed"]: r for r in rows if r["family"] == "random"}
    for r in rows:
        if r["family"] == "random":
            continue
        Href = 6 if r["signed"] else 8
        b = base[r["signed"]]
        fb = b["tail_obs"][Href] / b["npairs"]
        ff = r["tail_obs"][Href] / r["npairs"]
        ratio = ff / fb if fb else float("nan")
        lo, hi = (0.8, 1.25) if r["family"] == "crt_worst2" else (0.5, 2.0)
        sd = math.sqrt(fb * (1 - fb) / r["npairs"])
        z = (ff - fb) / sd if sd else 0.0
        out.append({"family": r["family"], "signed": r["signed"], "npairs": r["npairs"],
                    "Href": Href, "unc_frac": round(ff, 4), "base_frac": round(fb, 4),
                    "poisson_frac": round(r["tail_exp_frac"][Href], 4),
                    "ratio": round(ratio, 3), "z": round(z, 2), "pass": lo <= ratio <= hi,
                    "structured_adversary": ratio > 2 and z > 3})
    return out


def e4():
    rows = _load("e4")
    out = {"cells": [], "Hstar": []}
    by_n = {}
    for r in rows:
        by_n.setdefault(r["n"], []).append(r)
        if r["poisson_expected"] >= 20:
            ratio = r["uncovered"] / r["poisson_expected"]
            out["cells"].append({"n": r["n"], "H": r["H"], "obs": r["uncovered"],
                                 "exp": round(r["poisson_expected"], 1), "ratio": round(ratio, 3),
                                 "pass": 0.5 <= ratio <= 2})
    for n, rs in sorted(by_n.items()):
        rs.sort(key=lambda r: r["H"])
        Hs = next((r["H"] for r in rs if all(x["uncovered"] == 0 for x in rs if x["H"] >= r["H"])), None)
        Hp = next((r["H"] for r in rs if r["poisson_expected"] < 0.5), None)
        out["Hstar"].append({"n": n, "Hstar_obs": Hs, "Hstar_pred": Hp,
                             "pass": Hs is not None and Hp is not None and abs(Hs - Hp) <= 1,
                             "uncovered_clean": {r["H"]: r["uncovered_clean"] for r in rs}})
    return out


def e5():
    rows = _load("e5")
    out = []
    for r in rows:
        ra = r["uncovered_all4"] / r["poisson_all4"]
        rb = [o / e for o, e in zip(r["uncovered_per_block"], r["poisson_per_block"])]
        out.append({"H": r["H"], "Y": r["Y"], "all4_obs": r["uncovered_all4"],
                    "all4_exp": round(r["poisson_all4"], 1), "all4_ratio": round(ra, 3),
                    "block_ratios": [round(x, 3) for x in rb],
                    "pass_all4": 0.5 <= ra <= 2, "pass_blocks": all(0.8 <= x <= 1.25 for x in rb),
                    "examples": r["examples"][:3]})
    return out


def e6():
    return {"e6a": _load("e6a"), "e6b": _load("e6b")}


def main():
    s = {"e1": e1(), "e2": e2(), "e3": e3(), "e4": e4(), "e5": e5(), "e6": e6()}
    os.makedirs("experiments/analysis", exist_ok=True)
    with open("experiments/analysis/dchworst_summary.json", "w") as f:
        json.dump(s, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(json.dumps({k: v for k, v in s.items() if k != "e6"}, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))[:20000])


if __name__ == "__main__":
    main()
