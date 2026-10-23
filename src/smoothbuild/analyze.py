"""Analysis of E1-E3 against the pre-registered predictions (notes/smoothbuild_prereg.md).

    .venv/bin/python -m src.smoothbuild.analyze
Writes experiments/analysis/smoothbuild_summary.json and prints a table.
Floating point is used only for the reported ratios.
"""
from __future__ import annotations

import json
import statistics
from collections import defaultdict

RAW = "experiments/raw/smoothbuild_E{}.jsonl"
OUT = "experiments/analysis/smoothbuild_summary.json"


def load(k):
    with open(RAW.format(k)) as fh:
        return [json.loads(l) for l in fh if l.strip()]


def within(x, lo, hi):
    return lo <= x <= hi


def e1():
    rows = load(1)
    by = defaultdict(list)
    for r in rows:
        by[r["family"]].append(r)
    res = {"families": {}}
    mean_S = {}
    for fam, rs in by.items():
        mS = statistics.mean(r["S"] for r in rs)
        mean_S[fam] = mS
        ratio = statistics.mean(r["ratio_cheb"] for r in rs)
        d = {"mean_S": mS, "mean_ratio_cheb": ratio, "rank": rs[0].get("rank"),
             "ratios": [round(r["ratio_cheb"], 4) for r in rs]}
        if "ratio_arr" in rs[0]:
            d["mean_ratio_arr"] = statistics.mean(r["ratio_arr"] for r in rs)
        if fam.startswith("R:"):
            d["P1b"] = within(ratio, 0.90, 1.10)
        else:
            d["P1a"] = within(ratio, 0.80, 1.25)
        res["families"][fam] = d
    s4 = mean_S["S4:x4+x+1"]
    s6 = mean_S["S6:x6+x+1"]
    c = {
        "Phi5/S4": (mean_S["C4:Phi5"] / s4, 2.2, 4.0),
        "Phi8/S4": (mean_S["V4:Phi8"] / s4, 2.2, 4.0),
        "D4/S4": (mean_S["D4:x4-2"] / s4, 1.5, 2.6),
        "Phi7/S6": (mean_S["C6:Phi7"] / s6, 3.5, 6.5),
        "Phi9/S6": (mean_S["C6:Phi9"] / s6, 3.5, 6.5),
        "A4/S4": (mean_S["A4:x4+8x+12"] / s4, 0.75, 1.33),
    }
    res["P1c"] = {k: {"ratio": v[0], "band": [v[1], v[2]], "pass": within(*v)} for k, v in c.items()}
    return res


def e2():
    rows = load(2)
    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by[r["bits"]][r["method"]].append(r)
    res = {}
    for bits, ms in sorted(by.items()):
        h0 = sum(r["hits"] for r in ms["M0"]); n0 = sum(r["null_hits"] for r in ms["M0"])
        h1 = sum(r["hits"] for r in ms["M1"]); n1 = sum(r["null_hits"] for r in ms["M1"])
        t1 = sum(r["trials"] for r in ms["M1"])
        h2 = sum(r["hits"] for r in ms["M2"]); t2 = sum(r["trials"] for r in ms["M2"])
        th2 = sum(r["target_hits"] for r in ms["M2"]); tt2 = sum(r["targets"] for r in ms["M2"])
        h4 = sum(r["hits"] for r in ms["M4"]); e4 = sum(r["expected"] for r in ms["M4"])
        t0 = sum(r["trials"] for r in ms["M0"])
        p1 = h1 / t1
        p2 = h2 / t2
        pred_target = 1 - (1 - p2) ** 16
        d = {
            "M0_rate": h0 / t0, "M0_O/E": h0 / n0, "P2a": within(h0 / n0, 0.85, 1.18),
            "M1_rate": p1, "M1_O/E": h1 / n1, "P2b": within(h1 / n1, 0.75, 1.25),
            "M2_split_rate": p2, "M2/M1": p2 / p1, "P2c": within(p2 / p1, 0.90, 1.10),
            "M2_target_rate": th2 / tt2, "M2_target/pred": (th2 / tt2) / pred_target,
            "P2d": within((th2 / tt2) / pred_target, 0.90, 1.10),
            "M4_hits": h4, "M4_expected": e4, "M4_O/E": h4 / e4, "P2e": within(h4 / e4, 0.85, 1.18),
            "relations_per_test": {"M0": h0 / t0, "M1": p1 / 2, "M2": p2 / 2},
            "counts": {"M0": [h0, n0, t0], "M1": [h1, n1, t1], "M2": [h2, t2, th2, tt2]},
        }
        res[bits] = d
    return res


def e3():
    rows = load(3)
    by = defaultdict(list)
    for r in rows:
        by[r["bits"]].append(r)
    res = {}
    for bits, rs in sorted(by.items()):
        ratios = [r["T_over_J0"] for r in rs]
        tb = [r["Tbest_over_T"] for r in rs]
        # post hoc: fixed-in-advance y* strategy
        tys = []
        for r in rs:
            v = r["Ty"].get(str(r["ystar"]))
            tys.append(v / r["T"] if v else None)
        fin = [x for x in tys if x is not None]
        res[bits] = {
            "J0hat": rs[0]["J0hat"], "ystar": rs[0]["ystar"],
            "T_over_J0": [round(x, 4) for x in ratios],
            "median_T_over_J0": statistics.median(ratios), "P3a": within(statistics.median(ratios), 0.35, 0.80),
            "min_T_over_J0": min(ratios), "P3b": min(ratios) >= 0.25,
            "median_Tbest_over_T": statistics.median(tb), "P3c": within(statistics.median(tb), 1.0, 2.0),
            "posthoc_Tystar_over_T": [round(x, 4) if x else None for x in tys],
            "posthoc_median_Tystar_over_T": statistics.median(fin) if fin else None,
            "posthoc_Tystar_reached": f"{len(fin)}/{len(tys)}",
            "dep_max_prime_over_ystar_median": statistics.median(r["dep_max_prime"] / r["ystar"] for r in rs),
            "dep_size_median": statistics.median(r["dep_size"] for r in rs),
        }
    return res


def main():
    out = {"E1": e1(), "E2": e2(), "E3": e3()}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
