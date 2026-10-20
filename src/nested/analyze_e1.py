"""E1 analysis (notes/nested_prereg.md): P1-P4 with bootstrap CIs over parent nodes.

Usage: .venv/bin/python -m src.nested.analyze_e1 experiments/raw/nested_e1.jsonl experiments/raw/nested_e1_summary.json
"""
from __future__ import annotations

import json
import random
import sys
from collections import Counter


def _rates(parents, kind, field):
    e = b = n = 0
    for p in parents:
        for c in p["children"]:
            if c["kind"] != kind:
                continue
            if field == "A":
                e += c["d1"]; b += c["base_d1"]; n += 1
            else:
                if c.get("d1_A2") is None:
                    continue
                e += c["d1_A2"]; b += c["base_d1_A2"]; n += 1
    return e, b, n


def diff_ci(parents, kind, field, B=2000, seed=0):
    rng = random.Random(seed)
    e, b, n = _rates(parents, kind, field)
    point = (e - b) / n if n else float("nan")
    ds = []
    for _ in range(B):
        samp = [parents[rng.randrange(len(parents))] for _ in range(len(parents))]
        e2, b2, n2 = _rates(samp, kind, field)
        if n2:
            ds.append((e2 - b2) / n2)
    ds.sort()
    return {"n": n, "rate_endo": e / n if n else None, "rate_unif": b / n if n else None, "diff": point,
            "ci95": [ds[int(0.025 * len(ds))], ds[int(0.975 * len(ds)) - 1]]}


def main():
    rows = [json.loads(l) for l in open(sys.argv[1])]
    parents = [p for r in rows for p in r["parents"]]
    out = {"instances": len(rows), "parents": len(parents), "S": rows[0].get("S")}
    out["P1_alg_A"] = diff_ci(parents, "a", "A")
    out["P2_rat_A"] = diff_ci(parents, "r", "A")
    out["P4_rat_A2"] = diff_ci(parents, "r", "A2")
    # pooled matched-uniform rates per type
    pt = {}
    for k in ("a", "r"):
        e, b, n = _rates(parents, k, "A")
        pt[k] = b / n if n else 0.0
    multi = [p for p in parents if len(p["children"]) >= 2]
    obs = sum(all(c["d1"] for c in p["children"]) for p in multi) / max(1, len(multi))
    exp = 0.0
    for p in multi:
        pr = 1.0
        for c in p["children"]:
            pr *= pt[c["kind"]]
        exp += pr
    exp /= max(1, len(multi))
    out["P3"] = {"nodes": len(multi), "obs_all_good": obs, "expected_indep": exp,
                 "ratio": obs / exp if exp else None}
    nch = Counter(len(p["children"]) for p in parents)
    out["children_hist"] = dict(sorted(nch.items()))
    out["frac_nodes_gt6_children"] = sum(v for k, v in nch.items() if k > 6) / len(parents)
    out["degenerate_children"] = sum(p["degenerate"] for p in parents)
    out["dead"] = {r["seed"]: r["dead"] for r in rows if "dead" in r}
    out["time_total_s"] = sum(r["time"] for r in rows)
    win = 0.08
    out["verdict"] = {
        "P1": win >= max(abs(x) for x in out["P1_alg_A"]["ci95"]),
        "P2": win >= max(abs(x) for x in out["P2_rat_A"]["ci95"]),
        "P3": out["P3"]["ratio"] is not None and 0.75 <= out["P3"]["ratio"] <= 1.33,
        "P4": win >= max(abs(x) for x in out["P4_rat_A2"]["ci95"]),
    }
    json.dump(out, open(sys.argv[2], "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
