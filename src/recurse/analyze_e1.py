"""Analysis of recurse E1 against the frozen predictions P1-P5 (notes/recurse_prereg.md)."""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict


def main(path: str, out: str | None = None):
    recs = [json.loads(l) for l in open(path)]
    S = {"n": len(recs)}
    S["P1_correct"] = {o: [sum(r["runs"][o].get("ok", False) for r in recs), len(recs)] for o in ("pb", "rec")}
    S["P2_outer_identical"] = [sum(r["outer_identical"] for r in recs), len(recs)]
    big = [r for r in recs if r["bits"] >= 50 and r["spec"]["family"] in ("semi", "tri")]
    ex = [r for r in big if any(s["depth"] > 0 for s in r["runs"]["rec"].get("splits", []))]
    S["P3_exercised_ge50"] = [len(ex), len(big), len(ex) >= 0.8 * len(big)]
    S["P3_max_depth"] = max(r["runs"]["rec"].get("max_depth", 0) for r in recs)
    viol = sum(s["row_viol"] for r in recs for o in ("pb", "rec") for s in r["runs"][o].get("splits", []))
    S["P4_row_violations"] = viol
    dd = defaultdict(lambda: [0, 0])
    inner = defaultdict(lambda: [0, 0])
    for r in recs:
        for s in r["runs"]["rec"].get("splits", []):
            for f in s["fruitful"]:
                dd[s["omega"]][0] += 1; dd[s["omega"]][1] += bool(f)
                if s["depth"] > 0:
                    inner[s["omega"]][0] += 1; inner[s["omega"]][1] += bool(f)
    p5 = {}
    for name, tab in (("all", dd), ("inner_only", inner)):
        for om, (n, k) in sorted(tab.items()):
            p = 1 - 2 ** (1 - om)
            z = (k - n * p) / math.sqrt(n * p * (1 - p)) if n else None
            p5[f"{name}|omega={om}"] = {"tested": n, "fruitful": k, "rate": k / n if n else None, "pred": p,
                                        "z": z, "holds": (abs(z) <= 3) if z is not None else None}
    S["P5_DD"] = p5
    kinds = defaultdict(int)
    for r in recs:
        for s in r["runs"]["rec"].get("splits", []):
            kinds[f"depth{min(s['depth'], 2)}|{s['kind'].split('_')[0]}"] += 1
    S["split_kinds_rec"] = dict(kinds)
    S["time_ratio_rec_over_pb"] = {f"{r['spec']['family']}{r['bits']}s{r['spec']['seed']}":
                                   round(r["runs"]["rec"].get("seconds", 0) / max(r["runs"]["pb"].get("seconds", 1e-9), 1e-3), 1)
                                   for r in recs}
    txt = json.dumps(S, indent=1)
    print(txt)
    if out:
        open(out, "w").write(txt)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
