"""uprime E2 POST HOC (not pre-registered): within-field null for the conditional child mean.

The registered null (experiments.py) draws labels from a pool of FRESH fields, so R_CM mixes two effects: conditioning
on the parent (the quantity of interest) and field quality (fields that accept more level-1 points may also have more
good child points; multidescent E3 saw such height/Jensen effects).  This script reads the registered E2 output,
reconstructs each field exactly (same seeds and RNG streams as experiments.py), and draws for every child a random
simple label of the SAME field with matched size.  R_CM_within = mean conditional / mean within-field null.

Run: .venv/bin/python -m src.uprime.posthoc experiments/raw/uprime_e2.jsonl experiments/raw/uprime_e2_posthoc.json
"""
from __future__ import annotations

import json
import random
import sys

from src.common.instances import balanced_semiprime
from src.uprime import toy
from src.uprime.analyze import _boot
from src.uprime.experiments import CONFIG


def main(inp: str, out: str) -> None:
    rows = [json.loads(l) for l in open(inp)]
    cfg = rows[0]["cfg"]
    assert cfg == {k: v for k, v in CONFIG.items()}
    kids = [r for r in rows[1:] if r.get("kind") == "child"]
    ctx = toy.Ctx(cfg["B"], cfg["B2"])
    fields = {}
    for s in sorted({k["seed"] for k in kids}):
        N = balanced_semiprime(cfg["n"], s).N
        prng = random.Random(f"uprime-E2-pool|{s}")
        for _ in range(cfg["pool"]):
            toy.sample_field(N, prng, d=cfg["d"], M_bits=cfg["M_bits"])
        for fi in range(cfg["fields"]):
            rng = random.Random(f"uprime-E2|{s}|{fi}")
            fields[(s, fi)] = toy.sample_field(N, rng, d=cfg["d"], M_bits=cfg["M_bits"])
    hr = random.Random("uprime-E2-posthoc-within")
    res = []
    for k in kids:
        A = fields[(k["seed"], k["field"])]
        q = k["q"]
        lab = toy.random_label(A, q / cfg["match"], q * cfg["match"], hr, ctx)
        if lab is None:
            continue
        w = toy.child_count(A, lab[0], lab[1], cfg["S2"], ctx)
        res.append((k, w["good"]))
    mc = sum(k["cond"]["good"] for k, _ in res) / len(res)
    mw = sum(g for _, g in res) / len(res)
    parents = {}
    for k, g in res:
        parents.setdefault((k["seed"], k["field"], k["ell"]), []).append((k["cond"]["good"], g))
    pl = list(parents.values())

    def r(sample):
        return sum(a for grp in sample for a, _ in grp) / max(1, sum(b for grp in sample for _, b in grp))

    byf = {}
    for k, g in res:
        byf.setdefault((k["seed"], k["field"]), []).append((k["cond"]["good"], g))
    summary = {"children": len(res), "mean_cond": mc, "mean_within_null": mw, "R_CM_within": mc / mw,
               "CI_parents": _boot(pl, r), "CI_fields": _boot(list(byf.values()), r),
               "per_field": {f"{a}-{b}": r([v]) for (a, b), v in sorted(byf.items())}}
    json.dump(summary, open(out, "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
