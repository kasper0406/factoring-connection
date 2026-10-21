"""Analysis of E1–E3 (registered statistics plus the POST-HOC ones named in the results note).
Usage: .venv/bin/python -m src.siconstruct.analyze"""
from __future__ import annotations

import json
import os

from .core import ht_local_density, primes_upto
from .run_exp import RAW


def _load(name):
    return [json.loads(l) for l in open(os.path.join(RAW, f"siconstruct_{name}.jsonl"))]


def main():
    E1 = _load("E1")
    H = sum(r["hits"] for r in E1); P = sum(r["pred_hits"] for r in E1)
    print("E1 P1a pooled", round(H / P, 3), "P1b cells", round(min(r["ratio"] for r in E1), 3), round(max(r["ratio"] for r in E1), 3))
    byK = {}
    for r in E1:
        a = byK.setdefault(r["K"], [0, 0]); a[0] += r["hits"]; a[1] += r["pred_hits"]
    rk = {K: v[0] / v[1] for K, v in byK.items()}
    print("E1 P1c per-K", {K: round(v, 3) for K, v in rk.items()}, "max/min", round(max(rk.values()) / min(rk.values()), 3))
    E2 = _load("E2")
    m = [r for r in E2 if r["method"] == "mitm"]
    print("E2 P2a pooled", round(sum(r["hits"] for r in m) / sum(r["pred_hits"] for r in m), 3))
    print("E2 P2b", [round(r["cost_per_output"] / (5 * r["M"] / 64), 2) for r in m])
    for r in E2:
        if r["method"] == "brute" and "measured_density" in r:
            Pp = primes_upto(r["y"])
            avg = sum(ht_local_density(float(k * r["N"]), r["y"], Pp) for k in range(1, 9)) / 8
            print("E2 P2c", r["bits"], r["seed"], "registered", round(r["measured_density"] / r["ht_density"], 3),
                  "POST-HOC position-matched", round(r["measured_density"] / avg, 3))
    E3 = _load("E3")
    cells = {}
    for r in E3:
        Pp = primes_upto(r["y"])
        d = ht_local_density(float(r["N"]) / (1 << r["S0exp"]), r["y"], Pp)
        print("E3 P3c", r["bits"], r["seed"], r["S0exp"], "tests/out * delta(N/S0) =", round(r["per_output_tests"] * d, 2))
        c = cells.setdefault((r["bits"], r["S0exp"]), [0, 0.0]); c[0] += r["outputs"]; c[1] += r["pred_outputs"]
    for k, v in sorted(cells.items()):
        print("E3 P3a pooled cell", k, round(v[0] / v[1], 3))


if __name__ == "__main__":
    main()
