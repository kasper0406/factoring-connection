"""HARNESS: summarise Track A raw JSONL into experiments/analysis/trackA_summary.md.

Usage: .venv/bin/python -m src.trackA.analyze
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict

RAW = "experiments/raw"
OUT = "experiments/analysis/trackA_summary.md"


def load(name):
    path = os.path.join(RAW, name)
    if not os.path.exists(path):
        return []
    return [json.loads(l) for l in open(path) if l.strip()]


def descent_section(lines):
    recs = load("trackA_descent.jsonl")
    if not recs:
        return
    lines.append("## Descent / selectors (trackA_descent.jsonl)\n")
    # per level relative to log2 p
    both = defaultdict(lambda: [0, 0])
    rand_alive = defaultdict(lambda: [0, 0])
    exact_ok = defaultdict(lambda: [0, 0])
    budget_ok = defaultdict(lambda: [0, 0])
    evals = defaultdict(list)
    for d in recs:
        lp = math.log2(d["p"])
        alive = True
        for row in d["levels"]:
            rel = row["level"] - round(lp)
            if row["parent_nonzero"]:
                both[rel][0] += row["child0_nonzero"] and row["child1_nonzero"]
                both[rel][1] += 1
            tgt = {"random": rand_alive, "exact_scalar": exact_ok}.get(d["selector"],
                                                                       budget_ok)
            if row["parent_nonzero"]:
                tgt[rel][0] += row["chosen_nonzero"]
                tgt[rel][1] += 1
            if d["selector"] == "exact_scalar":
                evals[rel].append(row["evals"])
    lines.append("Rows: level k minus round(log2 p). `both` = fraction of live parents whose two "
                 "children are both nonzero (selection is free). Success = P(chosen child nonzero | "
                 "parent nonzero).\n")
    lines.append("| k - log2 p | both children nonzero | random success | exact success | "
                 "budget_1024 success | exact: mean evals |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for rel in sorted(both):
        def f(x):
            return f"{x[0]}/{x[1]}" if x[1] else "-"
        ev = evals.get(rel, [])
        lines.append(f"| {rel} | {f(both[rel])} | {f(rand_alive[rel])} | {f(exact_ok[rel])} | "
                     f"{f(budget_ok[rel])} | {sum(ev)/len(ev):.2f} |" if ev else
                     f"| {rel} | {f(both[rel])} | {f(rand_alive[rel])} | {f(exact_ok[rel])} | "
                     f"{f(budget_ok[rel])} | - |")
    # cost per level for exact
    lines.append("\nExact selector cost per level (q_mul_dimsum increment = sum over quotient "
                 "products of the algebra rank; one level evaluates rank-2M algebras):\n")
    lines.append("| bits | seed | level | M | rank | cost this level | ratio to previous |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for d in recs:
        if d["selector"] != "exact_scalar" or d["seed"] != 0:
            continue
        prev_c = 0
        prev_inc = None
        for row in d["levels"]:
            inc = row["cum_q_mul_dimsum"] - prev_c
            prev_c = row["cum_q_mul_dimsum"]
            if row["level"] % 3 == 0 or row["level"] == len(d["levels"]):
                ratio = f"{inc/prev_inc:.2f}" if prev_inc else "-"
                lines.append(f"| {d['bits']} | {d['seed']} | {row['level']} | {row['M']} | "
                             f"{row['rank']} | {inc} | {ratio} |")
            prev_inc = inc
    # discoveries
    lines.append("\nGCD discoveries from section values (exact selector):\n")
    lines.append("| bits | seed | p | q-p | first discovery level | 2^level | kinds |")
    lines.append("|---:|---:|---:|---:|---:|---:|---|")
    for d in recs:
        if d["selector"] != "exact_scalar":
            continue
        ds = d["discoveries"]
        kinds = defaultdict(int)
        for x in ds:
            kinds[x["kind"]] += 1
        first = min((x["level"] for x in ds), default=None)
        lines.append(f"| {d['bits']} | {d['seed']} | {d['p']} | {d['q']-d['p']} | {first} | "
                     f"{2**(first+1) if first is not None else '-'} | {dict(kinds)} |")
    lines.append("")


def circuits_section(lines):
    recs = load("trackA_circuits.jsonl")
    if not recs:
        return
    lines.append("## Circuit section growth (trackA_circuits.jsonl)\n")
    lines.append("Unique reachable nodes of the section circuit after k levels "
                 "(hash-consed, constants folded).\n")
    lines.append("| bits | seed | path | variant | nodes k=0..9 | mean growth ratio k=5..9 |")
    lines.append("|---:|---:|---|---|---|---:|")
    for d in recs:
        n = [r["nodes"] for r in d["rows"]]
        rat = [n[i] / n[i - 1] for i in range(5, len(n))]
        lines.append(f"| {d['bits']} | {d['seed']} | {d['path']} | "
                     f"{'karatsuba' if d['karatsuba'] else 'naive'} | {n} | "
                     f"{sum(rat)/len(rat):.2f} |")
    lines.append("")


def zero_section(lines):
    recs = load("trackA_zero_tests.jsonl")
    if not recs:
        return
    lines.append("## False-zero catalogue (trackA_zero_tests.jsonl)\n")
    agg = defaultdict(lambda: defaultdict(int))
    for d in recs:
        fam = d["family"]
        agg[fam]["instances"] += 1
        for r in d["levels"]:
            agg[fam]["nonzero_sections"] += r["nonzero_classes"]
            for k in ("false_zero_scalar", "false_zero_aks", "false_zero_monic2", "false_zero_all"):
                agg[fam][k] += r[k]
    lines.append("| family | instances | nonzero sections | scalar false zeros | "
                 "AKS(s<=4,a<=3) false zeros | random monic quadratic (4 trials) false zeros |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for fam, a in agg.items():
        lines.append(f"| {fam} | {a['instances']} | {a['nonzero_sections']} | "
                     f"{a['false_zero_scalar']} | {a['false_zero_aks']} | {a['false_zero_monic2']} |")
    lines.append("")


def theory_section(lines):
    recs = load("trackA_theory.jsonl")
    if not recs:
        return
    lines.append("## Theory checks (trackA_theory.jsonl)\n")
    a1 = [d for d in recs if d["check"] == "A1_exhaustive"]
    a2 = [d for d in recs if d["check"] == "A2_order"]
    lines.append(f"- A1 exhaustive (F_2, F_3, L=20): {sum(d['ok'] for d in a1)}/{len(a1)} cases "
                 f"satisfy min linear complexity >= predicted bound; tight in "
                 f"{sum(d['min_linear_complexity'] == d['predicted_lower'] for d in a1)} cases.")
    lines.append(f"- A2 Berlekamp-Massey order = M: {sum(d['ok'] for d in a2)}/{len(a2)}.")
    lines.append("\n| bits | seed | p | q-p | first M with a pure class | min(p-1,q-p) |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for d in recs:
        if d["check"] == "C_props":
            lines.append(f"| {d['bits']} | {d['seed']} | {d['p']} | {d['q_minus_p']} | "
                         f"{d['first_pure_M']} | {d['min_pm1_qmp']} |")
    lines.append("")


def masks_section(lines):
    recs = load("trackA_masks.jsonl")
    if not recs:
        return
    lines.append("## Mask costs (trackA_masks.jsonl)\n")
    lines.append("| bits | mask | rank | q_mul_dimsum | wall s |")
    lines.append("|---:|---|---:|---:|---:|")
    for d in recs:
        if d["log2_rank"] in (8, 14, 20):
            lines.append(f"| {d['bits']} | {d['kind']} | {d['rank']} | "
                         f"{d['counts'].get('q_mul_dimsum')} | {d['wall_s']:.3f} |")
    agree = [d.get("agree") for d in recs if d["kind"] == "prefix"]
    lines.append(f"\nPrefix mask nonzero <=> rank > p (L3): {sum(agree)}/{len(agree)} agree.\n")


def main():
    lines = ["# Track A experiment summary (generated by `python -m src.trackA.analyze`)\n"]
    theory_section(lines)
    circuits_section(lines)
    descent_section(lines)
    zero_section(lines)
    masks_section(lines)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
