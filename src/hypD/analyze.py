"""HARNESS: analysis of hypD E1-E3 raw JSONL against the pre-registered predictions.

Usage: .venv/bin/python -m src.hypD.analyze {e1,e2,e3} RAW.jsonl [RAW2.jsonl ...]
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from math import exp, sqrt


def load(paths):
    out = []
    for p in paths:
        with open(p) as fh:
            out += [json.loads(l) for l in fh if l.strip()]
    return out


def e1(recs):
    print("E1: Theorem S lower bound (P1) and coprime-pair model (P2)")
    fails, zs = 0, []
    pool = defaultdict(lambda: [0, 0.0])
    for r in recs:
        p_hat = r["hits"] / r["trials"]
        if p_hat < r["bound"] - 3 * sqrt(r["bound"] / r["trials"]):
            fails += 1
        mu = r["model"] * r["trials"]
        if mu > 0:
            zs.append((r["hits"] - mu) / sqrt(mu))
        key = (r["cls"], r["k"])
        pool[key][0] += r["hits"]
        pool[key][1] += mu
    print(f"  P1 records violating the bound: {fails} / {len(recs)}  (predicted 0)")
    for key, (h, mu) in sorted(pool.items()):
        print(f"  P2 {key}: hits {h}, model {mu:.1f}, ratio {h / mu:.3f}  (predicted in [0.8, 1.25])")
    ok = sum(1 for z in zs if abs(z) <= 3)
    print(f"  P2 per-record |z|<=3: {ok}/{len(zs)} = {ok / len(zs):.3f} (predicted >= 0.95)")
    mb = min(r["hits"] / r["trials"] / r["bound"] for r in recs if r["hits"] > 0)
    print(f"  min (empirical / Theorem-S bound) over records with hits: {mb:.3g}")


def e2(recs):
    print("E2: special-q counts N_l(A)")
    for r in sorted(recs, key=lambda x: (x["cls"], x["bits"], x["seed"])):
        Ns, nr, gA, cells = [], [], [], []
        pred_gamma = []
        for f in r["fields"]:
            pA = sum(f["N_l"]) / max(1, sum(f["nrat"]))
            gA.append(sum(1 for x in f["N_l"] if x > 0) / len(f["N_l"]))
            for n, k in zip(f["N_l"], f["nrat"]):
                Ns.append(n)
                nr.append(k)
                pred_gamma.append(1 - exp(-k * pA))
            cells.append([1 if x > 0 else 0 for x in f["N_l"]])
        M = len(Ns)
        mu = sum(Ns) / M
        gamma = sum(1 for x in Ns if x > 0) / M
        pc = (sum(x * (x - 1) for x in Ns) / M) / mu ** 2
        p_hat = sum(Ns) / sum(nr)
        pc_alg = sum(x * (x - 1) for x in Ns) / (sum(k * (k - 1) for k in nr) * p_hat ** 2)
        pg = sum(pred_gamma) / M
        nl = len(r["ells"])
        varg = sum((g - gamma) ** 2 for g in gA) / (len(gA) - 1)
        phi = varg / (gamma * (1 - gamma) / nl)
        frac_good = sum(1 for g in gA if g >= gamma / 2) / len(gA)
        both = tot = 0
        for c in cells:
            s = sum(c)
            both += s * (s - 1)
            tot += nl * (nl - 1)
        pair_ratio = (both / tot) / gamma ** 2
        print(f"  {r['cls']:5s} {r['bits']} s{r['seed']}: mu {mu:.3f} gamma {gamma:.3f} "
              f"PC {pc:.3f} PC_alg {pc_alg:.3f} [P3: 0.67..2.0] | gamma_pred {pg:.3f} "
              f"(P4 |diff|<=0.07: {abs(gamma - pg):.3f}) | overdisp {phi:.2f} [P5: 0.5..4] "
              f"Pr[g_A>=gamma/2] {frac_good:.2f} [>=0.8] pair {pair_ratio:.3f} [0.9..1.3]")


def e3(recs):
    print("E3: rank hypothesis RK")
    din = dtot = 0
    nq_in = nq_tot = q_in = q_tot = 0
    for r in sorted(recs, key=lambda x: (x["cls"], x["bits"], x["seed"])):
        valid = r["nrows"] >= 1.2 * r["ncols"]
        di = sum(r["dixon_in"])
        print(f"  {r['cls']:5s} {r['bits']} s{r['seed']}: cols {r['ncols']} rows {r['nrows']} "
              f"rank {r['rank']} valid {valid} dixon in-span {di}/{len(r['dixon_in'])}")
        if not valid:
            continue
        din += di
        dtot += len(r["dixon_in"])
        for s in r["single"]:
            if not s["qr"]:
                nq_tot += 1
                nq_in += s["inspan"]
            elif s["occ"] >= 3:
                q_tot += 1
                q_in += s["inspan"]
    print(f"  P6 Dixon vectors in span: {din}/{dtot} = {din / max(1, dtot):.3f} (predicted >= 0.8)")
    print(f"  P7a e_p, p not in Q_N, in span: {nq_in}/{nq_tot} = {nq_in / max(1, nq_tot):.3f} (<= 0.05)")
    print(f"  P7b e_p, p in Q_N, occ>=3, in span: {q_in}/{q_tot} = {q_in / max(1, q_tot):.3f} (>= 0.7)")


if __name__ == "__main__":
    which = sys.argv[1]
    recs = load(sys.argv[2:])
    dict(e1=e1, e2=e2, e3=e3)[which](recs)
