"""E1 harness (DIAGNOSTIC): instrumented one-layer N-exponent trees.  Pre-registration: notes/leakmap_prereg.md E1.

Usage:  .venv/bin/python -m src.leakmap.e1 [--quick]
Writes experiments/raw/leakmap_e1.jsonl (one record per (n, population, tree)) and
experiments/raw/leakmap_e1_p4.jsonl (random-base liar counts).
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from collections import Counter, defaultdict
from math import gcd
from typing import Dict, List, Tuple

import gmpy2

from . import netrees as nt


def primes_top2(b: int) -> List[int]:
    lo, hi = 3 << (b - 2), 1 << b
    out, x = [], int(gmpy2.next_prime(lo - 1))
    while x < hi:
        out.append(x)
        x = int(gmpy2.next_prime(x))
    return out


def population(b: int, kind: str, sample: int = 0, seed: str = "") -> List[Tuple[int, int]]:
    P = primes_top2(b)
    if kind == "iid":
        pairs = [(P[i], P[j]) for i in range(len(P)) for j in range(i + 1, len(P))]
        if sample and sample < len(pairs):
            rng = random.Random(seed)
            pairs = rng.sample(pairs, sample)
        return pairs
    if kind == "fermat":
        H = 1 << (b // 2)
        return [(P[i], P[j]) for i in range(len(P)) for j in range(i + 1, len(P)) if P[j] - P[i] <= H]
    raise ValueError(kind)


def targets(p: int, q: int, b: int) -> Dict[str, int]:
    return {"eps": int(p % 4 == 1), "p1mod3": int(p % 3 == 1), "arch": int(q - p < (1 << (b - 4)))}


def run_cell(b: int, kind: str, tname: str, pairs: List[Tuple[int, int]]) -> Dict:
    tree = nt.TREES[tname]()
    L = 1 << (b // 2)
    kappa: Dict[int, tuple] = {}
    p1_mismatch = 0
    nsplit = 0
    split_struct = 0
    split_acc = 0
    dev_count = Counter()          # deviations at x (x as a member of a pair)
    partners: Dict[int, List[int]] = defaultdict(list)
    dclass = Counter()
    examples = []
    nonsplit_dev = 0
    groups: Dict[str, Dict[tuple, Counter]] = {k: defaultdict(Counter) for k in ("eps", "p1mod3", "arch")}
    tcount = {k: Counter() for k in groups}
    for (p, q) in pairs:
        N = p * q
        partners[p].append(q)
        partners[q].append(p)
        out = nt.run_on_N(tree, N)
        tp = nt.diag_transcript(tree, N, p)
        tq = nt.diag_transcript(tree, N, q)
        for x in (p, q):
            if x not in kappa:
                kappa[x] = nt.formal_transcript(tree, x)
        dp, dq = tp != kappa[p], tq != kappa[q]
        if dp:
            dev_count[p] += 1
        if dq:
            dev_count[q] += 1
        # P1
        if out[0] == "split":
            nsplit += 1
            if tp == tq or out[1] not in (p, q):
                p1_mismatch += 1
            # STRUCTURAL: explained by a transcript disagreement (always, by the lemma) ; ACCIDENTAL cannot occur
            # for element tests (gcd of all coordinates), counted for rule 3.
            if tp != tq:
                split_struct += 1
            else:
                split_acc += 1
        else:
            if not (tp == tq == out[1]):
                p1_mismatch += 1
        ns = out[0] != "split"
        if ns and (dp or dq):
            nonsplit_dev += 1
        # classify the first deviating node at each side
        for x, tx, dx in ((p, tp, dp), (q, tq, dq)):
            if not dx:
                continue
            kx = kappa[x]
            nodes = nt.path_nodes(tree, kx)
            for i, nd in enumerate(nodes):
                if i >= len(tx) or tx[i] != kx[i]:
                    if isinstance(nd, nt.GateNode):
                        cls = nt.deviation_class(nd, N, x, L)
                        dclass[(cls, "nonsplit" if ns else "split")] += 1
                        if ns and len(examples) < 8:
                            examples.append(dict(p=p, q=q, side=x, node=nd.label, cls=cls,
                                                 orders=[c.order for c in nt.node_components(nd, x)]))
                    else:
                        dclass[("PURE?", "")] += 1
                    break
        tg = targets(p, q, b)
        for k, v in tg.items():
            tcount[k][v] += 1
            if ns:
                groups[k][out[1]][v] += 1
    tot = len(pairs)
    s = nsplit / tot
    pdev_ns = nonsplit_dev / tot
    p2 = {}
    for k in groups:
        base = max(tcount[k].values()) / tot
        corr_ns = sum(max(c.values()) for c in groups[k].values())
        success = (nsplit + corr_ns) / tot
        bound = base + 2 * s + 2 * pdev_ns
        gain = corr_ns / tot - (1 - s) * base
        p2[k] = dict(base=base, success=success, bound=bound, ok=success <= bound + 1e-12, nonsplit_gain=gain)
    # P3: per-x deviation frequency vs NE1 counting bound with exact mu over that x's partners
    p3_viol, p3_worst, p3_n = 0, 0.0, 0
    for x, ys in partners.items():
        kx = kappa[x]
        bound = 0.0
        for nd in nt.path_nodes(tree, kx):
            if not isinstance(nd, nt.GateNode) or nd.formal_zero(x):
                continue
            for comp in nt.node_components(nd, x):
                cnt = Counter(nd.E(x * y) % comp.order for y in ys)
                bound += max(cnt.values()) / len(ys)
        freq = dev_count[x] / len(ys)
        p3_n += 1
        if freq > bound + 1e-12:
            p3_viol += 1
        if bound > 0:
            p3_worst = max(p3_worst, freq / bound)
    return dict(b=b, n=2 * b, population=kind, tree=tname, pairs=tot, p1_mismatch=p1_mismatch, split_rate=s,
                splits_structural=split_struct, splits_accidental=split_acc,
                nonsplit_dev_rate=pdev_ns, dev_rate_side=sum(dev_count.values()) / (2 * tot),
                dev_classes={f"{a}|{c}": v for (a, c), v in sorted(dclass.items())},
                p2=p2, nonsplit_examples=examples, p3_primes=p3_n, p3_violations=p3_viol, p3_worst_ratio=p3_worst)


def run_p4(b: int, npairs: int, nbases: int, seed: str) -> Dict:
    """Random-base Fermat and MR liars (two-sided deviations) vs the order-pattern prediction."""
    rng = random.Random(seed)
    pairs = population(b, "iid")
    if npairs and npairs < len(pairs):
        pairs = rng.sample(pairs, npairs)
    obs_f = obs_mr = 0
    exp_f = exp_mr = 0.0
    for (p, q) in pairs:
        N = p * q
        phi = (p - 1) * (q - 1)
        g = gcd(p - 1, q - 1)
        h = g
        while h % 2 == 0:
            h //= 2
        k = min(((p - 1) & -(p - 1)).bit_length() - 1, ((q - 1) & -(q - 1)).bit_length() - 1)
        exp_f += nbases * g * g / phi
        exp_mr += nbases * h * h * (4 ** k + 2) / (3 * phi)
        t = ((N - 1) & -(N - 1)).bit_length() - 1
        m = (N - 1) >> t
        for _ in range(nbases):
            a = rng.randrange(1, N)
            if gcd(a, N) != 1:
                continue
            if pow(a, N - 1, N) == 1:
                obs_f += 1
            x = pow(a, m, N)
            ok = x == 1 or x == N - 1
            j = 0
            while not ok and j < t - 1:
                x = x * x % N
                j += 1
                if x == N - 1:
                    ok = True
            if ok:
                obs_mr += 1
    zf = (obs_f - exp_f) / max(exp_f, 1e-9) ** 0.5
    zm = (obs_mr - exp_mr) / max(exp_mr, 1e-9) ** 0.5
    return dict(b=b, n=2 * b, pairs=len(pairs), bases=nbases, fermat_obs=obs_f, fermat_exp=exp_f, fermat_z=zf,
                mr_obs=obs_mr, mr_exp=exp_mr, mr_z=zm, ok=abs(zf) <= 4 and abs(zm) <= 4)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default="experiments/raw/leakmap_e1.jsonl")
    ap.add_argument("--out4", default="experiments/raw/leakmap_e1_p4.jsonl")
    a = ap.parse_args(argv)
    cells = [(12, "iid", 0), (14, "iid", 0), (16, "iid", 200000), (12, "fermat", 0), (14, "fermat", 0),
             (16, "fermat", 0)]
    if a.quick:
        cells = [(10, "iid", 0), (10, "fermat", 0)]
    with open(a.out, "w") as f:
        for (b, kind, sample) in cells:
            pairs = population(b, kind, sample, seed=f"leakmapE1|{2 * b}")
            for tname in nt.TREES:
                t0 = time.time()
                rec = run_cell(b, kind, tname, pairs)
                rec["secs"] = round(time.time() - t0, 2)
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print(json.dumps({k: rec[k] for k in ("n", "population", "tree", "pairs", "p1_mismatch",
                                                      "split_rate", "nonsplit_dev_rate", "p3_violations",
                                                      "dev_classes", "secs")}), flush=True)
    with open(a.out4, "w") as f:
        p4cells = [(14, 20000, 20, "leakmapP4|28")] if not a.quick else [(10, 0, 50, "q")]
        if not a.quick:
            p4cells.append((10, 0, 2000, "leakmapP4|20-deviation"))
        for (b, npairs, nb, sd) in p4cells:
            rec = run_p4(b, npairs, nb, sd)
            f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
