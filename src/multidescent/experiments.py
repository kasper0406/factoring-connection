"""multidescent experiments E1-E3 (notes/multidescent_prereg.md).  DIAGNOSTIC / HARNESS.

E1  children's law: labels of the children of ALL accepted level-1 points, by dyadic window and sub-bucket.
E2  rational children in a pool of fresh fields vs log-matched uniform primes: per-target pool rates p_q.
E3  algebraic children under constrained resampling (the conditional field law given the parent's value):
    joint vs product goodness, per-child concentration; E3r: rational block along the line g + j*P0.

Run:  .venv/bin/python -m src.multidescent.experiments <E1|E2|E3> <config-name> <out.jsonl> [n_seeds] [first_seed]
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from math import gcd
from typing import Dict, List

from src.common.instances import balanced_semiprime
from src.multidescent import toy

CONFIGS: Dict[str, Dict] = {
    # E1: children's law (d = 3)
    "E1": dict(d=3, n=48, M_bits=20, lb=17, S1=16, B=2 ** 11, B2=2 ** 16, ells=300, fields=2),
    # E2: rational children, pooled fresh fields
    "E2": dict(d=3, n=48, M_bits=22, lb=16, S1=16, B=2 ** 11, B2=2 ** 16, ells=400, fields=2, pool=24, S2=36,
               match=1.05),
    # E3: algebraic children, fiber law (d = 3: fiber rank 2;  d = 5: fiber rank 4)
    "E3d3": dict(d=3, n=48, M_bits=20, lb=17, S1=16, B=2 ** 12, B2=2 ** 16, S2=20, Hh=4, R=81, parents=60,
                 fields=2, line=40),
    "E3d5": dict(d=5, n=48, M_bits=10, lb=17, S1=20, B=2 ** 13, B2=2 ** 16, S2=24, Hh=2, R=81, parents=40,
                 fields=2, line=40),
}


def _primes_in(lo: int, hi: int) -> List[int]:
    return [p for p in range(lo + 1, hi + 1) if toy.is_prime(p)]


def instance(seed: int, n: int) -> int:
    return balanced_semiprime(n, seed).N


# ------------------------------------------------------------------ E1
def run_e1(cfg: Dict, seeds: List[int], out: str) -> None:
    ctx = toy.Ctx(cfg["B"], cfg["B2"])
    ell_pool = _primes_in(2 ** cfg["lb"], 2 ** (cfg["lb"] + 1))
    with open(out, "w") as fh:
        for s in seeds:
            N = instance(s, cfg["n"])
            for fi in range(cfg["fields"]):
                rng = random.Random(f"multidescent-E1|{s}|{fi}")
                A = toy.sample_field(N, rng, d=cfg["d"], M_bits=cfg["M_bits"])
                # available labels per prime (for the uniform-on-labels reference), primes in (B, B2]
                nroots = {}
                for q in ctx.primes_B2:
                    if q > ctx.B:
                        nroots[q] = len(toy.algebraic_targets(A, q))
                n_nodes = 0
                for ell in rng.sample(ell_pool, cfg["ells"]):
                    for nd in toy.level1_nodes(A, ell, cfg["S1"], ctx, max_nodes=10 ** 9):
                        n_nodes += 1
                        fh.write(json.dumps({"seed": s, "field": fi, "ell": ell, "P0": [nd["a"], nd["b"]],
                                             "rat": nd["rat"], "alg": [q for q, _ in nd["alg"]],
                                             "deg": nd["degenerate"]}) + "\n")
                fh.write(json.dumps({"seed": s, "field": fi, "summary": True, "nodes": n_nodes,
                                     "nroots": {str(q): k for q, k in nroots.items() if k}}) + "\n")
                fh.flush()


# ------------------------------------------------------------------ E2
def run_e2(cfg: Dict, seeds: List[int], out: str) -> None:
    ctx = toy.Ctx(cfg["B"], cfg["B2"])
    ell_pool = _primes_in(2 ** cfg["lb"], 2 ** (cfg["lb"] + 1))
    allp = [p for p in ctx.primes_B2 if p > ctx.B]
    with open(out, "w") as fh:
        for s in seeds:
            N = instance(s, cfg["n"])
            for fi in range(cfg["fields"]):
                rng = random.Random(f"multidescent-E2|{s}|{fi}")
                A = toy.sample_field(N, rng, d=cfg["d"], M_bits=cfg["M_bits"])
                prng = random.Random(f"multidescent-E2-pool|{s}|{fi}")
                pool = [toy.sample_field(N, prng, d=cfg["d"], M_bits=cfg["M_bits"]) for _ in range(cfg["pool"])]
                cache: Dict = {}

                def prate(q: int) -> float:
                    if q not in cache:
                        good = 0
                        for P in pool:
                            t = toy.rational_target(P, q)
                            if t is not None and toy.good_points(P, t, cfg["S2"], ctx) > 0:
                                good += 1
                        cache[q] = good / len(pool)
                    return cache[q]

                for ell in rng.sample(ell_pool, cfg["ells"]):
                    nodes = toy.level1_nodes(A, ell, cfg["S1"], ctx, max_nodes=1)
                    if not nodes or not nodes[0]["rat"]:
                        continue
                    nd = nodes[0]
                    for q in nd["rat"]:
                        lo, hi = q / cfg["match"], q * cfg["match"]
                        cands = [p for p in allp if lo <= p <= hi and p != q]
                        qu = rng.choice(cands)
                        fh.write(json.dumps({"seed": s, "field": fi, "ell": ell, "q": q, "p_child": prate(q),
                                             "qu": qu, "p_unif": prate(qu), "m_rat": len(nd["rat"])}) + "\n")
                fh.flush()


# ------------------------------------------------------------------ E3
def _all_H(d: int, Hh: int) -> List[List[int]]:
    out = [[]]
    for _ in range(d - 1):
        out = [h + [v] for h in out for v in range(-Hh, Hh + 1)]
    return out


def fiber_from_h(A: Dict, P0, h: List[int]):
    a0, b0 = P0
    Gx = [-A["g0"], A["g1"]]
    L0 = [-a0, b0]
    add = toy.poly_mul(toy.poly_mul(Gx, L0), h)
    f = list(A["f"])
    for i, v in enumerate(add):
        f[i] += v
    return toy.finish_field(A["N"], A["g0"], A["g1"], A["w"], f, {"h": h})


def run_e3(cfg: Dict, seeds: List[int], out: str) -> None:
    ctx = toy.Ctx(cfg["B"], cfg["B2"])
    d = cfg["d"]
    ell_pool = _primes_in(2 ** cfg["lb"], 2 ** (cfg["lb"] + 1))
    Hs_all = _all_H(d, cfg["Hh"])
    with open(out, "w") as fh:
        for s in seeds:
            N = instance(s, cfg["n"])
            for fi in range(cfg["fields"]):
                rng = random.Random(f"multidescent-E3|{cfg['d']}|{s}|{fi}")
                A = toy.sample_field(N, rng, d=d, M_bits=cfg["M_bits"])
                got = 0
                for ell in rng.sample(ell_pool, len(ell_pool)):
                    if got >= cfg["parents"]:
                        break
                    nodes = toy.level1_nodes(A, ell, cfg["S1"], ctx, max_nodes=1)
                    if not nodes or not nodes[0]["alg"]:
                        continue
                    nd = nodes[0]
                    P0 = (nd["a"], nd["b"])
                    F0 = toy.F_val(A, *P0)
                    Hs = Hs_all if len(Hs_all) <= cfg["R"] else rng.sample(Hs_all, cfg["R"])
                    children = nd["alg"]
                    counts = {str(c): [] for c in children}
                    # null: rational targets matched (x/1.05 .. x*1.05) to the children's norms; valid in every
                    # fiber field (g is unchanged), so they share the fields' quality but not P0's congruence
                    nulls = []
                    for (q, _) in children:
                        cands = [p for p in ctx.primes_B2 if q / 1.05 <= p <= q * 1.05 and p != q
                                 and A["g1"] % p]
                        nulls.append(rng.choice(cands))
                    ncounts = {str(p): [] for p in nulls}
                    viol = 0
                    nf = 0
                    for h in Hs:
                        B = fiber_from_h(A, P0, h)
                        if B is None:
                            continue
                        # P3a: fiber algebra
                        if toy.F_val(B, *P0) != F0 or toy.hom_eval(B["f"], A["g0"], A["g1"]) != N * A["w"]:
                            viol += 1
                        bad_label = False
                        for (q, r) in children:
                            if (B["f"][d] * B["disc"]) % q == 0:
                                bad_label = True
                            elif r not in toy.roots_mod(B["f"], q):
                                viol += 1
                        if bad_label:
                            continue          # degenerate fiber field for this parent (disclosed; skipped)
                        nf += 1
                        for (q, r) in children:
                            counts[str((q, r))].append(toy.good_points(B, ("a", q, r), cfg["S2"], ctx))
                        for p in nulls:
                            ncounts[str(p)].append(toy.good_points(B, toy.rational_target(B, p), cfg["S2"], ctx))
                    # E3r: rational block along the line g + j*P0 (G(P0) unchanged); counts of G-smooth box points
                    line = {}
                    for (q, r) in children:
                        pts = toy.box_points(q, r, cfg["S2"])
                        vals = []
                        for j in range(-cfg["line"], cfg["line"] + 1):
                            g0, g1 = A["g0"] + j * P0[0], A["g1"] + j * P0[1]
                            if gcd(g0, g1) != 1 or g1 <= 1:
                                continue
                            assert g1 * P0[0] - g0 * P0[1] == A["g1"] * P0[0] - A["g0"] * P0[1]
                            xs = [g1 * a - g0 * b for (a, b) in pts]
                            xs = [x for x in xs if x != 0]
                            sp = toy.smooth_parts(xs, ctx.PB)
                            vals.append(sum(1 for x, y in zip(xs, sp) if abs(x) == y))
                        line[str((q, r))] = vals
                    got += 1
                    fh.write(json.dumps({"seed": s, "field": fi, "ell": ell, "P0": list(P0), "d": d,
                                         "children": [list(c) for c in children], "n_fiber": nf, "viol": viol,
                                         "counts": counts, "null_counts": ncounts, "line": line}) + "\n")
                    fh.flush()


def main() -> None:
    exp, name, out = sys.argv[1], sys.argv[2], sys.argv[3]
    nseeds = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    first = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    seeds = list(range(first, first + nseeds))
    cfg = CONFIGS[name]
    t = time.time()
    {"E1": run_e1, "E2": run_e2, "E3": run_e3}[exp](cfg, seeds, out)
    print(f"{exp} {name} done in {time.time() - t:.1f}s -> {out}")


if __name__ == "__main__":
    main()
