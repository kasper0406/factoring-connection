"""HARNESS: Track A (d) selector experiment on balanced semiprimes.

For each instance and selector, run `descend` (operational, sees only N) and
score every level with diagnostic truth (Lucas support, uses the factors):
  parent_nonzero, child0_nonzero, child1_nonzero, chosen_nonzero.
Also classify every gcd discovery made by the selector:
  STRUCTURAL  if the evaluated section class is hit by the support of exactly
              one prime (so the value is 0 mod the other prime), else ACCIDENTAL.

Usage:
  .venv/bin/python -m src.trackA.run_descent --bits 16 20 24 28 32 --seeds 5 \
      --max-rank 1048576 --out experiments/raw/trackA_descent.jsonl
"""
from __future__ import annotations

import argparse
import json
import time

from ..common import instances as I
from ..common.metrics import Metrics, append_jsonl
from .diagnostic import section_support
from .selectors import BudgetSelector, ExactScalarSelector, RandomSelector, descend


def pure_class(N, fs, r, M):
    sup = section_support(N, fs, r, M)
    owners = {l for i in sup for l in fs if i % l == 0}
    return len(owners) == 1  # nonzero and divisible by exactly one prime component


def run(bits_list, seeds, max_rank, out, selectors):
    for bits in bits_list:
        for seed in range(seeds):
            inst = I.balanced_semiprime(bits, seed)
            N, fs = inst.N, inst.factors
            p, q = fs
            for sel in selectors:
                m = Metrics()
                max_level = min(N.bit_length() + 1, max_rank.bit_length() - 1)
                t0 = time.perf_counter()
                r, M, logs, disc = descend(N, sel, seed=seed, max_level=max_level, m=m)
                wall = time.perf_counter() - t0
                rows = []
                rr = 0
                MM = 1
                alive = True
                for lg in logs:
                    c0 = section_support(N, fs, rr, 2 * MM)
                    c1 = section_support(N, fs, rr + MM, 2 * MM)
                    chosen = c1 if lg.bit else c0
                    rows.append({"level": lg.level, "M": lg.M, "bit": lg.bit,
                                 "parent_nonzero": alive,
                                 "child0_nonzero": bool(c0), "child1_nonzero": bool(c1),
                                 "chosen_nonzero": bool(chosen), "certified": lg.certified,
                                 "evals": lg.evals, "rank": lg.rank,
                                 "cum_q_mul_dimsum": lg.q_mul_dimsum})
                    alive = alive and bool(chosen)
                    rr += lg.bit * MM
                    MM *= 2
                for d in disc:
                    # STRUCTURAL iff the evaluated class is pure (hit by one prime's support)
                    d.kind = ("STRUCTURAL" if d.modulus and pure_class(N, fs, d.residue, d.modulus)
                              else "ACCIDENTAL")
                rec = m.to_record(experiment="trackA_descent", N=N, bits=bits, seed=seed,
                                  selector=sel.name, p=p, q=q, levels=rows,
                                  final_r=r, final_M=M, path_alive=alive,
                                  discoveries=[d.__dict__ for d in disc], wall_s=wall)
                append_jsonl(out, rec)
                first_dead = next((x["level"] for x in rows if not x["chosen_nonzero"]), None)
                print(f"bits={bits} seed={seed} {sel.name:14s} levels={len(rows)} alive={alive} "
                      f"first_dead={first_dead} log2p={p.bit_length()} disc={len(disc)} "
                      f"cost={m.get('q_mul_dimsum')} wall={wall:.2f}s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="+", default=[16, 20, 24])
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--max-rank", type=int, default=1 << 20)
    ap.add_argument("--out", default="experiments/raw/trackA_descent.jsonl")
    ap.add_argument("--selectors", nargs="+", default=["random", "exact", "budget1024"])
    a = ap.parse_args()
    sels = []
    for s in a.selectors:
        if s == "random":
            sels.append(RandomSelector())
        elif s == "exact":
            sels.append(ExactScalarSelector())
        elif s.startswith("budget"):
            sels.append(BudgetSelector(int(s[6:])))
    run(a.bits, a.seeds, a.max_rank, a.out, sels)


if __name__ == "__main__":
    main()
