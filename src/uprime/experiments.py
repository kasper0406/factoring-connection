"""uprime E2: mixed-smoothness level-1 acceptance and the conditional child mean (notes/uprime_prereg.md).
DIAGNOSTIC / HARNESS.

For each seed (48-bit balanced semiprime) and field (design (a'''), d = 3) and each sampled special-q prime ell:
  - level-1 box (reduced basis, radius S1): counts rat / alg2 / both / acc and the removal reasons (P2a, P2b);
  - the FIRST accepted point in box order (the search order), and each of its algebraic children (q, r):
      child count N_q in the square box of radius V2 = floor(sqrt(q) S2) of {a = r b mod q}   (conditional law),
      and a null: a fresh field A' from a per-seed pool, a uniformly random simple label (q', r') of A' with
      q' in [q/1.05, 1.05 q], the same box rule                                                (product main term);
    plus the irregularity flag 4 S2^2 lambda1^2 < q (P2d-P2f).

Run: .venv/bin/python -m src.uprime.experiments <out.jsonl> [n_seeds] [first_seed]
"""
from __future__ import annotations

import json
import random
import sys
import time
from typing import Dict, List

from src.common.instances import balanced_semiprime
from src.uprime import toy

CONFIG = dict(d=3, n=48, M_bits=20, lb=17, S1=16, B=2 ** 11, B2=2 ** 16, S2=24, ells=300, fields=2, pool=4,
              match=1.05)


def run(out: str, seeds: List[int], cfg: Dict = CONFIG) -> None:
    ctx = toy.Ctx(cfg["B"], cfg["B2"])
    ell_pool = [p for p in range(2 ** cfg["lb"] + 1, 2 ** (cfg["lb"] + 1)) if toy.is_prime(p)]
    with open(out, "w") as fh:
        fh.write(json.dumps({"header": True, "cfg": cfg, "seeds": seeds}) + "\n")
        for s in seeds:
            N = balanced_semiprime(cfg["n"], s).N
            prng = random.Random(f"uprime-E2-pool|{s}")
            pool = [toy.sample_field(N, prng, d=cfg["d"], M_bits=cfg["M_bits"]) for _ in range(cfg["pool"])]
            for fi in range(cfg["fields"]):
                rng = random.Random(f"uprime-E2|{s}|{fi}")
                A = toy.sample_field(N, rng, d=cfg["d"], M_bits=cfg["M_bits"])
                t0 = time.time()
                for ell in rng.sample(ell_pool, cfg["ells"]):
                    res = toy.classify_level1(A, ell, cfg["S1"], ctx)
                    nodes = res.pop("nodes", [])
                    fh.write(json.dumps({"kind": "box", "seed": s, "field": fi, "ell": ell, **res}) + "\n")
                    if not nodes:
                        continue
                    nd = nodes[0]
                    for (q, r) in nd["children"]:
                        cc = toy.child_count(A, q, r, cfg["S2"], ctx)
                        k = rng.randrange(len(pool))
                        lab = toy.random_label(pool[k], q / cfg["match"], q * cfg["match"], rng, ctx)
                        nul = toy.child_count(pool[k], lab[0], lab[1], cfg["S2"], ctx) if lab else None
                        fh.write(json.dumps({"kind": "child", "seed": s, "field": fi, "ell": ell,
                                             "P0": [nd["a"], nd["b"]], "q": q, "m": len(nd["children"]),
                                             "cond": cc, "null": nul, "null_q": lab[0] if lab else None,
                                             "null_pool": k}) + "\n")
                fh.flush()
                print(f"seed {s} field {fi} done in {time.time() - t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    run(out, list(range(first, first + n)))
