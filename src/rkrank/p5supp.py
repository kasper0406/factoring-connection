"""HARNESS (supplementary, disclosed deviation; notes/rkrank_results.md §5).

The frozen harness computed the P5 functionals only at level 3.0 (`structured = top` with top = max of the
configured levels), whereas the pre-registration says "at the top level reached".  This script recomputes,
for every instance of experiments/configs/rkrank_e1.json, the identical relation set (same seeds), and
evaluates the registered P5 functionals at the top level actually reached.  It also checks the independence
assumption behind excess_peel: rank(Gen ∪ {e_p : p peeled}) == omega + n_peel.
Usage: .venv/bin/python -m src.rkrank.p5supp CONFIG.json OUT.jsonl [workers]
"""
from __future__ import annotations

import json
import random
import sys
from multiprocessing import Pool

from ..nfsfruit.instances import composite
from . import gf2
from .harness import peel_core
from .toynfs import ToyField, gen_functionals


def job(args):
    cls, bits, seed, cfg = args
    N, fs = composite(cls, bits, seed)
    p = cfg["params"][str(bits)]
    T = ToyField(N, p["B"], cfg["nchar"], seed, cfg["H"])
    target = int(max(cfg["levels"]) * T.ncols) + 1
    rows, pairs, blast = T.sieve(p["A"], p["bmax"], target)
    rng = random.Random(f"rkrank-e1|{cls}|{bits}|{seed}")
    order = list(range(len(rows)))
    rng.shuffle(order)
    rows = [rows[i] for i in order]
    gen = gen_functionals(T, fs)
    reached = [L for L in cfg["levels"] if int(L * T.ncols) <= len(rows)]
    D = T.D_rat
    ideal_mask = 0
    for key, c in T.col.items():
        if key[0] == "a":
            ideal_mask |= 1 << c
    out = dict(cls=cls, bits=bits, seed=seed, top=max(reached), indep=[])
    for L in reached:
        sub = rows[:int(L * T.ncols)]
        core = peel_core(sub, ideal_mask)
        core_rat = 0
        for i in core:
            core_rat |= sub[i]
        peeled = [1 << T.col[("r", q)] for q in T.ps if not (core_rat >> T.col[("r", q)]) & 1]
        out["indep"].append((L, gf2.rank(gen + peeled) == len(fs) + len(peeled)))
        if L == max(reached):
            # identical functional construction to harness.level_record (same rng stream position:
            # the frozen harness drew the dense/sparse functionals from rng after the shuffle; replicate)
            def delta(phi):
                return sum(gf2.dot(r, phi) for r in sub) / len(sub)
            rat_primes = (1 << D) - 2
            dense = [rng.getrandbits(T.ncols) for _ in range(10)]
            out["delta"] = dict(liouville_joint=delta(rat_primes | ideal_mask),
                                liouville_rat=delta(rat_primes), liouville_alg=delta(ideal_mask),
                                dense=[delta(x) for x in dense])
    return out


def main():
    cfg = json.load(open(sys.argv[1]))
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    jobs = [(c, b, s, cfg) for c in cfg["classes"] for b in cfg["bits"] for s in cfg["seeds"]]
    with Pool(workers) as pool, open(sys.argv[2], "w") as out:
        for rec in pool.imap_unordered(job, jobs):
            out.write(json.dumps(rec) + "\n")
            out.flush()


if __name__ == "__main__":
    main()
