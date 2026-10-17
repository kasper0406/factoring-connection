"""Experiment B2c: exact expectation over the hidden cofactor of the single-sketch
(scalar) vanishing count, for fixed p.

HARNESS/diagnostic.  For a prime p and every unit e mod (p-1) (the residue of the
other prime q mod p-1; uniform over units for random q by Dirichlet), count
x in F_p \\ U (U = {0,-1} and primitive cube roots) with (1+x)^e = 1 + x^e, split into
   FROB : x^(e-1) = 1 = (1+x)^(e-1)   (Frobenius coincidence, e == 1 mod L_x)
   SPOR : the rest.
Also G_1 roots (g - 1 with g = gcd(p-1, e-1)) minus the public root -1/2 when present.
Outputs mean, median, P(count > 0) over e.  Cost O(p * phi(p-1)).

Usage: .venv/bin/python -m src.trackBD.exp_expect --out experiments/raw/trackBD_expect.jsonl
"""
from __future__ import annotations

import argparse
import json
import random
from math import gcd
from multiprocessing import Pool

from src.common import instances as I


def run(p):
    U = {0, p - 1} | {x for x in range(p) if (x * x + x + 1) % p == 0}
    xs = [x for x in range(p) if x not in U]
    fr, sp, g1 = [], [], []
    for e in range(1, p - 1):
        if gcd(e, p - 1) != 1:
            continue
        f = s = 0
        for x in xs:
            if (pow(1 + x, e, p) - 1 - pow(x, e, p)) % p == 0:
                if pow(x, e - 1, p) == 1 and pow(1 + x, e - 1, p) == 1:
                    f += 1
                else:
                    s += 1
        fr.append(f)
        sp.append(s)
        g1.append(gcd(p - 1, e - 1) - 2)  # nonpublic G_1 roots (e-1 even => -1/2 is a root)

    def st(v):
        v2 = sorted(v)
        return {"mean": sum(v) / len(v), "median": v2[len(v) // 2], "P_pos": sum(1 for a in v if a) / len(v),
                "max": v2[-1]}
    return {"p": p, "n_e": len(fr), "frob": st(fr), "spor": st(sp), "G1_nonpublic": st(g1),
            "spor_mod6_nonzero": sum(1 for a in sp if a % 6)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="+", default=[8, 10, 11])
    ap.add_argument("--count", type=int, default=12)
    ap.add_argument("--procs", type=int, default=12)
    ap.add_argument("--out", default="experiments/raw/trackBD_expect.jsonl")
    a = ap.parse_args()
    rng = random.Random("trackBD_expect")
    ps = sorted({I.random_prime(rng, b) for b in a.bits for _ in range(a.count)})
    with Pool(a.procs) as pool, open(a.out, "w") as fh:
        for rec in pool.imap_unordered(run, ps):
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
