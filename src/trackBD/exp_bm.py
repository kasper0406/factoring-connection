"""Experiment D5: Berlekamp-Massey over Z/N on projected sketch sequences.

HARNESS script.  Operational sequences (from N only):
  S1  s_c = F_N(c),                         c = 0..T-1     (weights c, r = 1)
  S2  s_c = [X^1] F_N mod (X^2 - c),        c = 1..T
  S3  s_m = [X^0] F_N mod (X^m - 1),        m = 1..T
BM over Z/N (linalg.bm_mod_N) stops at the first nonzero nonunit discrepancy.
Diagnostic: linear complexities over F_p, F_q (linalg.bm_field) and the event class
  STRUCT  the connection polynomial at the event annihilates the whole available
          sequence mod p (the p-side linear complexity has been reached)
  ACC     a transient discrepancy vanishing mod p.
Prediction for S1 (Prop. D5 in notes): LC_p = e_p := (q mod (p-1)), LC_q = p
(for q > p), so a STRUCT event needs T >~ 2 min(e_p, p).

Usage: .venv/bin/python -m src.trackBD.exp_bm --bits 16 20 24 --inst 30 --procs 8
"""
from __future__ import annotations

import argparse
import json
import sys
from multiprocessing import Pool

from src.common import instances as I
from src.trackBD.linalg import bm_field, bm_mod_N
from src.trackBD.sketch import binomial_modulus, scalar_F, sketch


def seqs(N: int, T: int):
    out = {}
    out["S1"] = [scalar_F(N, c) for c in range(T)]
    out["S2"] = [sketch(N, binomial_modulus(2, c))[1] for c in range(1, T + 1)]
    out["S3"] = [sketch(N, binomial_modulus(m, 1))[0] for m in range(1, T + 1)]
    return out


def persistent(C, seq, p, upto):
    L = len(C) - 1
    for n in range(L, upto):
        v = sum(C[i] * seq[n - i] for i in range(len(C))) % p
        if v:
            return False
    return True


def run(args):
    bits, seed, T = args
    inst = I.balanced_semiprime(bits, seed)
    N = inst.N
    p, q = inst.labels["factors"]
    recs = []
    for name, s in seqs(N, T).items():
        res = bm_mod_N(s, N)
        lp = bm_field(s, p)["L"]
        lq = bm_field(s, q)["L"]
        rec = {"bits": bits, "seed": seed, "N": N, "p": p, "q": q, "seq": name, "T": T,
               "LC_p": lp, "LC_q": lq, "e_p": q % (p - 1), "factor": res.get("factor"),
               "step": res.get("step")}
        if res.get("factor"):
            d = res["factor"]
            ell = p if d % p == 0 else q
            C = [x % ell for x in res["C"]]
            rec["revealed"] = "p" if ell == p else "q"
            rec["class"] = "STRUCT" if persistent(C, s, ell, T) else "ACC"
        recs.append(rec)
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="+", default=[16, 20, 24])
    ap.add_argument("--inst", type=int, default=30)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--out", default="experiments/raw/trackBD_bm.jsonl")
    a = ap.parse_args()
    Tmap = {16: 768, 20: 2560, 24: 4096}
    jobs = [(b, s, Tmap.get(b, 2048)) for b in a.bits for s in range(a.inst)]
    with Pool(a.procs) as pool, open(a.out, "w") as fh:
        for i, recs in enumerate(pool.imap_unordered(run, jobs)):
            for r in recs:
                fh.write(json.dumps(r) + "\n")
            print(f"{i}/{len(jobs)}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
