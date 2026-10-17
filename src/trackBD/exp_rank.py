"""Experiment B3: rank-split tests for sketch-derived matrices over Z/N.

HARNESS script.  For seeded balanced semiprimes N = p q and public moduli f,
build (operationally, from N and f only)
  M1  multiplication-by-S matrix on R_N[X]/(f), S = F_N mod f
  M2  Hankel matrix of u_k = ell(X^k S mod f), ell seeded public random
  M3  Q_N : columns X^{iN} mod f ;  M3b  Q_N - I
  M5  multiplication-by-G1 matrix (G1 = (1+X)^(N-1) - X^(N-1)), Track D comparison
and run column-gcd elimination (linalg.elim_profile).  Diagnostic: ranks and pivot
profiles mod p and mod q (flint nmod_mat), predicted rank of M1 from Theorem B1
(r - deg gcd(f mod p, F_q mod p)), and classification STRUCT (profiles differ) /
ACC (factor found although profiles agree: impossible by Prop. B3, counted as check).

Usage:
  .venv/bin/python -m src.trackBD.exp_rank --bits 16 20 24 28 32 --inst 60 --procs 12 \
      --out experiments/raw/trackBD_rank.jsonl
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from multiprocessing import Pool

import flint

from src.common import instances as I
from src.trackBD.linalg import berlekamp_QN, elim_profile, hankel_proj, mult_matrix
from src.trackBD.sketch import binomial_modulus, sketch, sketch_G
from src.trackBD.exp_tiny import random_family


def moduli(seed: int):
    rng = random.Random(f"trackBD_rank_f|{seed}")
    fs = []
    for r in (2, 3, 4, 6, 8, 12, 16, 24):
        c = rng.randint(2, 50)
        fs.append((f"binom_{r}_{c}", binomial_modulus(r, c)))
    for m in (5, 7, 9, 12, 15, 21):
        fs.append((f"cyclo_{m}", [int(x) for x in flint.fmpz_poly.cyclotomic(m).coeffs()]))
    for i, f in enumerate(random_family(8, seed=777 + seed)):
        fs.append((f"rand_{i}", f))
    return fs


def rank_profile_mod(M, p):
    A = flint.nmod_mat([[x % p for x in row] for row in M], p)
    R = A.rref()[0] if isinstance(A.rref(), tuple) else A.rref()
    prof = []
    rows = R.nrows()
    for i in range(rows):
        for j in range(R.ncols()):
            if int(R[i, j]) != 0:
                prof.append(j)
                break
    return len(prof), prof


def run(args):
    bits, seed = args
    inst = I.balanced_semiprime(bits, seed)
    N = inst.N
    p, q = inst.labels["factors"]
    out = []
    for name, f in moduli(seed):
        r = len(f) - 1
        S = sketch(N, f)
        G = sketch_G(N, f, 1)
        rng = random.Random(f"trackBD_rank_ell|{N}|{name}")
        ell = [rng.randrange(N) for _ in range(r)]
        mats = {"M1": mult_matrix(N, f, S), "M2": hankel_proj(N, f, S, ell, r),
                "M3": berlekamp_QN(N, f), "M3b": berlekamp_QN(N, f, True),
                "M5": mult_matrix(N, f, G)}
        # Theorem B1 prediction for M1 rank mod p, q
        pred = {}
        for a, b in ((p, q), (q, p)):
            fa = flint.nmod_poly([x % a for x in f], a)
            Fb = (flint.nmod_poly([1, 1], a).pow_mod(b, fa) - 1 - flint.nmod_poly([0, 1], a).pow_mod(b, fa)) % fa
            Fbp = Fb.pow_mod(a, fa)
            pred[a] = r - fa.gcd(Fbp).degree() if not Fbp.is_zero() else 0
        for key, M in mats.items():
            res = elim_profile(M, N, seed=seed)
            rp, pp = rank_profile_mod(M, p)
            rq, pq_ = rank_profile_mod(M, q)
            split = pp != pq_
            rec = {"bits": bits, "seed": seed, "N": N, "f": name, "r": r, "mat": key,
                   "rank_p": rp, "rank_q": rq, "profile_split": split,
                   "factor": res.get("factor"), "acc_nonunits": res.get("acc_nonunits", 0)}
            if key == "M1":
                rec["pred_rank_p"], rec["pred_rank_q"] = pred[p], pred[q]
            if res.get("factor"):
                rec["class"] = "STRUCT" if split else "ACC"
            out.append(rec)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="+", default=[16, 20, 24, 28, 32])
    ap.add_argument("--inst", type=int, default=60)
    ap.add_argument("--procs", type=int, default=12)
    ap.add_argument("--out", default="experiments/raw/trackBD_rank.jsonl")
    a = ap.parse_args()
    jobs = [(b, s) for b in a.bits for s in range(a.inst)]
    t0 = time.time()
    with Pool(a.procs) as pool, open(a.out, "w") as fh:
        for i, recs in enumerate(pool.imap_unordered(run, jobs)):
            for rec in recs:
                fh.write(json.dumps(rec) + "\n")
            if i % 50 == 0:
                print(f"{i}/{len(jobs)} {time.time() - t0:.0f}s", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
