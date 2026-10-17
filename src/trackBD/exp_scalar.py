"""Experiment B2/D2: degree-1 sketches (scalar F_N(c), G_k(c)) on seeded balanced
semiprimes and diagnostic families; measured rates vs the exact criterion.

HARNESS script.  Operational part per instance: draw K public random c in
[2, N-2] from rng("trackBD_scalar|N|seed"), compute
   F_N(c) = (1+c)^N - 1 - c^N,  G_k(c) = (1+c)^(N-k) - c^(N-k)  (k=1,2)
and a null residue u (uniform mod N), then gcd with N.

Diagnostic part: classify each hit and run the exact census of roots over F_p
(diagnose.scalar_roots_census-like, here inlined & extended) for p, q <= CENSUS_MAX.

Usage:
  .venv/bin/python -m src.trackBD.exp_scalar --bits 16 20 24 28 32 36 40 --inst 200 --K 4000 \
     --procs 16 --out experiments/raw/trackBD_scalar.jsonl
  .venv/bin/python -m src.trackBD.exp_scalar --family diag_bigg --bits 24 32 --inst 50 ...
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from collections import Counter
from math import gcd
from multiprocessing import Pool

from src.common import instances as I
from src.trackBD.sketch import scalar_F, scalar_G

CENSUS_MAX = 1 << 20


def census(p: int, m: int, k_list=(1, 2)) -> dict:
    """Exact counts over c in F_p.  F-roots split into
    triv  : c in {0,-1} or c^2+c+1 = 0 (universal Cauchy-Mirimanoff roots),
    frob  : c^(m-1) = 1 = (1+c)^(m-1)    (Frobenius coincidence),
    spor  : the rest (sporadic).
    G_k roots: (1+c)^(N-k) == c^(N-k) mod p, N = p m, exponent N-k == m-k (mod p-1)
    for c not in {0,-1}."""
    triv = frob = spor = 0
    Gk = {k: 0 for k in k_list}
    half = (-pow(2, -1, p)) % p
    G1_half = 0
    for c in range(p):
        a = pow(1 + c, m, p)
        b = (1 + pow(c, m, p)) % p
        if a == b:
            if c in (0, p - 1) or (c * c + c + 1) % p == 0:
                triv += 1
            elif pow(c, m - 1, p) == 1 and pow(1 + c, m - 1, p) == 1:
                frob += 1
            else:
                spor += 1
        if c not in (0, p - 1):
            for k in k_list:
                e = (p * m - k) % (p - 1)
                if pow(1 + c, e, p) == pow(c, e, p):
                    Gk[k] += 1
                    if k == 1 and c == half:
                        G1_half += 1
    return {"triv": triv, "frob": frob, "spor": spor,
            **{f"G{k}": v for k, v in Gk.items()}, "G1_half": G1_half,
            "g": gcd(p - 1, m - 1)}


def classify_F(p: int, m: int, c: int) -> str:
    """TD_x: p divides the public integer c (x=0), c+1 (x=m1), c^2+c+1 (x=phi3)."""
    c %= p
    if c == 0:
        return "TD_0"
    if c == p - 1:
        return "TD_m1"
    if (c * c + c + 1) % p == 0:
        return "TD_phi3"
    if pow(c, m - 1, p) == 1 and pow(1 + c, m - 1, p) == 1:
        return "FROB"
    return "SPOR"


def make_diag_bigg(bits: int, seed: int, s: int) -> I.Instance:
    """DIAGNOSTIC family: q = 1 + t (p-1)/s prime, so gcd(p-1, q-1) >= (p-1)/s."""
    rng = random.Random(f"trackBD_diag_bigg|{bits}|{s}|{seed}")
    while True:
        p = I.random_prime(rng, bits // 2)
        step = (p - 1) // s
        lo = (1 << (bits - bits // 2 - 1)) // step + 1
        for _ in range(2000):
            t = rng.randint(lo, 2 * lo)
            q = 1 + t * step
            if q != p and I.is_prime(q):
                N = p * q
                return I._mk(N, "diag_bigg", {"bits": bits, "s": s}, seed, [p, q], q=q)


def run_instance(args):
    fam, bits, seed, K, extra = args
    t0 = time.perf_counter()
    if fam == "balanced":
        inst = I.balanced_semiprime(bits, seed)
    elif fam == "diag_bigg":
        inst = make_diag_bigg(bits, seed, extra)
    elif fam == "carmichael":
        inst = I.carmichael_small(seed)
    else:
        raise ValueError(fam)
    N = inst.N
    primes = sorted(set(inst.factors))
    # ---- operational -------------------------------------------------
    rng = random.Random(f"trackBD_scalar|{N}|{seed}")
    hits = Counter()
    for _ in range(K):
        c = rng.randint(2, N - 2)
        u = rng.randrange(N)
        vals = {"F": scalar_F(N, c), "G1": scalar_G(N, c, 1), "G2": scalar_G(N, c, 2), "NULL": u}
        for name, v in vals.items():
            d = gcd(v, N)
            if 1 < d < N:
                # ---- diagnostic classification --------------------------
                for p in primes:
                    if d % p == 0:
                        m = N // p
                        if name == "F":
                            lab = classify_F(p, m, c)
                        elif name.startswith("G"):
                            cp = c % p
                            lab = ("TD" if cp in (0, p - 1) else
                                   "TD_half" if (2 * cp + 1) % p == 0 else "ORDER")
                        else:
                            lab = "NULL"
                        side = "p" if p == primes[0] else "q"
                        hits[f"{name}|{lab}|{side}"] += 1
                hits[f"{name}|any"] += 1
    rec = {"family": fam, "bits": bits, "seed": seed, "N": N, "extra": extra, "K": K,
           "factors": inst.factors, "hits": dict(hits)}
    if len(primes) == 2 and primes[1] <= CENSUS_MAX:
        p, q = primes
        rec["census"] = {"p": census(p, N // p), "q": census(q, N // q)}
    rec["wall_s"] = round(time.perf_counter() - t0, 3)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", default="balanced")
    ap.add_argument("--bits", type=int, nargs="+", default=[16, 20, 24, 28, 32, 36, 40])
    ap.add_argument("--inst", type=int, default=200)
    ap.add_argument("--K", type=int, default=4000)
    ap.add_argument("--s", type=int, nargs="+", default=[0])
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--out", default="experiments/raw/trackBD_scalar.jsonl")
    a = ap.parse_args()
    jobs = [(a.family, b, sd, a.K, s) for b in a.bits for s in a.s for sd in range(a.inst)]
    if a.family == "carmichael":
        jobs = [("carmichael", 0, i, a.K, 0) for i in range(a.inst)]
    t0 = time.time()
    with Pool(a.procs) as pool, open(a.out, "a") as fh:
        for i, rec in enumerate(pool.imap_unordered(run_instance, jobs, chunksize=1)):
            fh.write(json.dumps(rec) + "\n")
            if i % 100 == 0:
                print(f"{i}/{len(jobs)} {time.time() - t0:.0f}s", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
