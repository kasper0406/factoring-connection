"""HARNESS for hypD experiments E1-E3 (pre-registered in notes/hypD_prereg.md).

Usage:  .venv/bin/python -m src.hypD.harness {e1,e2,e3} CONFIG.json OUT.jsonl [workers]
Instances come from src/nfsfruit/instances.composite (seeded, harness-only factor lists).
"""
from __future__ import annotations

import json
import random
import sys
from math import isqrt
from multiprocessing import Pool
from typing import Dict

import gmpy2

from ..nfsfruit.instances import composite
from ..nfsfruit.smooth import FactorBase, primes_upto
from . import split as sp
from . import specialq as sq
from . import rank as rk


# ------------------------------------------------------------------ E1: Theorem S
def e1_job(args) -> Dict:
    cls, bits, seed, k, g, trials = args
    N, fs = composite(cls, bits, seed)
    rng = random.Random(f"hypD-e1|{cls}|{bits}|{seed}|{k}|{g}")
    Y = isqrt(N // 2)                                   # Y^2 < N
    lo = 1 << ((Y.bit_length() - 1) // k)               # J = (lo, 2 lo]
    J = [p for p in primes_upto(2 * lo) if p > lo]
    G = sorted(p for p in J if rng.random() < g)
    member = sp.k_prime_member(G, k)
    inS = lambda u: (Y // 2 <= u <= Y) and member(u)
    S = sp.count_S_window(G, k, Y)
    res = sp.sample_rate(N, Y, inS, trials, seed)
    phiN = 1
    for p in fs:
        phiN *= p - 1
    model = S * S / phiN                                # heuristic: coprime pairs, chars cancel
    bound = sp.theorem_s_bound(S, Y, N, len(fs))
    return dict(exp="e1", cls=cls, bits=bits, seed=seed, k=k, g=g, N=N, omega=len(fs),
                Y=Y, J=[lo, 2 * lo], G_size=len(G), J_size=len(J), S=S, model=model,
                bound=bound, **res)


# ------------------------------------------------------------------ E2: special-q counts
def e2_job(args) -> Dict:
    cls, bits, seed, cfg = args
    N, fs = composite(cls, bits, seed)
    fb = FactorBase(cfg["B"])
    rng = random.Random(f"hypD-e2|{cls}|{bits}|{seed}")
    lo = 1 << cfg["lbits"]
    J = [p for p in primes_upto(2 * lo) if p > lo]
    ells = sorted(rng.sample(J, cfg["n_ell"]))
    fields = []
    for fi in range(cfg["n_fields"]):
        m = sq.random_m(N, rng)
        c = sq.base_m_digits(N, m)
        H = cfg["H"]
        r0, r1 = rng.randrange(-H, H + 1), rng.randrange(-H, H + 1)
        nrat, counts = [], []
        for l in ells:
            pts = sq.rational_smooth_points(l, m, cfg["S"], fb)
            nrat.append(len(pts))
            counts.append(sq.count_Nl(N, c, m, r0, r1, pts, fb))
        fields.append(dict(m=m, r0=r0, r1=r1, nrat=nrat, N_l=counts))
    return dict(exp="e2", cls=cls, bits=bits, seed=seed, N=N, ells=ells, cfg=cfg, fields=fields)


# ------------------------------------------------------------------ E3: rank hypothesis RK
def e3_job(args) -> Dict:
    cls, bits, seed, cfg = args
    N, fs = composite(cls, bits, seed)
    T = rk.ToyNFS(N, cfg["B"], cfg["nchar"], seed, cfg["H"])
    target = int(cfg["oversample"] * T.ncols)
    T.sieve(cfg["A"], cfg["Bb"], target)
    RS = rk.RowSpace(T.rows)
    # occurrences of each rational prime (odd exponent) in rows
    occ = {p: 0 for p in T.ps}
    for r in T.rows:
        for p in T.ps:
            if (r >> T.col[("r", p)]) & 1:
                occ[p] += 1
    dix, trials = rk.dixon_vectors(N, T.fb, cfg["n_dixon"], seed, cfg["dixon_max_trials"])
    dix_in = [RS.contains(T.rational_vector(1, e)) for e in dix]
    # single-prime vectors e_p: harness-only label  p in Q_N  (Legendre +1 mod every factor)
    sp_rec = []
    for p in T.ps:
        if N % p == 0:
            continue
        qr = all(gmpy2.legendre(p, q) == 1 for q in fs)
        sp_rec.append(dict(p=p, occ=occ[p], qr=qr,
                           inspan=RS.contains(1 << T.col[("r", p)])))
    return dict(exp="e3", cls=cls, bits=bits, seed=seed, N=N, ncols=T.ncols, nrows=len(T.rows),
                rank=RS.rank, target=target, dixon_trials=trials, dixon_in=dix_in,
                single=sp_rec, m=T.m, r=[T.r0, T.r1])


def main():
    exp, cfgp, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    cfg = json.load(open(cfgp))
    jobs = []
    for cls in cfg["classes"]:
        for bits in cfg["bits"]:
            for seed in range(cfg["seeds"]):
                if exp == "e1":
                    for k in cfg["k"]:
                        for g in cfg["g"]:
                            jobs.append((cls, bits, seed, k, g, cfg["trials"]))
                else:
                    jobs.append((cls, bits, seed, cfg["params"]))
    fn = dict(e1=e1_job, e2=e2_job, e3=e3_job)[exp]
    with Pool(workers) as pool, open(outp, "w") as out:
        for rec in pool.imap_unordered(fn, jobs):
            out.write(json.dumps(rec) + "\n")
            out.flush()


if __name__ == "__main__":
    main()
