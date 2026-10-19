"""HARNESS for rkrank experiment E1 (pre-registered in notes/rkrank_prereg.md).

DIAGNOSTIC: uses the factors of N only to build the Gen functionals (Λ/2)^⊥ and to label results.
Usage: .venv/bin/python -m src.rkrank.harness CONFIG.json OUT.jsonl [workers]
"""
from __future__ import annotations

import json
import random
import sys
import time
from math import gcd
from multiprocessing import Pool
from typing import Dict, List

from ..nfsfruit.instances import composite
from . import gf2
from . import sampler as smp
from .toynfs import ToyField, gen_functionals


def peel_core(rows: List[int], ideal_mask: int) -> List[int]:
    """Indices of rows surviving iterated removal of rows that contain an ideal column
    (bit in ideal_mask) occurring (oddly) in exactly one remaining row (algebraic 2-core)."""
    alive = set(range(len(rows)))
    while True:
        cnt: Dict[int, int] = {}
        who: Dict[int, int] = {}
        for i in alive:
            v = rows[i] & ideal_mask
            while v:
                lb = v & -v
                c = lb.bit_length() - 1
                cnt[c] = cnt.get(c, 0) + 1
                who[c] = i
                v ^= lb
        kill = {who[c] for c, k in cnt.items() if k == 1}
        if not kill:
            return sorted(alive)
        alive -= kill


def level_record(T: ToyField, rows: List[int], gen: List[int], vallee, structured, rng) -> Dict:
    D = T.D_rat
    E = gf2.Echelon(rows)
    Wb = E.basis_below(D)
    high = ((1 << T.ncols) - 1) ^ ((1 << D) - 1)
    rankA = gf2.rank(r & high for r in rows)
    dimW = len(Wb)
    assert dimW == E.rank - rankA
    gen_ok = gf2.all_annihilate(gen, Wb)
    sign_ok = all(b & 1 == 0 for b in Wb)
    forced = gen + [1]
    forced_rank = gf2.rank(forced)
    in_perp_forced = gf2.annihilated_dim(forced, Wb) if gf2.rank(forced) == len(forced) else None
    dim_perp = D - dimW
    resid = dim_perp - in_perp_forced if in_perp_forced is not None else None
    occ = {p: 0 for p in T.ps}
    for r in rows:
        v = r & ((1 << D) - 2)
        while v:
            lb = v & -v
            c = lb.bit_length() - 1
            occ[T.ps[c - 1]] += 1
            v ^= lb
    n_occ = {k: sum(1 for p in T.ps if occ[p] == k) for k in range(4)}
    ck = {}
    for k in (1, 2, 3, 5):
        Pk = [T.col[("r", p)] for p in T.ps if occ[p] >= k]
        mask = gf2.mask_of(Pk)
        lam = len(Pk) - gf2.restrict_rank(gen, mask)
        ck[k] = lam - gf2.dim_in_coords(Wb, mask)
    # peeling predictor: rational primes occurring in no row of the algebraic 2-core
    ideal_mask = 0
    for key, c in T.col.items():
        if key[0] == "a":
            ideal_mask |= 1 << c
    core = peel_core(rows, ideal_mask)
    core_rat = 0
    for i in core:
        core_rat |= rows[i]
    n_peel = sum(1 for p in T.ps if not (core_rat >> T.col[("r", p)]) & 1)
    const_real = []
    for j in range(len(T.real)):
        c = T.col[("rs", j)]
        vals = {(r >> c) & 1 for r in rows}
        const_real.append(vals.pop() if len(vals) == 1 else None)
    rec = dict(const_real=const_real, M=len(rows), rankR=E.rank, rankA=rankA, dimW=dimW, dim_perp=dim_perp,
               forced_rank=forced_rank, gen_in_perp=gen_ok, sign_in_perp=sign_ok,
               resid=resid, n_occ=n_occ, ck=ck, n_peel=n_peel, core_rows=len(core))
    if vallee is not None:
        P3 = {T.col[("r", p)] for p in T.ps if occ[p] >= 3}
        vres = []
        WE = gf2.Echelon(Wb)
        for (z, v) in vallee:
            supp = {c for c in range(1, D) if (v >> c) & 1}
            vres.append(dict(neg=z < 0, inP3=supp <= P3, inW=WE.contains(v)))
        rec["vallee"] = vres
    if structured:
        def delta(phi):
            return sum(gf2.dot(r, phi) for r in rows) / len(rows)
        rat_primes = ((1 << D) - 2)
        liou = rat_primes | ideal_mask
        dense = [rng.getrandbits(T.ncols) for _ in range(10)]
        sparse = [1 << T.col[("r", p)] for p in rng.sample(T.ps, 10)]
        rec["delta"] = dict(liouville_joint=delta(liou), liouville_rat=delta(rat_primes),
                            liouville_alg=delta(ideal_mask),
                            dense=[delta(x) for x in dense], sparse=[delta(x) for x in sparse])
    return rec


def job(args) -> Dict:
    cls, bits, seed, cfg = args
    t0 = time.time()
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
    # small squares (linearised cells), keep B-smooth ones; harness records gcd accidents
    h = smp.h_two_thirds(N)
    vallee, accidental = [], 0
    for (x, z) in smp.small_squares(N, p["cells"], seed, h):
        if gcd(x, N) != 1 or gcd(z, N) != 1:
            accidental += 1
            continue
        v = T.rational_vector(z)
        if v is not None:
            vallee.append((z, v))
    levels = []
    for L in cfg["levels"]:
        M = int(L * T.ncols)
        if M > len(rows):
            break
        top = L == max(cfg["levels"])
        levels.append(dict(level=L, **level_record(T, rows[:M], gen, vallee,
                                                   structured=top, rng=rng)))
    return dict(exp="rkrank_e1", cls=cls, bits=bits, seed=seed, N=N, omega=len(fs),
                ncols=T.ncols, D_rat=T.D_rat, D_alg=T.D_alg, n_real=len(T.real),
                nchar=len(T.chars), m=T.m, r=[T.r0, T.r1], nrel=len(rows), b_last=blast,
                n_vallee=len(vallee), accidental_gcd=accidental, levels=levels,
                secs=time.time() - t0)


def main():
    cfgp, outp = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    cfg = json.load(open(cfgp))
    jobs = [(c, b, s, cfg) for c in cfg["classes"] for b in cfg["bits"] for s in cfg["seeds"]]
    with Pool(workers) as pool, open(outp, "w") as out:
        for rec in pool.imap_unordered(job, jobs):
            out.write(json.dumps(rec) + "\n")
            out.flush()


if __name__ == "__main__":
    main()
