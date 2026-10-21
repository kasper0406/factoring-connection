"""abcsunit pre-registered experiments E1-E4 (see notes/abcsunit_prereg.md).

HARNESS code: it uses the factors p, q for the diagnostic quantities |H| and the quadratic
characters.  Operational pieces (smooth enumeration, bucketing, relations, F2 algebra, square
extraction) see only N, y, X.

Usage:  .venv/bin/python -m src.abcsunit.run_experiment [E123|E4] [--smoke]
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from math import gcd
from typing import Dict, List

import numpy as np

from src.common.instances import balanced_semiprime
from src.abcsunit import core, model

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CFG = os.path.join(ROOT, "experiments", "configs", "abcsunit.json")
RAW = os.path.join(ROOT, "experiments", "raw")


def iroot_floor(a: int, b: int) -> int:
    """floor(a ** (1/b)) exactly."""
    x = int(round(a ** (1.0 / b)))
    while x ** b > a:
        x -= 1
    while (x + 1) ** b <= a:
        x += 1
    return x


def X_of(N: int, Knum: int, Kden: int) -> int:
    return iroot_floor(N ** Knum, Kden)


def y_of(nbits: int, kappa: float) -> int:
    # y = round((n ln 2)^kappa): depends only on the public bit length
    return int(round((nbits * math.log(2)) ** kappa))


def legendre_table(p: int) -> np.ndarray:
    t = np.zeros(p, dtype=np.int8)
    sq = (np.arange(1, p, dtype=np.int64) ** 2) % p
    t[:] = -1
    t[0] = 0
    t[np.unique(sq)] = 1
    return t


def char_vector(primes, p) -> int:
    bits = 0
    for j, l in enumerate(primes):
        if pow(l, (p - 1) // 2, p) == p - 1:
            bits |= 1 << j
    return bits


def f2_rank_of(vecs) -> int:
    r, _ = core.f2_dependencies([v for v in vecs if v], 0)
    return r


def run_cell(n: int, seed: int, kappa: float, Kfrac, cfg: Dict) -> Dict:
    t0 = time.time()
    inst = balanced_semiprime(n, seed)
    N = inst.N
    p, q = inst.factors
    y = y_of(n, kappa)
    P = core.primes_upto(y)
    Knum, Kden = Kfrac
    X = X_of(N, Knum, Kden)
    rec = {"n": n, "seed": seed, "kappa": kappa, "K": f"{Knum}/{Kden}", "N": N, "y": y, "pi_y": len(P),
           "X": X}
    try:
        vals = core.smooth_numbers(X, P, cap=cfg["psi_cap"])
    except core.TooMany:
        rec["skipped"] = "psi_cap"
        return rec
    Psi = int(vals.size)
    rec["Psi"] = Psi
    # ---------------- E1: congruent pair counts ----------------
    pairs, occupied, cmax = core.congruent_pair_count(vals, N)
    csz = core.residue_class_sizes(vals, N).astype(object)
    sumsq = int(sum(int(c) * int(c) for c in csz))
    H = core.subgroup_order_diag(P, p, q)  # DIAGNOSTIC
    phi = (p - 1) * (q - 1)
    # registered models (floating point: model evaluation only)
    vals_ext = vals  # density windows use only smooth numbers <= X (left window near X)
    M0 = model.pair_model_fast(vals, N, X, H, P, vals_ext, with_series=False)
    M1 = model.pair_model_fast(vals, N, X, H, P, vals_ext, with_series=True)
    M1k = model.pair_model_fast(vals, N, X, H, P, vals_ext, with_series=True,
                                k_filter=lambda k: core.is_smooth(k, P))
    del vals_ext
    Pnaive = Psi * (Psi - 1) / (2 * H)
    # quadratic characters (DIAGNOSTIC)
    lp, lq = legendre_table(p), legendre_table(q)
    chp = lp[(vals % p)].astype(np.int64)
    chq = lq[(vals % q)].astype(np.int64)
    Tp, Tq, Tpq = int(chp.sum()), int(chq.sum()), int((chp * chq).sum())
    cp, cq = char_vector(P, p), char_vector(P, q)
    quad = []
    for name, T, vec in (("p", Tp, cp), ("q", Tq, cq), ("pq", Tpq, cp ^ cq)):
        if vec != 0:  # character nontrivial on H
            quad.append((name, T))
    corr = 1.0 + sum((T / Psi) ** 2 for _, T in quad)
    rec.update({"pairs": pairs, "occupied": occupied, "cmax": cmax, "sumsq": sumsq, "H": H, "phi": phi,
                "H_eq_phi": H == phi, "propAB_holds": sumsq * H >= Psi * Psi,
                "M0": M0, "M1": M1, "M1_smoothk": M1k, "Pnaive": Pnaive,
                "R0": pairs / M0 if M0 > 0 else None, "R1": pairs / M1 if M1 > 0 else None,
                "T_quad": {nm: T for nm, T in quad}, "quad_corr": corr,
                "R1b": pairs / (M1 * corr) if M1 > 0 else None})
    # ---------------- relations ----------------
    buckets = core.congruent_buckets(vals, N)
    del chp, chq
    # all pairs: coprime counts, triply smooth, abc quality
    kmax = X // N + 1
    ksm = core.smooth_numbers(kmax, P)
    ksmooth = np.zeros(kmax + 1, dtype=bool)
    ksmooth[ksm] = True
    radk = np.ones(kmax + 1, dtype=np.int64)  # rad over all primes (sieve), exact
    for l in core.primes_upto(kmax):
        radk[l::l] *= l
    ncop = 0
    ntrip = 0
    ntrip_all = 0
    nq1 = 0
    nq11 = 0
    qmax = 0.0
    abcv_viol = 0
    trip_examples: List = []
    epsq = cfg["abc_eps"]
    Pset = P
    radcache: Dict[int, int] = {}

    def rad_s(x: int) -> int:
        r = radcache.get(x)
        if r is None:
            r = core.rad_smooth(x, Pset)
            if len(radcache) < 2_000_000:
                radcache[x] = r
        return r

    forest: List = []  # primitive relations forming a spanning forest of each residue class
    for b in buckets:
        L = len(b)
        parent = list(range(L))

        def find(t):
            while parent[t] != t:
                parent[t] = parent[parent[t]]
                t = parent[t]
            return t

        for j in range(1, L):
            v = b[j]
            for i in range(j):
                u = b[i]
                k = (v - u) // N
                if gcd(u, v) == 1:
                    ri, rj = find(i), find(j)
                    if ri != rj:
                        parent[ri] = rj
                        forest.append((u, v, k))
                if ksmooth[k]:
                    ntrip_all += 1
                if gcd(u, v) != 1:
                    continue
                ncop += 1
                if ksmooth[k]:
                    ntrip += 1
                    if len(trip_examples) < 20:
                        trip_examples.append((u, v, k))
                # abc quality of (u, -v, kN): exact radical
                g = gcd(int(radk[k]), N)
                R = rad_s(u) * rad_s(v) * int(radk[k]) * (N // g)
                qv = math.log(v) / math.log(R)
                qmax = max(qmax, qv)
                if v > R:
                    nq1 += 1
                if qv > 1 + epsq:
                    nq11 += 1
                    # Prop ABCV necessary condition (with K_eps = 1): (i) k*g/rad(k) > N^eps/2 or (ii) u/rad(u) > N/2
                    ci = (k * g / int(radk[k])) > (N ** epsq) / 2
                    cii = (u / rad_s(u)) > N / 2
                    if not (ci or cii):
                        abcv_viol += 1
    rec.update({"coprime_pairs": ncop, "triply": ntrip, "triply_all": ntrip_all,
                "R3": ntrip_all / M1k if M1k > 0 else None,
                "q_gt1": nq1, "q_gt1eps": nq11, "qmax": qmax, "abcv_violations": abcv_viol,
                "evertse_bound_log10": math.log10(3) + (2 * (len(P) + 2) + 3) * math.log10(7),
                "triply_examples": trip_examples})
    # ---------------- E2: rank and fruitfulness on primitive spanning-forest relations ----------------
    chain = sorted(forest, key=lambda t: (t[1], t[0]))
    rec["forest_size"] = len(chain)
    R = 3 * len(P) + 30
    chain = chain[:R]
    rows, exps = [], []
    for u, v, k in chain:
        bts, e = core.relation_row(u, v, P)
        rows.append(bts)
        exps.append(e)
    rank, deps = core.f2_dependencies(rows, len(P))
    span = f2_rank_of([cp, cq])
    # operational: seeded random combinations of the dependency basis
    rng = random.Random(f"abcsunit|E2|{n}|{seed}|{kappa}|{Knum}/{Kden}")
    gs = []
    if deps:
        for _ in range(cfg["max_deps"]):
            m = 0
            while m == 0:
                for d in deps:
                    if rng.random() < 0.5:
                        m ^= d
            gs.append(core.square_from_dependency(m, chain, exps, P, N))
    fruit = sum(1 for g in gs if 1 < g < N)
    # DIAGNOSTIC: dimension of the image of the z-map (dependency -> square root of 1) on the basis
    zbits = []
    for m in deps:
        tot = [0] * len(P)
        V = 1
        for i in range(len(chain)):
            if m >> i & 1:
                tot = [a + b for a, b in zip(tot, exps[i])]
                V = V * chain[i][1] % N
        W = 1
        for l, t in zip(P, tot):
            W = W * pow(l, t // 2, N) % N
        z = W * pow(V, -1, N) % N
        zbits.append((0 if z % p == 1 else 1) | ((0 if z % q == 1 else 2)))
    zdim = f2_rank_of(zbits)
    rec.update({"relations_used": len(chain), "R_target": R, "rank": rank, "char_span": span,
                "rank_pred": len(P) - span, "n_deps": len(deps), "deps_tested": len(gs), "fruitful": fruit,
                "z_image_dim": zdim, "factor_found": any(1 < g < N for g in gs),
                "factor_class": "STRUCTURAL" if fruit else None})
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


def run_E4(cfg: Dict, out) -> None:
    for nprim in cfg["E4"]["nprimes"]:
        P = core.primes_upto(60)[:nprim]
        for hb in cfg["E4"]["Hbits"]:
            H = 1 << hb
            vals = core.smooth_numbers(H, P)
            if vals.size > cfg["E4"]["psi_cap"]:
                out.write(json.dumps({"exp": "E4", "nprimes": nprim, "Hbits": hb, "skipped": "psi_cap"}) + "\n")
                continue
            t0 = time.time()
            tr = core.primitive_smooth_triples(vals, H)
            rec = {"exp": "E4", "nprimes": nprim, "Hbits": hb, "Psi": int(vals.size), "count": len(tr),
                   "vkm": cfg["E4"]["vkm"][nprim - 1], "max_c": max(t[2] for t in tr) if tr else None,
                   "evertse_bound": 3 * 7 ** (2 * nprim + 3), "seconds": round(time.time() - t0, 1)}
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print(rec, flush=True)


def main():
    with open(CFG) as f:
        cfg = json.load(f)
    which = sys.argv[1] if len(sys.argv) > 1 else "E123"
    smoke = "--smoke" in sys.argv
    os.makedirs(RAW, exist_ok=True)
    if which == "E123":
        path = os.path.join(RAW, "abcsunit_E123_smoke.jsonl" if smoke else "abcsunit_E123.jsonl")
        cells = cfg["E123"]["cells"]
        seeds = cfg["E123"]["seeds"]
        if smoke:
            cells = cfg["E123"]["smoke_cells"]
            seeds = seeds[:1]
        with open(path, "w") as out:
            for c in cells:
                for s in seeds:
                    rec = run_cell(c["n"], s, c["kappa"], c["K"], cfg["E123"])
                    rec["exp"] = "E123"
                    out.write(json.dumps(rec) + "\n")
                    out.flush()
                    print({k: rec.get(k) for k in ("n", "seed", "kappa", "K", "Psi", "pairs", "R1", "R1b",
                                                    "H_eq_phi", "triply", "R3", "rank", "rank_pred",
                                                    "fruitful", "deps_tested", "q_gt1eps", "seconds")},
                          flush=True)
    elif which == "E4":
        path = os.path.join(RAW, "abcsunit_E4.jsonl")
        with open(path, "w") as out:
            run_E4(cfg, out)


if __name__ == "__main__":
    main()
