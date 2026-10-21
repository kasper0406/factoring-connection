"""E2 (pre-registered, notes/schnorrcvp_prereg.md): relation yield of exact / list / approximate CVP outputs
in Schnorr-type prime-number lattices, against the residual-smoothness heuristics H_res and H_res'.

Usage: .venv/bin/python -m src.schnorrcvp.run_e2 [--quick]
Writes experiments/raw/schnorrcvp_E2.jsonl (one line per lattice) and schnorrcvp_E2_factor.jsonl.
"""
from __future__ import annotations

import json
import math
import sys
import time

import mpmath as mp

from src.common.instances import balanced_semiprime

from .core import first_primes, lambda_vector, squares_dependency_factor
from .lattice import analyse_vector, babai, build_lattice, list_cvp, reduce_basis
from .psi import PsiBracket

BITS = [20, 24, 28, 32, 36, 40]
NS = [12, 20, 30, 40]
VARIANTS = [  # (label, weight, c_exp or None for Schnorr's N' = N^{1/(n+1)})
    ("S21perm", "perm", None),
    ("S21sqrt", "sqrtperm", None),
    ("A0.5", "sqrtlnpperm", 0.5),
    ("A1.0", "sqrtlnpperm", 1.0),
    ("A1.25", "sqrtlnpperm", 1.25),
]
NSEEDS_N = 3
NLAT = 10
K_LIST = 20


class Res:
    """Residual-smoothness predictions from one bracket table (model numerics, float)."""

    def __init__(self, P, lnmax):
        self.P = P
        self.pb = PsiBracket(P, 200, lnmax)
        self.lnmax = lnmax
        self.cache = {}

    def d(self, lnx, h=0.1):
        key = round(lnx * 100)
        if key not in self.cache:
            x = max(mp.mpf(key) / 100 - h / 2, mp.mpf(0))
            self.cache[key] = float(self.pb.density(x, h))
        return self.cache[key]

    def naive(self, w):
        return self.d(math.log(abs(w)))

    def refined(self, w, Q):
        """H_res': |w| is coprime to Q = supp(uv); CRT-uniform elsewhere (independence-product form)."""
        lx = math.log(abs(w))
        base = self.d(lx)
        if base == 0.0:
            return 0.0
        f = 1.0
        for q in Q:
            frac_div = self.d(lx - math.log(q)) / (q * base)  # fraction of smooth m near x divisible by q
            f *= max(0.0, 1.0 - frac_div) / (1.0 - 1.0 / q)
        return base * f


def supp(x, P):
    return [p for p in P if x % p == 0]


def main(quick: bool = False):
    bits_l = [20, 28] if quick else BITS
    ns = [12, 20] if quick else NS
    nseeds = 1 if quick else NSEEDS_N
    nlat = 2 if quick else NLAT
    out = open("experiments/raw/schnorrcvp_E2.jsonl", "w")
    fout = open("experiments/raw/schnorrcvp_E2_factor.jsonl", "w")
    for bits in bits_l:
        for s in range(nseeds):
            N = balanced_semiprime(bits, s).N
            pool = {}
            for n in ns:
                P = first_primes(n)
                lnN = math.log(N)
                res = Res(P, 6 * lnN + 40)
                for (lab, wt, c) in VARIANTS:
                    for ls in range(nlat):
                        t0 = time.time()
                        B, t, d, _, meta = build_lattice(N, n, c_exp=c, weight=wt, schnorr_scale=(c is None), seed=ls)
                        reduce_basis(B)
                        bab_lll = babai(B, t)
                        L = list_cvp(B, t, K_LIST)
                        Bk = B.__copy__()
                        reduce_basis(Bk, beta=10)
                        bab_bkz = babai(Bk, t)
                        recs = []
                        for kind, vec in [("babai_lll", bab_lll), ("babai_bkz10", bab_bkz)] + \
                                [("list%d" % i, v) for i, v in enumerate(L)]:
                            a = analyse_vector(vec, d, P, N, t)
                            u, v, w = a["u"], a["v"], a["w"]
                            if w == 0:
                                pn = pr = None
                            elif math.log(abs(w)) > res.lnmax - 1:
                                pn = pr = None
                            else:
                                pn = res.naive(w)
                                pr = res.refined(w, supp(u * v, P))
                            big = max(u.bit_length(), v.bit_length(), abs(w).bit_length()) > 4000
                            recs.append({"kind": kind, "u": None if big else str(u), "v": None if big else str(v),
                                         "w": None if big else str(w), "huge": big, "rel": a["rel"],
                                         "dist2": None if big else str(a["dist2"]), "nnz": a["nnz"], "p_naive": pn, "p_ref": pr,
                                         "log2w": (w.bit_length() if big else (math.log2(abs(w)) if w else None)),
                                         "log2v": (v.bit_length() if big else math.log2(v))})
                            if a["rel"]:
                                lv = lambda_vector(u, w, first_primes(max(ns)))
                                pool[(u, w)] = lv
                        out.write(json.dumps({"bits": bits, "seed": s, "N": str(N), "n": n, "variant": lab,
                                              "lat_seed": ls, "outputs": recs,
                                              "t": round(time.time() - t0, 3)}) + "\n")
                        out.flush()
                print(bits, s, n, "pool", len(pool), flush=True)
            # end-to-end correctness check (toy sizes): congruence of squares from pooled relations
            Pmax = first_primes(max(ns))
            vecs = [v for v in pool.values() if v is not None]
            g, ctr = squares_dependency_factor(vecs, Pmax, N) if len(vecs) > len(Pmax) + 1 else (None, {})
            fout.write(json.dumps({"bits": bits, "seed": s, "N": str(N), "distinct_relations": len(vecs),
                                   "fb_size": len(Pmax) + 1, "factor": (str(g) if g else None),
                                   "divides": (N % g == 0 if g else None), "ctr": ctr}) + "\n")
            fout.flush()


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
