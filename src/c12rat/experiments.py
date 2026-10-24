"""c12rat experiments E1, E2, E3 (notes/c12rat_prereg.md).  DIAGNOSTIC.

E1: label shapes and thin pairs in child boxes (exact; Lemmas LS, MULT, TF).
E2: Lemma BT-S (smooth numbers in an AP in a short interval), exact inequality on random instances.
E3: cell (Cauchy–Schwarz) decomposition of the child-box pair covariances, rational block, under the all-g and
    the primitive-g laws (e3.py).

Usage: python -m src.c12rat.experiments {E1|E2|E3} out.jsonl
"""
from __future__ import annotations

import json
import sys
import time
from math import gcd

import numpy as np

from src.pc3 import toy as pc3toy
from . import bts, e3, geom

PARAMS = dict(
    seed=20261024,
    # E1
    E1_qs_shape=(10007, 100003), E1_sigmas=(1 / 256, 1 / 64, 1 / 16, 1 / 4),
    E1_q_thin=10007, E1_n_thin=1000, E1_S=10, E1_Ks=(2, 8, 32),
    # E2
    E2_X=50_000_000, E2_y=50, E2_n=2000,
    # E3
    E3_X=20_000_000, E3_y=150, E3_M=1200, E3_S=8, E3_h=0.5, E3_Kthin=16.0, E3_qlo=500,
    E3_n_reg=12, E3_n_skew=4, E3_sigma_skew=1 / 16,
)


def run_E1(out):
    P = PARAMS
    rng = np.random.default_rng(P["seed"])
    for q in P["E1_qs_shape"]:
        sig = np.array([geom.shape(q, r)[0] for r in range(q)])
        rec = dict(kind="shape", q=q)
        for s in P["E1_sigmas"]:
            rec["F_%g" % s] = float((sig <= s).mean())
            rec["bound_%g" % s] = geom.shape_cdf_bound(s, q)
        out.write(json.dumps(rec) + "\n")
        out.flush()
    q = P["E1_q_thin"]
    S = P["E1_S"]
    for _ in range(P["E1_n_thin"]):
        r = int(rng.integers(0, q))
        sigma, t2 = geom.shape(q, r)
        ij, ab = geom.box_arrays(q, r, S)
        iu, ju, Delta, det, kappa = geom.pair_stats(ij, ab)
        assert np.all(det == q * Delta)
        n = len(ij)
        tf = geom.thin_fraction(kappa, P["E1_Ks"])
        rec = dict(kind="thin", q=q, r=r, sigma=sigma, t2=t2, n=n, mult_max=geom.mult_max(Delta),
                   mult_bound=5 * S * S + S)
        for K in P["E1_Ks"]:
            rec["f_%g" % K] = tf[float(K)]
            rec["fb_%g" % K] = geom.thin_bound(S, n, t2, K)
        out.write(json.dumps(rec) + "\n")
    out.flush()


def run_E2(out):
    P = PARAMS
    rng = np.random.default_rng(P["seed"] + 2)
    X, y = P["E2_X"], P["E2_y"]
    smooth, _ = pc3toy.build_tables(X, y, y + 1)
    cs = None
    for _ in range(P["E2_n"]):
        # DEVIATION (disclosed in notes/c12rat_results.md §10.2): the frozen code drew infeasible Z > X - 10^6 and
        # crashed after 5 instances; (D, Z/D) are now redrawn until the instance fits in the table.
        while True:
            D = int(rng.integers(10, 10_001))
            Zterms = int(np.exp(rng.uniform(np.log(100), np.log(10_000))))
            Z = Zterms * D
            if Z < X - 1_000_002:
                break
        while True:
            a = int(rng.integers(1, D))
            if gcd(a, D) == 1:
                break
        X0 = int(rng.integers(1_000_000, X - Z - 1))
        W = max(1, Z // D)
        if W > X0:
            continue
        L = bts.lhs(smooth, X0, Z, a, D)
        R = bts.rhs(smooth, X0, Z, D, y, W)
        RA = bts.rhs_analytic(smooth, Z, D, y, W)
        dens = float(smooth[X0 + 1:X0 + Z + 1].mean())
        phiD = D * np.prod([1 - 1 / p for p in pc3toy.factor_small(D)])
        expect = Z / phiD * dens
        out.write(json.dumps(dict(D=D, a=a, X0=X0, Z=Z, W=W, lhs=L, rhs=R, rhs_an=RA, expect=expect)) + "\n")
    out.flush()


def run_E3(out):
    P = PARAMS
    rng = np.random.default_rng(P["seed"] + 3)
    smooth, _ = pc3toy.build_tables(P["E3_X"], P["E3_y"], P["E3_y"] + 1)
    ps = pc3toy.primes_upto(2 * P["E3_qlo"])
    ps = ps[ps > P["E3_qlo"]]
    labs = []
    for _ in range(P["E3_n_reg"]):
        q = int(rng.choice(ps))
        labs.append((q, int(rng.integers(0, q)), "uniform"))
    k = 0
    while k < P["E3_n_skew"]:
        q = int(rng.choice(ps))
        r = int(rng.integers(0, q))
        if geom.shape(q, r)[0] < P["E3_sigma_skew"]:
            labs.append((q, r, "skew"))
            k += 1
    g_all = e3.all_g(P["E3_M"])
    g_prim = e3.all_primitive_g(P["E3_M"])
    for (q, r, kind) in labs:
        sigma, t2 = geom.shape(q, r)
        for law, g in (("all", g_all), ("prim", g_prim)):
            t0 = time.time()
            rec = e3.label_run(q, r, P["E3_S"], g, smooth, P["E3_y"], P["E3_h"], P["E3_Kthin"],
                               prim_model=(law == "prim"))
            rec = {kk: (float(v) if isinstance(v, (np.floating,)) else v) for kk, v in rec.items()}
            rec.update(kind=kind, law=law, sigma=sigma, t2=t2, secs=round(time.time() - t0, 1))
            out.write(json.dumps(rec) + "\n")
            out.flush()


if __name__ == "__main__":
    which, path = sys.argv[1], sys.argv[2]
    with open(path, "w") as out:
        {"E1": run_E1, "E2": run_E2, "E3": run_E3}[which](out)
