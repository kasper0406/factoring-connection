"""e6: verification of Theorem R (DCL-wc succeeds on every covered pair).  HARNESS.
Pre-registered in notes/dchworst_prereg.md, addendum A (frozen before these runs).

  e6a  exhaustive: every pair of S_n, n = 26, 28; H = 3; Y = max(3, 2^floor(n/6) // 2);
       4 blocks of primes = 3 (mod 4) in (2^k Y, 2^(k+1) Y].
       designs: "pos" (D' weights V(H), no (1,1)) and "sgn" (signed + (1,1)).
       modes: "full" (Lemma O, T = 2^ceil(n/3), B0 = 30) and "force2" (base 2, no Lemma O).
       D' (dcl_det) is run on the "pos" design for comparison.
  e6b  adversarial families at n = 40, 48, 56 (q | 2^k - 1; p, q | 2^k - 1, 2^k' - 1;
       equal-order pairs p = 2D + 1, q = 3D + 1 with ord_p(2) = ord_q(2) = D = l0 * o).

Run:  PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments_wc e6a|e6b
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from math import isqrt
from multiprocessing import Pool

import gmpy2
import numpy as np

from src.common.metrics import Metrics
from src.dchworst import coverage as cv
from src.dchworst import wcalg as wc
from src.dchworst.experiments import blocks
from src.dcltheorem.algorithm import dcl_det
from src.dcltheorem.params import Design

RAW = "experiments/raw"


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return os.path.join(RAW, f"dchworst_{name}.jsonl")


def covered_mask(P, Q, bl, H, design):
    """Harness coverage: some block, some band prime l and weight with l | Bp - Aq."""
    cov = np.zeros(len(P), dtype=bool)
    for b in bl:
        if not b:
            continue
        h = cv.sweep(P, Q, b, H, design == "sgn", record_every=[])["h"]
        c = h <= H
        if design == "sgn":                       # the diagonal weight (1, 1)
            for l in b:
                c |= ((Q - P) % l) == 0
        cov |= c
    return cov


def _job(args):
    n, H, Y, design, mode, pairs = args
    bl = blocks(Y)
    T = 1 << (-(-n // 3))
    out = []
    for p, q in pairs:
        N = p * q
        met = Metrics()
        r = wc.dcl_wc(N, bl, H, T, 30, signed=(design == "sgn"), diag=(design == "sgn"),
                      met=met, force_base=(2 if mode == "force2" else None))
        f = r["factor"]
        rec = {"p": p, "q": q, "ok": f is not None, "proper": f is None or (1 < f < N and N % f == 0),
               "status": r["status"], "eq": r["log"].get("equal_order"),
               "Ngiants": r["log"].get("N_giants", 0),
               "explore_evals": met.get("explore_evals"), "explore_frontier": met.get("explore_frontier"),
               "leaf_candidates": met.get("leaf_candidates")}
        if design == "pos" and mode == "full":
            rd = dcl_det(N, Design(n=n, H=H, Y=Y, K=4))
            rec["dprime_ok"] = rd["factor"] is not None
        out.append(rec)
    return out


def e6a():
    path = _out("e6a")
    if os.path.exists(path):
        os.remove(path)
    for n in (26, 28):
        P, Q = cv.all_pairs(n)
        H = 3
        Y = max(3, (1 << (n // 6)) // 2)
        bl = blocks(Y)
        for design in ("pos", "sgn"):
            cov = covered_mask(P, Q, bl, H, design)
            pairs = list(zip(P.tolist(), Q.tolist()))
            for mode in ("full", "force2"):
                t0 = time.time()
                chunks = [pairs[i:i + 2000] for i in range(0, len(pairs), 2000)]
                with Pool(6) as pool:
                    res = pool.map(_job, [(n, H, Y, design, mode, c) for c in chunks])
                recs = [x for c in res for x in c]
                agg = {"n": n, "H": H, "Y": Y, "design": design, "mode": mode,
                       "band_sizes": [len(b) for b in bl], "npairs": len(recs),
                       "covered": int(cov.sum()),
                       "covered_fail": 0, "uncovered_ok": 0, "uncovered": int((~cov).sum()),
                       "improper": 0, "status": {}, "equal_order_runs": 0,
                       "max_frontier_over_Ngiants": 0.0, "explore_evals": 0, "Ngiants": 0,
                       "covered_fail_examples": []}
                if design == "pos" and mode == "full":
                    agg["dprime_fail_covered"] = 0
                    agg["dprime_fail_covered_wc_ok"] = 0
                for i, rc in enumerate(recs):
                    c = bool(cov[i])
                    agg["status"][rc["status"]] = agg["status"].get(rc["status"], 0) + 1
                    agg["improper"] += (not rc["proper"])
                    agg["equal_order_runs"] += rc["eq"] is not None
                    agg["explore_evals"] += rc["explore_evals"]
                    agg["Ngiants"] += rc["Ngiants"]
                    if rc["Ngiants"]:
                        agg["max_frontier_over_Ngiants"] = max(
                            agg["max_frontier_over_Ngiants"], rc["explore_frontier"] / rc["Ngiants"])
                    if c and not rc["ok"]:
                        agg["covered_fail"] += 1
                        if len(agg["covered_fail_examples"]) < 20:
                            agg["covered_fail_examples"].append([rc["p"], rc["q"], rc["status"]])
                    if not c and rc["ok"]:
                        agg["uncovered_ok"] += 1
                    if "dprime_ok" in rc and c and not rc["dprime_ok"]:
                        agg["dprime_fail_covered"] += 1
                        agg["dprime_fail_covered_wc_ok"] += rc["ok"]
                agg["seconds"] = round(time.time() - t0, 1)
                with open(path, "a") as f:
                    f.write(json.dumps(agg) + "\n")
                print(json.dumps({k: v for k, v in agg.items() if k != "covered_fail_examples"}),
                      flush=True)


# ---------------------------------------------------------------------------
def _fam_pairs(n, rng):
    """Adversarial families (harness: sympy factorisation of 2^k - 1 is diagnostic)."""
    import sympy
    lo, hi = 1 << (n - 1), 1 << n
    plo, phi = isqrt(lo // 2), isqrt(hi) + 1
    qhi = isqrt(2 * hi) + 1
    cand = set()
    for k in range(2, 241):
        for f in sympy.factorint((1 << k) - 1, limit=10 ** 6).keys():
            if plo // 2 <= f <= qhi and sympy.isprime(f):
                cand.add(int(f))
    cand = sorted(cand)
    fam = {"small2_q": [], "small2_both": [], "eqord": []}
    for q in cand:
        for _ in range(4):
            p = int(gmpy2.next_prime(rng.randint(q // 2, q)))
            a, b = sorted((p, q))
            if b <= 2 * a and lo <= a * b < hi:
                fam["small2_q"].append((a, b))
    for i, a in enumerate(cand):
        for b in cand[i + 1:]:
            if b <= 2 * a and lo <= a * b < hi:
                fam["small2_both"].append((a, b))
    # equal order: p = j D + 1, q = k D + 1 with ord_p(2) = ord_q(2) = D = l0 * o, l0 a band
    # prime and o small (so gamma = 2^L has the same small order o mod p and q)
    import sympy
    Y = max(3, (1 << (n // 6)) // 2)
    band_primes = [l for b in blocks(Y) for l in b]

    def exact(D, x):
        return pow(2, D, x) == 1 and all(pow(2, D // r, x) != 1 for r in sympy.factorint(D))
    tries = 0
    while len(fam["eqord"]) < 40 and tries < 3 * 10 ** 6:
        tries += 1
        l0 = rng.choice(band_primes)
        o = rng.randint(1, 60)
        D = l0 * o
        j = 2 * rng.randint(max(1, plo // (2 * D)), max(1, phi // (2 * D)))
        p = j * D + 1
        if not (plo <= p <= phi) or pow(2, D, p) != 1 or not gmpy2.is_prime(p) or not exact(D, p):
            continue
        for k in range(j + 2, 2 * j + 1, 2):
            q = k * D + 1
            if q <= 2 * p and lo <= p * q < hi and pow(2, D, q) == 1 and gmpy2.is_prime(q) \
                    and exact(D, q):
                fam["eqord"].append((p, q))
                break
    return fam


def e6b():
    path = _out("e6b")
    if os.path.exists(path):
        os.remove(path)
    for n in (40, 48, 56):
        rng = random.Random(f"dchworst-e6b|{n}")
        fam = _fam_pairs(n, rng)
        H = 3
        Y = max(3, (1 << (n // 6)) // 2)
        bl = blocks(Y)
        for name, pairs in fam.items():
            if not pairs:
                continue
            P = np.array([a for a, _ in pairs], dtype=np.int64)
            Q = np.array([b for _, b in pairs], dtype=np.int64)
            for design in ("pos", "sgn"):
                cov = covered_mask(P, Q, bl, H, design)
                for mode in ("full", "force2"):
                    recs = _job((n, H, Y, design, mode, pairs))
                    agg = {"n": n, "family": name, "design": design, "mode": mode,
                           "npairs": len(recs), "covered": int(cov.sum()),
                           "covered_fail": sum(1 for i, r in enumerate(recs) if cov[i] and not r["ok"]),
                           "uncovered_ok": sum(1 for i, r in enumerate(recs) if not cov[i] and r["ok"]),
                           "improper": sum(1 for r in recs if not r["proper"]),
                           "equal_order_runs": sum(1 for r in recs if r["eq"] is not None),
                           "status": {}}
                    for r in recs:
                        agg["status"][r["status"]] = agg["status"].get(r["status"], 0) + 1
                    if design == "pos" and mode == "full":
                        agg["dprime_fail_covered"] = sum(
                            1 for i, r in enumerate(recs) if cov[i] and not r["dprime_ok"])
                    with open(path, "a") as f:
                        f.write(json.dumps(agg) + "\n")
                    print(json.dumps(agg), flush=True)


if __name__ == "__main__":
    {"e6a": e6a, "e6b": e6b}[sys.argv[1]]()
