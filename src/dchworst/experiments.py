"""Pre-registered experiments for dchworst (notes/dchworst_prereg.md).  HARNESS code.

  e1  exhaustive worst case of h_J over ALL pairs of S_n (n = 32, 36), nested
      prefixes of the band (Y, 2Y], Y = 2^9, 2^11, 2^13, signed and positive weights.
  e2  sampled tail at the D'-scaling Y = 2^(n/6), four blocks (Y, 16Y], n = 60, 72, 84.
  e3  structured (adversarial) families at n = 72, block 0, vs a random baseline.
  e4  exhaustive D'-design coverage (positive weights, 4 blocks, Y = 2^(n/6) H^(-1/3)),
      n = 24, 28, 32, H = 2..10: uncovered and uncovered-in-every-clean-block counts.
  e5  random search for pairs uncovered in all four D' blocks, n = 64 (H = 7, 9).

Run:  PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments e1|e2|e3|e4|e5|summarize
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from math import isqrt
from multiprocessing import Pool
from typing import Dict, List

import gmpy2
import numpy as np

from src.dchworst import coverage as cv
from src.dcltheorem.params import icbrt, pow2_frac, primes_upto

RAW = "experiments/raw"
CFG = "experiments/configs"
ANA = "experiments/analysis"


def _out(name: str) -> str:
    os.makedirs(RAW, exist_ok=True)
    return os.path.join(RAW, f"dchworst_{name}.jsonl")


def _write(path: str, rec: Dict) -> None:
    with open(path, "a") as f:
        f.write(json.dumps(rec) + "\n")


def record_points(J: int) -> List[int]:
    pts = {1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, J}
    return sorted(x for x in pts if x <= J)


def h_pred(npairs: int, surv: List[float]) -> int:
    """Smallest H with npairs * Pr_Poisson[h > H] < 1 (len(surv) = Hcap + 1)."""
    for H, s in enumerate(surv):
        if npairs * s < 1:
            return H
    return len(surv)


# ---------------------------------------------------------------------------
def e1():
    path = _out("e1")
    if os.path.exists(path):
        os.remove(path)
    for n in (32, 36):
        P, Q = cv.all_pairs(n)
        for logY in (9, 11, 13):
            Y = 1 << logY
            pr = cv.band(Y)
            Hcap = isqrt(2 * Y) + 1
            for signed in (True, False):
                t0 = time.time()
                out = cv.sweep(P, Q, pr, Hcap, signed, record_every=record_points(len(pr)))
                for rec in out["prefix"]:
                    i = rec.pop("argmax")
                    surv = rec["poisson_surv"]
                    hist = rec["hist"]
                    tail_obs = [int(sum(hist[H + 1:])) for H in range(Hcap + 1)]
                    tail_exp = [len(P) * s for s in surv]
                    _write(path, {
                        "n": n, "Y": Y, "signed": signed, "Hcap": Hcap, "npairs": int(len(P)),
                        "J": rec["J"], "l_max": rec["l_max"], "hstar": rec["hstar"],
                        "worst_pair": [int(P[i]), int(Q[i])],
                        "hpred": h_pred(len(P), surv),
                        "tail_obs": tail_obs, "tail_exp": tail_exp,
                    })
                print(f"e1 n={n} Y=2^{logY} signed={signed} J={len(pr)} "
                      f"h*={out['prefix'][-1]['hstar']} ({time.time()-t0:.0f}s)", flush=True)


# ---------------------------------------------------------------------------
# samplers
# ---------------------------------------------------------------------------
def _sample_chunk(args):
    n, seed, count = args
    rng = random.Random(f"dchworst-pair|n={n}|{seed}")
    plo = isqrt(1 << (n - 2))
    phi = isqrt((1 << n) - 1) + 1
    qlo = isqrt(1 << (n - 1))
    qhi = isqrt(1 << (n + 1)) + 1
    lo, hi = 1 << (n - 1), 1 << n
    out = []
    while len(out) < count:
        p = rng.randint(plo, phi)
        q = rng.randint(qlo, qhi)
        if not (p < q <= 2 * p) or not (lo <= p * q < hi):
            continue
        if gmpy2.is_prime(p) and gmpy2.is_prime(q):
            out.append((p, q))
    return out


def sample_pairs(n: int, count: int, seed0: int, procs: int = 30):
    per = 2000
    tasks = [(n, seed0 * 100000 + k, min(per, count - k * per))
             for k in range((count + per - 1) // per)]
    with Pool(procs) as pool:
        chunks = pool.map(_sample_chunk, tasks)
    pairs = [pq for c in chunks for pq in c]
    P = np.array([p for p, _ in pairs], dtype=np.int64)
    Q = np.array([q for _, q in pairs], dtype=np.int64)
    return P, Q


def _sweep_job(args):
    P, Q, pr, Hcap, signed = args
    return cv.sweep(P, Q, pr, Hcap, signed, record_every=[len(pr)])


def blocks(Y: int, K: int = 4) -> List[List[int]]:
    ps = primes_upto((1 << K) * Y)
    return [[l for l in ps if (1 << k) * Y < l <= (1 << (k + 1)) * Y and l % 4 == 3]
            for k in range(K)]


# ---------------------------------------------------------------------------
def e2():
    path = _out("e2")
    if os.path.exists(path):
        os.remove(path)
    plan = {60: 200000, 72: 100000, 84: 50000}
    for n, cnt in plan.items():
        P, Q = sample_pairs(n, cnt, seed0=n)
        Y = pow2_frac(n, 1, 6)
        bl = blocks(Y)
        Hcap = 64
        for signed in (True, False):
            jobs = [(P, Q, bl[k], Hcap, signed) for k in range(4)]
            jobs.append((P, Q, [l for b in bl for l in b], Hcap, signed))
            with Pool(5) as pool:
                outs = pool.map(_sweep_job, jobs)
            for k, out in enumerate(outs):
                rec = out["prefix"][-1]
                hist = rec["hist"]
                surv = rec["poisson_surv"]
                i = int(np.argmax(out["h"]))
                _write(path, {
                    "n": n, "Y": Y, "block": k if k < 4 else "all4", "signed": signed,
                    "npairs": int(len(P)), "J": rec["J"], "hstar": rec["hstar"],
                    "worst_pair": [int(P[i]), int(Q[i])],
                    "hpred": h_pred(len(P), surv),
                    "tail_obs": [int(sum(hist[H + 1:])) for H in range(Hcap + 1)],
                    "tail_exp": [len(P) * s for s in surv],
                })
            print(f"e2 n={n} signed={signed} h*(all4)={outs[4]['prefix'][-1]['hstar']}", flush=True)


# ---------------------------------------------------------------------------
# structured families (e3)
# ---------------------------------------------------------------------------
def _in_Sn(n, p, q):
    return p < q <= 2 * p and (1 << (n - 1)) <= p * q < (1 << n)


def _fam_chunk(args):
    fam, n, seed, count, extra = args
    rng = random.Random(f"dchworst-fam|{fam}|n={n}|{seed}")
    plo = isqrt(1 << (n - 2))
    phi = isqrt((1 << n) - 1) + 1
    out = []
    M30 = 6469693230
    tries = 0
    while len(out) < count:
        tries += 1
        if tries > 10 ** 7:
            break
        p = int(gmpy2.next_prime(rng.randint(plo, phi)))
        if fam == "twin_gap":                # q = next prime after p
            q = int(gmpy2.next_prime(p))
        elif fam == "near_2p":               # q = largest prime < 2p
            q = int(gmpy2.prev_prime(2 * p))
        elif fam == "near_3_2":              # q = next prime after 3p/2
            q = int(gmpy2.next_prime(3 * p // 2))
        elif fam == "near_5_4":
            q = int(gmpy2.next_prime(5 * p // 4))
        elif fam in ("cong_p_M30", "cong_mp_M30"):
            r = p % M30 if fam == "cong_p_M30" else (-p) % M30
            kmin = (p - r) // M30 + 1
            kmax = (2 * p - r) // M30
            if kmax < kmin:
                continue
            q = r + M30 * rng.randint(kmin, kmax)
            if not gmpy2.is_prime(q):
                continue
        elif fam == "crt_worst2":            # force the worst class at 2 band primes
            (l1, t1), (l2, t2) = extra
            L = l1 * l2
            # q == p * t^{-1} (mod l) <=> t_l = p q^{-1} == t
            r1 = p * pow(t1, -1, l1) % l1
            r2 = p * pow(t2, -1, l2) % l2
            r = (r1 * l2 * pow(l2, -1, l1) + r2 * l1 * pow(l1, -1, l2)) % L
            kmin = (p - r) // L + 1
            kmax = (2 * p - r) // L
            q = r + L * rng.randint(kmin, kmax)
            if not gmpy2.is_prime(q):
                continue
        else:
            raise ValueError(fam)
        if _in_Sn(n, p, q):
            out.append((p, q))
    return out


FAMILIES = ["twin_gap", "near_2p", "near_3_2", "near_5_4", "cong_p_M30", "cong_mp_M30",
            "crt_worst2"]


def e3():
    path = _out("e3")
    if os.path.exists(path):
        os.remove(path)
    n = 72
    Y = pow2_frac(n, 1, 6)
    b0 = blocks(Y)[0]
    Hcap = 64
    # worst classes (positive weights) at the two smallest band primes
    extra = []
    for l in b0[:2]:
        tab = cv.hmin_table(l, Hcap, False)
        t = int(np.argmax(tab[2:]) + 2)
        extra.append((l, t))
    fams = {"random": sample_pairs(n, 100000, seed0=7200)}
    for fam in FAMILIES:
        tasks = [(fam, n, s, 1000, extra) for s in range(20)]
        with Pool(20) as pool:
            ch = pool.map(_fam_chunk, tasks)
        pairs = [x for c in ch for x in c]
        fams[fam] = (np.array([p for p, _ in pairs], dtype=np.int64),
                     np.array([q for _, q in pairs], dtype=np.int64))
    for fam, (P, Q) in fams.items():
        for signed in (False, True):
            out = cv.sweep(P, Q, b0, Hcap, signed, record_every=[len(b0)])
            rec = out["prefix"][-1]
            hist = rec["hist"]
            _write(path, {"n": n, "Y": Y, "family": fam, "signed": signed,
                          "npairs": int(len(P)), "J": len(b0),
                          "tail_obs": [int(sum(hist[H + 1:])) for H in range(Hcap + 1)],
                          "tail_exp_frac": rec["poisson_surv"],
                          "hstar": rec["hstar"]})
        print(f"e3 {fam} npairs={len(P)}", flush=True)


# ---------------------------------------------------------------------------
def dprime_Y(n: int, H: int) -> int:
    return max(3, -(-pow2_frac(n, 1, 6) // max(1, icbrt(H))))


def e4():
    path = _out("e4")
    if os.path.exists(path):
        os.remove(path)
    for n in (24, 28, 32):
        P, Q = cv.all_pairs(n)
        for H in range(2, 11):
            Y = dprime_Y(n, H)
            bl = blocks(Y)
            covered_any = np.zeros(len(P), dtype=bool)
            covered_clean = np.zeros(len(P), dtype=bool)
            logsurv = 0.0
            for k, b in enumerate(bl):
                if not b:
                    continue
                out = cv.sweep(P, Q, b, H, False, record_every=[len(b)])
                cov = out["h"] <= H
                Lk = 1
                clean = np.ones(len(P), dtype=bool)
                for l in b:
                    clean &= ((Q - 1) % l) != 0
                covered_any |= cov
                covered_clean |= cov & clean
                logsurv += math.log(max(1e-300, out["prefix"][-1]["poisson_surv"][H]))
            _write(path, {"n": n, "H": H, "Y": Y, "band_sizes": [len(b) for b in bl],
                          "minkowski_regime": bool(H * H >= Y), "npairs": int(len(P)),
                          "uncovered": int((~covered_any).sum()),
                          "uncovered_clean": int((~covered_clean).sum()),
                          "poisson_expected": len(P) * math.exp(logsurv)})
            print(f"e4 n={n} H={H} Y={Y} unc={(~covered_any).sum()} "
                  f"pred={len(P)*math.exp(logsurv):.2f}", flush=True)


def e5():
    path = _out("e5")
    if os.path.exists(path):
        os.remove(path)
    n = 64
    P, Q = sample_pairs(n, 1000000, seed0=6400)
    for H in (7, 9):
        Y = dprime_Y(n, H)
        bl = blocks(Y)
        with Pool(4) as pool:
            outs = pool.map(_sweep_job, [(P, Q, b, H, False) for b in bl])
        unc = np.ones(len(P), dtype=bool)
        logsurv = 0.0
        per_block = []
        for out in outs:
            u = out["h"] > H
            per_block.append(int(u.sum()))
            unc &= u
            logsurv += math.log(out["prefix"][-1]["poisson_surv"][H])
        idx = np.nonzero(unc)[0]
        _write(path, {"n": n, "H": H, "Y": Y, "npairs": int(len(P)),
                      "uncovered_per_block": per_block,
                      "poisson_per_block": [len(P) * o["prefix"][-1]["poisson_surv"][H]
                                            for o in outs],
                      "uncovered_all4": int(unc.sum()),
                      "poisson_all4": len(P) * math.exp(logsurv),
                      "examples": [[int(P[i]), int(Q[i])] for i in idx[:20]]})
        print(f"e5 H={H} all4 uncovered {unc.sum()} vs {len(P)*math.exp(logsurv):.2f}", flush=True)


if __name__ == "__main__":
    {"e1": e1, "e2": e2, "e3": e3, "e4": e4, "e5": e5}[sys.argv[1]]()
