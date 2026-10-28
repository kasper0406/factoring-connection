"""HARNESS: experiments E1-E4 of notes/cuspcoef_prereg.md.

  PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e4|e3|e1|e2 [--quick]
"""
from __future__ import annotations

import json
import math
import os
import random
import resource
import sys
import time
from math import isqrt

import gmpy2
import numpy as np

from src.cuspcoef import decode as D
from src.cuspcoef.family import family
from src.cuspcoef import harness as H

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")


def out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, f"cuspcoef_{name}.jsonl"), "w")


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20  # macOS: bytes


def primes_upto(B):
    s = bytearray([1]) * (B + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(B) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(B + 1) if s[i]]


# ------------------------------------------------------------------ E4: 37a1 image checks

def e4(quick=False):
    f = out("e4")
    E37 = [0, 0, 1, -1, 0]
    P = primes_upto(20000 if quick else 200000)
    ells = [l for l in primes_upto(47) if l >= 3]
    data = {}
    for p in P:
        if p in (2, 37):
            continue
        data[p] = H.ap_curve(E37, p)
    for l in ells:
        pairs = {(data[p] % l, p % l) for p in data if p != l}
        rec = {"l": l, "pairs_seen": len(pairs), "pairs_total": l * (l - 1),
               "all": len(pairs) == l * (l - 1)}
        f.write(json.dumps(rec) + "\n")
        print(rec)
    viol = sum(1 for p, a in data.items() if a % 2 == 1 and gmpy2.kronecker(37, p) != 1)
    conv_fail = sum(1 for p, a in data.items() if a % 2 == 0 and gmpy2.kronecker(37, p) == 1)
    rec = {"entanglement_violations": viol, "converse_failures": conv_fail, "nprimes": len(data)}
    f.write(json.dumps(rec) + "\n")
    print(rec)


# ------------------------------------------------------------------ E3: Gekeler spread constant

def e3(quick=False):
    f = out("e3")
    rng = random.Random(3003)
    for bits in (16, 20) if quick else (16, 20, 24):
        for i in range(3 if quick else 10):
            p = H.random_prime(rng, bits)
            T = isqrt(4 * p - 1)
            vals = [H.hurwitz(4 * p - t * t) for t in range(0, T + 1)]   # 12 H
            tot = vals[0] + 2 * sum(vals[1:])
            ok = tot == 12 * 2 * p
            mx = max(vals) / 12.0
            Bemp = math.sqrt(p) * mx / (2 * p)
            bound = (2 + math.log(4 * p)) * (math.exp(0.5772156649) * math.log(math.log(4 * p)) + 1) / math.pi
            rec = {"bits": bits, "p": p, "KH_ok": ok, "B_emp": Bemp, "bound": bound,
                   "t_argmax": int(np.argmax(vals))}
            f.write(json.dumps(rec) + "\n")
            print(rec)


# ------------------------------------------------------------------ E1: products of K curves -> factoring

def e1(quick=False):
    f = out("e1" if not quick else "e1_quick")
    F = family()
    Ks = (2000, 8000, 32000) if not quick else (500, 2000)
    Kmax = Ks[-1]
    ninst = 20 if not quick else 3
    for nbits in (32, 40, 48) if not quick else (32,):
        X = 1 << (nbits // 2 - 1)
        R = isqrt(8 * X)
        t0 = time.time()
        spf = D.SPF(R * R)
        print(f"n={nbits}: SPF built in {time.time() - t0:.1f}s, rss {rss_mb():.0f} MB", flush=True)
        for i in range(ninst):
            seed = 1000 + i
            N, p, q = H.instance(seed, nbits)
            curves = F[:Kmax]
            A = [H.ap_AB(a, b, p) * H.ap_AB(a, b, q) for (a, b) in curves]
            ells = D.decoding_primes(X)
            s_true = p + q
            runs = [("full", K, A[:K]) for K in Ks]
            rng = random.Random(10 ** 6 + seed)
            A_null = [rng.randint(-R, R) * rng.randint(-R, R) for _ in range(Kmax)]
            runs.append(("null", Kmax, A_null))
            runs.append(("cong", Kmax, A))
            for mode, K, Avals in runs:
                t1 = time.time()
                res = D.decode(N, nbits, curves[:K], Avals, mode="cong" if mode == "cong" else "full", spf=spf)
                ranks = [rl.index(s_true % l) if (s_true % l) in rl else None for l, rl in zip(res["ells"], res["rank_lists"])]
                rec = {"n": nbits, "seed": seed, "mode": mode, "K": K, "success": res["p"] in (p, q),
                       "tried": res["tried"], "ells": res["ells"], "true_ranks": ranks, "n_used": res["n_used"],
                       "mean_cand": res["mean_cand"], "accidental": len(res["accidental"]),
                       "accidental_proper": [g for g in res["accidental"] if 1 < g < N],
                       "time": time.time() - t1}
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print(n_fmt(rec), flush=True)
        del spf
    print(f"peak rss {rss_mb():.0f} MB")


def n_fmt(r):
    return (f"n={r['n']} seed={r['seed']} {r['mode']:4s} K={r['K']:5d} ok={r['success']!s:5s} tried={r['tried']:4d} "
            f"ranks={r['true_ranks']} cand={r['mean_cand']} acc={r['accidental']} t={r['time']:.1f}")


# ------------------------------------------------------------------ E2: selection of easy products

def smooth_cofactor(A: int, primorial: int) -> int:
    m = abs(A)
    while True:
        g = int(gmpy2.gcd(m, primorial % m if m > 1 else 1)) if m > 1 else 1
        if g == 1:
            return m
        m //= g


def accept(A: int, primorial: int) -> bool:
    if A == 0:
        return False
    c = smooth_cofactor(A, primorial)
    return c == 1 or H.is_prime(c)


def e2(quick=False):
    f = out("e2" if not quick else "e2_quick")
    F = family()[:32000]
    nsamp = 600 if not quick else 40
    nsim = 100000 if not quick else 2000
    for nbits in (64, 128, 192, 256) if not quick else (64, 128):
        B = nbits * nbits
        primorial = 1
        for r in primes_upto(B):
            primorial *= r
        primorial = gmpy2.mpz(primorial)
        rng = random.Random(7_000_000 + nbits)
        acc_real = 0
        insts = []
        t0 = time.time()
        for i in range(nsamp):
            N, p, q = H.instance(7_000_000 + 1000 * nbits + i, nbits)
            a, b = F[rng.randrange(len(F))]
            Aval = H.ap_AB(a, b, p) * H.ap_AB(a, b, q)
            acc_real += accept(Aval, primorial)
            insts.append((p, q))
        t_real = time.time() - t0
        acc_sim = 0
        rs = random.Random(8_000_000 + nbits)
        for j in range(nsim):
            p, q = insts[j % len(insts)]
            hp, hq = 2 * isqrt(p), 2 * isqrt(q)
            Aval = rs.randint(-hp, hp) * rs.randint(-hq, hq)
            acc_sim += accept(Aval, primorial)
        rec = {"n": nbits, "B": B, "real_acc": acc_real, "real_n": nsamp, "real_rate": acc_real / nsamp,
               "sim_acc": acc_sim, "sim_n": nsim, "sim_rate": acc_sim / nsim, "t_real": t_real,
               "t_total": time.time() - t0}
        f.write(json.dumps(rec) + "\n")
        f.flush()
        print(rec, flush=True)


if __name__ == "__main__":
    exp = sys.argv[1]
    quick = "--quick" in sys.argv
    {"e1": e1, "e2": e2, "e3": e3, "e4": e4}[exp](quick)
