"""HARNESS for nativeana experiments E1-E3 (notes/nativeana_prereg.md).  Uses PARI (via cypari2) for ground truth,
which factors discriminants internally; never imported by operational code.

Usage:  PYTHONPATH=. .venv/bin/python -m src.nativeana.harness E1|E2|E3|summary
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import gmpy2
import numpy as np

from src.common.instances import balanced_semiprime
from src.nativeana import chi as C

RAW = Path("experiments/raw")
_pari = None


def pari():
    global _pari
    if _pari is None:
        import cypari2
        _pari = cypari2.Pari()
        _pari.allocatemem(2 * 10 ** 8)
    return _pari


def class_number(D: int) -> int:
    """h(D) for D < 0 (order class number).  qfbclassno, cross-checked with quadclassunit for |D| > 2e10."""
    P = pari()
    h = int(P.qfbclassno(D))
    if abs(D) > 2 * 10 ** 10:
        h2 = int(P(f"quadclassunit({D}).no"))
        assert h == h2, (D, h, h2)
    return h


def hurwitz(n: int) -> Fraction:
    return Fraction(str(pari().qfbhclassno(n)))


def instances(n: int, count: int, mod4: int | None):
    out, skipped, seed = [], [], 0
    while len(out) < count:
        inst = balanced_semiprime(n, seed)
        if mod4 is None or inst.N % 4 == mod4:
            out.append(inst)
        else:
            skipped.append(seed)
        seed += 1
    return out, skipped


def slope(xs, ys):
    lx = np.log(np.asarray(xs, float))
    ly = np.log(np.asarray(ys, float))
    return float(np.polyfit(lx, ly, 1)[0])


# ------------------------------------------------------------------------------------------------ E1
def run_E1(sizes=(28, 36, 44), count=20):
    path = RAW / "nativeana_E1.jsonl"
    with open(path, "w") as f:
        for n in sizes:
            insts, skipped = instances(n, count, 1)
            f.write(json.dumps({"kind": "skips", "n": n, "skipped_seeds": skipped}) + "\n")
            for inst in insts:
                t0 = time.time()
                N = inst.N
                D = -4 * N
                Q = 4 * N
                h = class_number(D)
                sq = math.isqrt(Q)
                afe_top = int(6.6 * math.sqrt(Q / math.pi)) + 2      # erfc(6.6) < 1e-19
                Tafe = [1 << j for j in range(4, (8 * sq).bit_length() + 1)]
                Tps = [1 << j for j in range(4, min(20, int(math.log2(sq / 4))) + 1)]
                xs = [1 << j for j in range(4, 21)]
                Tmax = max(afe_top, 6 * max(Tps), max(xs))
                chi = C.chi_table(D, Tmax)
                ep = C.euler_product_h(D, chi, xs)
                ps = C.partial_sum_h(D, chi, Tps)
                sm = C.smooth_sum_h(D, chi, Tps)
                af = C.afe_h(D, chi[:afe_top + 1], Tafe)
                rs = C.random_summation_h(D, [1 << j for j in range(6, 17)], f"nativeana|RS|{N}")
                rec = {"kind": "inst", "n": n, "seed": inst.seed, "N": N, "p": inst.factors[0], "q": inst.factors[1],
                       "h": h, "sqrtQ": math.sqrt(Q),
                       "EP": {x: abs(v - h) for x, v in ep.items()},
                       "PS": {T: abs(v - h) for T, v in ps.items()},
                       "SM": {T: abs(v - h) for T, v in sm.items()},
                       "AFE": {T: abs(v - h) for T, v in af.items()},
                       "RS_B": rs["B"], "RS": {K: abs(v - h) for K, v in rs["est"].items()},
                       "secs": time.time() - t0}
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print(n, inst.seed, f"{rec['secs']:.1f}s", flush=True)
                del chi


# ------------------------------------------------------------------------------------------------ E2
def chi_features(D: int, ns):
    return [int(gmpy2.kronecker(D, n)) for n in ns]


def run_E2(Q=1 << 30, ntrain=6000, ntest=3000):
    path = RAW / "nativeana_E2.jsonl"
    Kmax = 512
    odd_ns = list(range(1, 2 * Kmax, 2))
    odd_primes = [int(x) for x in C.primes_upto(4000) if x > 2][:520]
    total = ntrain + ntest
    fams = {}
    # family A: any m = 1 mod 4 in [Q, 2Q]
    ms = []
    i = 0
    while len(ms) < total:
        r = random.Random(f"nativeana|E2|A|{i}")
        m = Q + r.randrange(Q)
        m -= (m - 1) % 4
        if m >= Q:
            ms.append(m)
        i += 1
    fams["A"] = (ms, [])
    # family B: balanced 31-bit semiprimes, N = 1 mod 4, N in [Q, 2Q]
    ms, skipped, seed = [], [], 0
    while len(ms) < total:
        inst = balanced_semiprime(31, seed)
        if inst.N % 4 == 1 and Q <= inst.N <= 2 * Q:
            ms.append(inst.N)
        else:
            skipped.append(seed)
        seed += 1
    fams["B"] = (ms, skipped)
    with open(path, "w") as f:
        for fam, (ms, skipped) in fams.items():
            t0 = time.time()
            X = np.array([chi_features(-4 * m, odd_ns) for m in ms], dtype=np.float64)
            Xp = np.array([chi_features(-4 * m, odd_primes) for m in ms], dtype=np.float64)
            L = np.array([math.pi * class_number(-4 * m) / (2 * math.sqrt(m)) for m in ms])
            u = (np.array(ms, dtype=np.float64) - Q) / Q
            tr, te = slice(0, ntrain), slice(ntrain, total)
            res = {"kind": "family", "fam": fam, "skipped": len(skipped), "meanL": float(L.mean()),
                   "sdL": float(L.std()), "rules": {}}
            for K in [2 ** j for j in range(1, 10)]:
                F = X[:, :K]
                one = np.ones((total, 1))
                for rule, A in (("i", np.hstack([one, F])), ("ii", np.hstack([one, u[:, None], F, F * u[:, None]]))):
                    w, *_ = np.linalg.lstsq(A[tr], L[tr], rcond=None)
                    rms_te = float(np.sqrt(np.mean((A[te] @ w - L[te]) ** 2)))
                    rms_tr = float(np.sqrt(np.mean((A[tr] @ w - L[tr]) ** 2)))
                    res["rules"][f"{rule}_{K}"] = {"test": rms_te, "train": rms_tr}
            logL = np.log(L)
            for K in [2, 4, 8, 16, 32, 64, 128, 256, 512]:
                Fp = Xp[:, :K]
                A = np.hstack([np.ones((total, 1)), Fp, Fp ** 2])
                w, *_ = np.linalg.lstsq(A[tr], logL[tr], rcond=None)
                res["rules"][f"log_{K}"] = {"test": float(np.sqrt(np.mean((A[te] @ w - logL[te]) ** 2))),
                                            "train": float(np.sqrt(np.mean((A[tr] @ w - logL[tr]) ** 2))),
                                            "x": odd_primes[K - 1]}
            res["secs"] = time.time() - t0
            f.write(json.dumps(res) + "\n")
            f.flush()
            print(fam, f"{res['secs']:.1f}s", flush=True)


# ------------------------------------------------------------------------------------------------ E3
def det_counts(ell: int, k: int):
    """#{M in M_2(Z/ell^k): det M = a} for all a, by exact convolution of product counts."""
    m = ell ** k
    x = np.arange(m)
    prod = (x[:, None] * x[None, :]) % m
    P = np.bincount(prod.ravel(), minlength=m).astype(np.int64)   # #{(a,d): ad = r}
    # det = ad - bc: count[a] = sum_r P[r] * P[r - a]
    cnt = np.array([int(np.dot(P, np.roll(P, a))) for a in range(m)])
    return cnt


def run_E3(sizes=(16, 20, 24), count=20):
    path = RAW / "nativeana_E3.jsonl"
    gl_u, gl_w = np.polynomial.legendre.leggauss(4000)
    wbump = lambda u: np.where(np.abs(u) < 1, np.exp(-1.0 / np.maximum(1e-300, 1 - u * u)), 0.0)
    with open(path, "w") as f:
        # P3.2
        checks = []
        for ell, k in ((3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1), (7, 2)):
            if ell ** (4 * k) > 2 * 10 ** 7:
                continue
            m = ell ** k
            cnt = det_counts(ell, k)
            bad = 0
            tested = 0
            for a in range(m):
                v = 0
                aa = a
                if aa == 0:
                    continue  # v >= k
                while aa % ell == 0:
                    aa //= ell
                    v += 1
                pred = Fraction(ell ** (3 * k)) * (1 - Fraction(1, ell * ell)) * sum(Fraction(1, ell ** j) for j in range(v + 1))
                tested += 1
                bad += int(pred != cnt[a])
            checks.append({"ell": ell, "k": k, "tested": tested, "bad": bad})
        f.write(json.dumps({"kind": "P3.2", "checks": checks}) + "\n")
        print("P3.2", checks, flush=True)
        for n in sizes:
            insts, _ = instances(n, count, None)
            for inst in insts:
                t0 = time.time()
                N = inst.N
                p, q = inst.factors
                tmax = math.isqrt(4 * N - 1)
                if tmax * tmax == 4 * N:
                    tmax -= 1
                S0 = S2 = S4 = Fraction(0)
                Sw = Aw = 0.0
                sqN = math.sqrt(N)
                for t in range(-tmax, tmax + 1):
                    Hn = hurwitz(4 * N - t * t)
                    S0 += Hn
                    S2 += (t * t - N) * Hn
                    S4 += (t ** 4 - 3 * t * t * N + N * N) * Hn
                    u = t / (2 * sqN)
                    wu = float(wbump(np.array(u)))
                    Sw += wu * float(Hn)
                    Aw += wu * math.sqrt(4 * N - t * t) / math.pi
                sigma = (1 + p) * (1 + q)
                Aw *= sigma / N
                r = math.sqrt(p / q)
                Jw = float(2 / math.pi * np.sum(gl_w * wbump(gl_u) * np.sqrt(1 - gl_u ** 2) / (1 - 2 * gl_u * r + r * r)))
                rec = {"kind": "inst", "n": n, "seed": inst.seed, "N": N, "p": p, "q": q,
                       "P3.1": [S0 == 2 * (N + q), S2 == -2 * (1 + p ** 3), S4 == -2 * (1 + p ** 5)],
                       "Sw": Sw, "Aw": Aw, "Rw": Sw - Aw, "Rw_over_sqrtN": (Sw - Aw) / sqN,
                       "hyp_pred": -2 * p * Jw, "secs": time.time() - t0}
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print(n, inst.seed, rec["P3.1"], f"{rec['Rw_over_sqrtN']:.3f}", f"{rec['hyp_pred']/sqN:.3f}", flush=True)


if __name__ == "__main__":
    {"E1": run_E1, "E2": run_E2, "E3": run_E3}[sys.argv[1]]()
