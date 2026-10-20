"""HARNESS: experiments E1-E4 of notes/periodfind_prereg.md.

Usage: .venv/bin/python -m src.periodfind.run_experiments {e1b,e1,e2,e3,e4} [--small]
Raw JSONL: experiments/raw/periodfind_<e>.jsonl.

Orbit shortcuts (harness speed only; identical values, checked in tests/test_periodfind.py):
  Jacobi(a^t) = Jacobi(a)^t;  (a^t)^{-1} = a^{(r-t) mod r};  [a^{tM} = 1] = [r/gcd(r,M) | t].
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time

import gmpy2
import numpy as np

from src.common.instances import balanced_semiprime
from src.periodfind import observables as ob
from src.periodfind import spectra as sp
from src.periodfind import generic as gen
from src.periodfind.sft import CountingOracle, sft, recover_order

RAW = "experiments/raw"


def pick_bases(N, p, q, n, rng, count=3, min_shift=6, max_r=None, max_tries=20000):
    out = []
    tries = 0
    while len(out) < count and tries < max_tries:
        tries += 1
        a = rng.randrange(2, N - 1)
        if math.gcd(a, N) != 1:
            continue
        info = sp.order_info(a, p, q)
        if info["r"] >= 2 ** (n - min_shift) and (max_r is None or info["r"] <= max_r):
            out.append((a, info))
    return out


def panel_values(name, a, N, r, orb, rng):
    t = np.arange(r, dtype=np.int64)
    ja = gmpy2.jacobi(a, N)
    if name == "P5_jacobi":
        return np.where((t & 1) & (ja == -1), -1.0, 1.0)
    if name == "P7_jacobi_interval":
        return np.where((t & 1) & (ja == -1), -1.0, 1.0) * ob.obs_interval(orb, N)
    if name == "P9_lsb_inverse":
        return ob.obs_lsb(orb[(-t) % r], N)
    if name == "P11_oa":
        d = r // math.gcd(r, ob.lcm_upto(30))
        return (t % d == 0).astype(np.float64)
    return ob.evaluate(name, orb, N, rng)


def weil_bound(name, N, r):
    if name in ("P2_additive", "P6_shifted_jacobi"):
        return N / r ** 2  # |F^(j)| <= sqrt(N)/r, |F| = 1 a.e.
    return None


def e1(small=False, bits_list=None, mode="w", seeds=None):
    bits_list = bits_list or ([16] if small else [22, 26])
    seeds = seeds or (range(2) if small else range(10))
    path = f"{RAW}/periodfind_e1{'_small' if small else ''}.jsonl"
    with open(path, mode) as fh:
        for n in bits_list:
            for seed in seeds:
                inst = balanced_semiprime(n, seed)
                N = inst.N
                p, q = inst.factors
                rng = random.Random(f"periodfind-e1|{n}|{seed}")
                # DEVIATION (disclosed): at n = 26 also r <= 2^23, to respect the 2 GB RSS cap;
                # seeds where no such base exists within 20000 tries are recorded as skipped.
                bases = pick_bases(N, p, q, n, rng, max_r=(2 ** 23 if n >= 26 else None))
                if len(bases) < 3:
                    fh.write(json.dumps({"n": n, "seed": seed, "N": N, "skipped": True,
                                         "bases_found": len(bases)}) + "\n")
                for a, info in bases:
                    r = info["r"]
                    orb = sp.orbit(a, N, r)
                    ja = gmpy2.jacobi(a, N)
                    for name in ob.PANEL:
                        t0 = time.time()
                        F = panel_values(name, a, N, r, orb, rng)
                        rho_H = float(np.mean(F > 0)) if name == "P10_smooth" else None  # post hoc
                        spec = sp.spectrum(F)
                        del F
                        met = sp.metrics(spec, r, ja)
                        del spec
                        rec = {"n": n, "seed": seed, "N": N, "a": a, "r": r, "i": (p - 1) * (q - 1) // r,
                               "g": info["g"], "obs": name, **met, "rho_H": rho_H, "weil_bound": weil_bound(name, N, r),
                               "sec": round(time.time() - t0, 2)}
                        fh.write(json.dumps(rec) + "\n")
                        fh.flush()
                    del orb
                    print(n, seed, a, r, flush=True)


def e1b():
    path = f"{RAW}/periodfind_e1b.jsonl"
    rng = random.Random("periodfind-e1b")
    pairs = [(23, 47), (29, 53), (31, 59), (37, 43), (11, 101), (41, 47)]
    with open(path, "w") as fh:
        for p, q in pairs:
            N = p * q
            for _ in range(3):
                while True:
                    a = rng.randrange(2, N)
                    if math.gcd(a, N) == 1:
                        break
                for name in ob.PANEL:
                    if name == "P1_comb":
                        gfun = lambda y, N, c=pow(a, 5, N): ob.obs_comb(y, N, c)
                    elif name == "P2_additive":
                        gfun = lambda y, N: ob.obs_additive(y, N, 7)
                    else:
                        gfun = (lambda nm: (lambda y, N: ob.evaluate(nm, y, N, None)))(name)
                    if name == "P10_smooth":
                        gfun = lambda y, N: ob.obs_smooth(y, N, centre=False)
                    err = sp.transfer_check(N, p, q, a, gfun)
                    fh.write(json.dumps({"N": N, "p": p, "q": q, "a": a, "obs": name, "max_err": err}) + "\n")
    print("e1b done")


def e2():
    path = f"{RAW}/periodfind_e2.jsonl"
    n = 20
    with open(path, "w") as fh:
        for seed in range(10):
            inst = balanced_semiprime(n, seed)
            N = inst.N
            p, q = inst.factors
            rng = random.Random(f"periodfind-e2|{seed}")
            a, info = pick_bases(N, p, q, n, rng, count=1)[0]
            r = info["r"]
            orb = sp.orbit(a, N, r)
            ja = gmpy2.jacobi(a, N)
            m = 2 * N.bit_length() + 2
            t = np.arange(r)
            h = 1
            while True:
                h = rng.randrange(1, N)
                if math.gcd(h, N) == 1:
                    break
            tables = {
                "control_tone": np.exp(2j * np.pi * t / r),  # HARNESS-built e(log_a(y)/r)
                "jacobi": np.where((t & 1) & (ja == -1), -1.0, 1.0),
                "additive": ob.obs_additive(orb, N, h),
                "interval": ob.obs_interval(orb, N),
                "comb": ob.obs_comb(orb, N, int(orb[rng.randrange(r)])),
            }
            for fname, tab in tables.items():
                for tau in (0.25, 1 / 16):
                    rr = np.uint64(r)
                    oracle = CountingOracle(lambda x, tab=tab, rr=rr: tab[(x % rr).astype(np.int64)])
                    res = sft(oracle, m, tau, seed=hash((seed, fname, tau)) & 0xFFFF)
                    Q = 1 << m
                    rec_r = recover_order(res["out"], Q, a, N) if res["out"] else 0
                    near = [round(k * r / Q, 3) for k in res["out"]]
                    fh.write(json.dumps({"seed": seed, "N": N, "a": a, "r": r, "jac_a": ja, "f": fname,
                                         "tau": tau, "m": m, "out_k": res["out"], "out_jr": near,
                                         "queries": res["queries"], "aborted": res["aborted"],
                                         "norm2_est": res["norm2_est"], "recovered_r": rec_r,
                                         "recovered_ok": rec_r == r}) + "\n")
                    fh.flush()
                    print(seed, fname, tau, res["out"][:3], res["queries"], rec_r == r, flush=True)


def e3():
    path = f"{RAW}/periodfind_e3.jsonl"
    n = 18
    with open(path, "w") as fh:
        for seed in range(6):
            inst = balanced_semiprime(n, seed)
            N = inst.N
            p, q = inst.factors
            rng = random.Random(f"periodfind-e3|{seed}")
            a, info = pick_bases(N, p, q, n, rng, count=1)[0]
            r, rp, rq, g = info["r"], info["r_p"], info["r_q"], info["g"]
            orb = sp.orbit(a, N, r)
            divs = [d for d in range(1, 13) if rp % d == 0] + [d for d in range(1, 13) if rq % d == 0]
            exps = ([rng.randrange(2, N) for _ in range(3)] + [rp // d for d in divs if rp % d == 0][:3]
                    + [rq // d for d in divs if rq % d == 0][:3] + [r // d for d in range(1, 7) if r % d == 0])
            exps = [e for e in exps if e > 0]
            for k in range(40):
                tree = gen.random_tree(3, rng, N, exps)
                res = gen.check_tree(tree, orb, N, p, q, r, g)
                fh.write(json.dumps({"seed": seed, "N": N, "r": r, "g": g, "tree": k, **res}) + "\n")
    print("e3 done")


def e4():
    path = f"{RAW}/periodfind_e4.jsonl"
    n = 22
    with open(path, "w") as fh:
        for seed in range(5):
            inst = balanced_semiprime(n, seed)
            N = inst.N
            p, q = inst.factors
            rng = random.Random(f"periodfind-e4|{seed}")
            a, info = pick_bases(N, p, q, n, rng, count=1)[0]
            r = info["r"]
            orb = sp.orbit(a, N, r)
            while True:
                h = rng.randrange(1, N)
                if math.gcd(h, N) == 1:
                    break
            F = ob.obs_additive(orb, N, h)
            Fh = np.fft.fft(F) / r
            P = np.abs(Fh) ** 2
            P[0] = 0
            jstar = int(np.argmax(P))
            nprng = np.random.default_rng(seed)
            for ratio in (1 / 16, 1 / 4, 1, 4):
                mm = max(1, int(ratio * r))
                trials = 40

                def est(j):
                    tt = nprng.integers(0, r, size=mm)
                    return abs(np.mean(F[tt] * np.exp(-2j * np.pi * ((j * tt) % r) / r))) ** 2
                null = sorted(est(int(nprng.integers(1, r))) for _ in range(200))
                thr = null[int(0.99 * len(null)) - 1]
                det = sum(est(jstar) > thr for _ in range(trials)) / trials
                fh.write(json.dumps({"seed": seed, "N": N, "r": r, "ratio": ratio, "m": mm,
                                     "peak_rel_mass": float(P[jstar]), "r_times_peak": float(P[jstar] * r),
                                     "null_q99": thr, "detect_rate": det}) + "\n")
                fh.flush()
                print(seed, ratio, det, flush=True)


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    which = sys.argv[1]
    small = "--small" in sys.argv
    {"e1": lambda: e1(small), "e1_26": lambda: e1(False, [26], "a", seeds=range(int(sys.argv[2]), 10)), "e1b": e1b, "e2": e2, "e3": e3, "e4": e4}[which]()
