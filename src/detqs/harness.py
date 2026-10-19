"""HARNESS (diagnostic; uses primality of p, q for instance generation; never
imported by operational code).  Pre-registered experiments E1-E3 of
notes/detqs_prereg.md.

    PYTHONPATH=. .venv/bin/python -m src.detqs.harness e1 [n ...]
    PYTHONPATH=. .venv/bin/python -m src.detqs.harness e2 [n ...]
    PYTHONPATH=. .venv/bin/python -m src.detqs.harness summary

Parameters (frozen): y(n) = 2^floor(n/5); T(n) from FROZEN_T; samples from FROZEN_M.
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from typing import Dict, List, Tuple

import gmpy2
import numpy as np

from .localmodel import _vdist, rho, u_blocks
from .qs import qs_run
from .sieve import factor_base, primes_upto, smooth_mask

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
ANA = os.path.join(ROOT, "experiments", "analysis")

FROZEN_T = {40: 16, 48: 18, 56: 18, 64: 19, 72: 20, 80: 21}   # log2 T
E1_M = {40: (2000, 1000), 48: (2000, 1000), 56: (2000, 1000), 64: (2000, 1000), 72: (1000, 500)}
E2_M = {40: 1000, 48: 1000, 56: 1000, 64: 1000, 72: 400, 80: 200}
DMAX = 30


def y_of(n: int) -> int:
    return 1 << (n // 5)


def T_of(n: int) -> int:
    return 1 << FROZEN_T[n]


# ----------------------------------------------------------------- instances
def rand_int(n: int, i: int) -> int:
    rng = random.Random(f"detqs|int|{n}|{i}")
    return rng.randrange(1 << (n - 1), 1 << n)


def rand_semiprime(n: int, i: int, tag: str = "sp") -> Tuple[int, int, int]:
    """Uniform over S_n = {(p,q): p < q <= 2p primes, 2^{n-1} <= pq < 2^n} (rectangle rejection)."""
    rng = random.Random(f"detqs|{tag}|{n}|{i}")
    lo, hi = 1 << (n - 1), 1 << n
    A = int(gmpy2.isqrt(lo // 2))
    B = int(gmpy2.isqrt(hi)) + 1
    while True:
        p = rng.randrange(A, B)
        q = rng.randrange(A, 2 * B)
        if not (p < q <= 2 * p):
            continue
        N = p * q
        if not (lo <= N < hi):
            continue
        if gmpy2.is_prime(p) and gmpy2.is_prime(q):
            return N, p, q


# ----------------------------------------------------------------- local model, fast form
class LModel:
    """L(N) of localmodel.model_LN, with per-(l, u) caches (same numbers, faster)."""

    def __init__(self, n: int, y: int, T: int, nblocks: int = 8, E: int = 12):
        self.y, self.E = y, E
        self.logy = math.log(y)
        self.us, self.ws = u_blocks(n, y, T, nblocks)
        self.ps = primes_upto(y)
        self.g = {}   # (l, b) -> list g(e)
        self.den = {}
        for b, u in enumerate(self.us):
            r0 = rho(u)
            for l in self.ps:
                g = [rho(u - e * math.log(l) / self.logy) / r0 for e in range(E + 1)]
                self.g[(l, b)] = g
                self.den[(l, b)] = sum((1 - 1 / l) * l ** (-e) * g[e] for e in range(E)) + l ** (-E) * g[E]
        self.lam = {}  # (l, b, chi) -> log lambda for odd l not dividing N
        for b in range(len(self.us)):
            for l in self.ps:
                if l == 2:
                    continue
                for chi in (1, -1):
                    Nfake = self._rep(l, chi)
                    self.lam[(l, b, chi)] = math.log(self._lf(l, b, _vdist(l, Nfake, E)))
        self.cache2: Dict[Tuple[int, int], List[float]] = {}

    @staticmethod
    def _rep(l: int, chi: int) -> int:
        a = 1
        while True:
            c = 1 if pow(a, (l - 1) // 2, l) == 1 else -1
            if a % l and c == chi:
                return a
            a += 1

    def _lf(self, l: int, b: int, d: List[float]) -> float:
        g = self.g[(l, b)]
        return sum(d[e] * g[e] for e in range(self.E + 1)) / self.den[(l, b)]

    def L(self, N: int) -> float:
        tot = 0.0
        dl: Dict[int, List[float]] = {}
        chis: Dict[int, int] = {}
        for l in self.ps:
            if l == 2 or N % l == 0:
                EE = self.E
                while l ** EE > 131072:
                    EE -= 1
                key = (l, N % (l ** EE))
                if key not in self.cache2:
                    self.cache2[key] = _vdist(l, N, self.E)
                dl[l] = self.cache2[key]
            else:
                chis[l] = 1 if pow(N % l, (l - 1) // 2, l) == 1 else -1
        for b, w in enumerate(self.ws):
            s = 0.0
            for l, c in chis.items():
                s += self.lam[(l, b, c)]
            for l, d in dl.items():
                s += math.log(max(self._lf(l, b, d), 1e-300))
            tot += w * math.exp(s)
        return tot / sum(self.ws)


# ----------------------------------------------------------------- E1 + E3
def _autocorr(f: np.ndarray) -> np.ndarray:
    T = f.size
    F = np.fft.rfft(f.astype(np.float64), 2 * T)
    a = np.fft.irfft(F * np.conj(F), 2 * T)[:T]
    r = np.rint(a)
    assert np.max(np.abs(a - r)) < 0.01
    return r.astype(np.int64)


def run_e1(n: int) -> None:
    y, T = y_of(n), T_of(n)
    Mi, Ms = E1_M[n]
    model = LModel(n, y, T)
    out = open(os.path.join(RAW, f"detqs_e1_n{n}.jsonl"), "w")
    A = np.zeros(T, dtype=np.int64)   # sum_N autocorr (integers only, E3)
    F = np.zeros(T, dtype=np.int64)   # sum_N f
    cnt = 0
    t0 = time.time()
    jobs = [("int", i) for i in range(Mi)] + [("sp", i) for i in range(Ms)]
    for kind, i in jobs:
        if kind == "int":
            N = rand_int(n, i)
            p = q = None
        else:
            N, p, q = rand_semiprime(n, i)
        r = int(gmpy2.isqrt(N))
        if r * r == N:
            out.write(json.dumps({"n": n, "kind": kind, "i": i, "N": str(N), "square": True}) + "\n")
            continue
        fb = factor_base(N, y)
        mask, m, Q, fb = smooth_mask(N, y, T, fb)
        Z = int(mask.sum())
        rec = {"n": n, "kind": kind, "i": i, "N": str(N), "Z": Z, "fb": len(fb.primes),
               "fb_divides_N": fb.divides_N[:3], "L": model.L(N), "N_mod8": N % 8}
        if kind == "sp":
            rec["blum"] = bool(p % 4 == 3 and q % 4 == 3)
        out.write(json.dumps(rec) + "\n")
        if kind == "int":
            f = mask.astype(np.int64)
            A += _autocorr(f)
            F += f
            cnt += 1
    out.close()
    # E3 aggregate
    G = _autocorr(F)  # sum_t F(t) F(t+k)
    np.savez_compressed(os.path.join(RAW, f"detqs_e3_n{n}.npz"), A=A, G=G, cnt=cnt, T=T)
    print(f"E1/E3 n={n} done in {time.time()-t0:.1f}s", flush=True)


# ----------------------------------------------------------------- E2
def run_e2(n: int) -> None:
    y, T = y_of(n), T_of(n)
    out = open(os.path.join(RAW, f"detqs_e2_n{n}.jsonl"), "w")
    t0 = time.time()
    for i in range(E2_M[n]):
        N, p, q = rand_semiprime(n, i, tag="e2")
        res = qs_run(N, y, T, dmax=DMAX)
        blum = bool(p % 4 == 3 and q % 4 == 3)
        cls = "none"
        if res.factor is not None:
            cls = "STRUCTURAL" if res.how == "dependency" else "ACCIDENTAL"
            assert res.factor in (p, q) or res.how == "square"
        jac_match = None
        if blum and res.dep_outcomes:
            jac_match = all((j == -1) == bool(o) for j, o in zip(res.dep_jacobi, res.dep_outcomes))
        rec = {"n": n, "i": i, "N": str(N), "blum": blum, "Z": res.Z, "fb": res.fb_size,
               "how": res.how, "class": cls, "first_fruitful": res.first_fruitful,
               "deps_tested": res.deps_tested, "relations_used": res.relations_used,
               "dep_sizes": res.dep_sizes, "outcomes": res.dep_outcomes, "jacobi_match": jac_match}
        out.write(json.dumps(rec) + "\n")
    out.close()
    print(f"E2 n={n} done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(ANA, exist_ok=True)
    cmd = sys.argv[1]
    ns = [int(a) for a in sys.argv[2:]]
    if cmd == "e1":
        for n in ns or sorted(E1_M):
            run_e1(n)
    elif cmd == "e2":
        for n in ns or sorted(E2_M):
            run_e2(n)
    elif cmd == "summary":
        from .analysis import summarize
        summarize()
