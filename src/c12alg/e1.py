"""c12alg E1 (pre-registered in notes/c12alg_prereg.md).  DIAGNOSTIC code, no N is factored.

E1a  exact total variation of Lemma TP: Law(Phi(c)) vs Law(Phi(c) + z), z uniform on R_{d-2,d-1}
     (= the explicit measure pi * Unif(Gamma cap Pi) of Lemma TP), against the bound d(1 + 4V^2)/(2H + 1).
E1b  Monte Carlo joint smoothness of (F(P1), F(P2)) under design DW for fixed g, against the product of
     marginals, the lattice gcd-class model J_gcd and the archimedean-corrected model J_full.

Usage:
  .venv/bin/python -m src.c12alg.e1 E1a experiments/raw/c12alg_e1a.jsonl
  .venv/bin/python -m src.c12alg.e1 E1b experiments/raw/c12alg_e1b.jsonl
"""
from __future__ import annotations

import json
import random
import sys
import time
from fractions import Fraction
from math import gcd, log
from typing import Dict, List, Tuple

import numpy as np

from . import dw
from src.thmA.numfield import hom_eval


# ----------------------------------------------------------------------------- E1a
def _tv_int(k1, c1, n1, k2, c2, n2) -> Fraction:
    allk = np.union1d(k1, k2)
    a = np.zeros(len(allk), dtype=np.int64)
    b = np.zeros(len(allk), dtype=np.int64)
    a[np.searchsorted(allk, k1)] = c1
    b[np.searchsorted(allk, k2)] = c2
    num = int(np.abs(a * n2 - b * n1).sum())          # |a n2 - b n1| < 2^63 for the sizes used
    return Fraction(num, 2 * n1 * n2)


def _rand_pair(rng: random.Random, V: int) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    while True:
        P1 = (rng.randint(-V, V), rng.randint(1, V))
        P2 = (rng.randint(-V, V), rng.randint(1, V))
        if (gcd(*P1) == 1 and gcd(*P2) == 1 and dw.det2(P1, P2) != 0 and max(map(abs, P1 + P2)) == V
                and P1[1] * P2[1] >= 2):
            return P1, P2


E1A_CONFIGS = [  # (d, V, H)
    (3, 2, 40), (3, 2, 80), (3, 2, 120),
    (3, 3, 80), (3, 3, 120),
    (4, 2, 40),
]
E1A_PAIRS = 12


def run_e1a(out: str, seed: int = 20261024):
    rng = random.Random(f"c12alg|E1a|{seed}")
    with open(out, "w") as fh:
        for (d, V, H) in E1A_CONFIGS:
            for i in range(E1A_PAIRS):
                P1, P2 = _rand_pair(rng, V)
                t0 = time.time()
                k1, c1 = dw.law_phi(P1, P2, d, H)
                n1 = (2 * H + 1) ** d
                R = dw.representatives(P1, P2, d)
                k2, c2 = dw.shift_mix(k1, c1, R)
                n2 = n1 * len(R)
                tv = _tv_int(k1, c1, n1, k2, c2, n2)
                xmax = max(max(abs(x) for x in dw.short_preimage(P1, P2, d, z)) for z in R)
                bound = dw.tv_lemma_bound(P1, P2, d, H)
                rec = {"d": d, "V": V, "H": H, "P1": P1, "P2": P2, "Delta": abs(dw.det2(P1, P2)),
                       "nR": len(R), "tv": float(tv), "tv_num": tv.numerator, "tv_den": tv.denominator,
                       "bound": float(bound), "xmax": xmax, "xmax_bound": 1 + 4 * V * V,
                       "seconds": round(time.time() - t0, 2)}
                fh.write(json.dumps(rec) + "\n"); fh.flush()
                print("E1a", d, V, H, P1, P2, float(tv), float(bound), flush=True)


# ----------------------------------------------------------------------------- E1b
def primes_upto(n: int) -> List[int]:
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def smooth_mask(x: np.ndarray, primes: List[int]) -> np.ndarray:
    """Exact y-smoothness of |x| (x != 0), vectorised trial division (int64)."""
    r = np.abs(x).astype(np.int64)
    for p in primes:
        m = (r % p) == 0
        while m.any():
            r[m] //= p
            m = (r % p) == 0
            m &= r > 0
    return (r == 1)


def dickman_table(umax: float = 12.0, h: float = 1e-3) -> Tuple[np.ndarray, float]:
    """Dickman rho on a grid (Euler on u rho'(u) = -rho(u-1)); floating point, used only for the archimedean proxy."""
    n = int(umax / h) + 1
    rho = np.ones(n)
    k = int(1 / h)
    for i in range(k + 1, n):
        u = i * h
        rho[i] = rho[i - 1] - h * rho[i - 1 - k] / u
    return rho, h


E1B = {"d": 3, "M": 8, "V": 3, "Wm": 2 ** 14, "H": 2 ** 19, "y": 512, "samples": 2_000_000,
       "fields": 6, "filter_div": 8, "chunk": 250_000}


def toy_points(V: int) -> List[Tuple[int, int]]:
    return [(a, b) for b in range(1, V + 1) for a in range(-V, V + 1) if gcd(a, b) == 1]


def run_e1b(out: str, seed: int = 20261024):
    cfg = E1B
    d, M, V, Wm, H, y = cfg["d"], cfg["M"], cfg["V"], cfg["Wm"], cfg["H"], cfg["y"]
    rng = random.Random(f"c12alg|E1b|{seed}")
    nrng = np.random.default_rng(int.from_bytes(__import__("hashlib").sha256(f"c12alg|E1b|np|{seed}".encode()).digest()[:8], "big"))
    primes = primes_upto(y)
    rho, hstep = dickman_table()
    pts = toy_points(V)
    with open(out, "w") as fh:
        for fi in range(cfg["fields"]):
            # toy N ~ M^(d+1), coprime to small primes; g primitive in [M/2, M]^2 with gcd(g1, N) = 1
            while True:
                N = rng.randrange(2 ** 12, 2 ** 13) | 1
                if all(N % p for p in primes[:20]):
                    break
            while True:
                g0, g1 = rng.randint(M // 2, M), rng.randint(max(2, M // 2), M)
                if gcd(g0, g1) == 1 and gcd(g1, N) == 1:
                    break
            g = (g0, g1)
            a = dw.fhat(N, d, g0, g1)
            keep = []
            for P in pts:
                t = abs(dw.G_of(g, P))
                if cfg["filter_div"] * t >= M * max(abs(P[0]), P[1]) and t > 0 and gcd(t, N) == 1:
                    keep.append(P)
            fh_vals = np.array([hom_eval(a, P[0], P[1]) for P in keep], dtype=np.int64)
            Gs = np.array([dw.G_of(g, P) for P in keep], dtype=np.int64)
            Vs = np.array([dw.v_of(P, d) for P in keep], dtype=np.int64)          # (k, d)
            k = len(keep)
            S1 = np.zeros(k); S12 = np.zeros((k, k))
            # gcd-class statistics per pair are accumulated on the fly: store per pair per delta | Delta
            dets = {(i, j): abs(dw.det2(keep[i], keep[j])) for i in range(k) for j in range(i + 1, k)}
            divs = {D: [e for e in range(1, D + 1) if D % e == 0] for D in set(dets.values())}
            cls_cnt = {}          # (i, D, delta) -> [count, smooth count]
            Zs = np.zeros(k); Z12 = np.zeros((k, k))
            n_done = 0
            while n_done < cfg["samples"]:
                n = min(cfg["chunk"], cfg["samples"] - n_done)
                w = nrng.integers(Wm + 1, 2 * Wm + 1, size=n, dtype=np.int64)
                c = nrng.integers(-H, H + 1, size=(n, d), dtype=np.int64)
                F = w[:, None] * fh_vals[None, :] + Gs[None, :] * (c @ Vs.T)        # (n, k)
                assert np.abs(F).max() < 2 ** 62
                sm = smooth_mask(F.ravel(), primes).reshape(F.shape) & (F != 0)
                smf = sm.astype(np.float64)
                S1 += smf.sum(0)
                S12 += smf.T @ smf
                # archimedean proxy Z = rho(log|F| / log y)
                u = np.log(np.maximum(np.abs(F), 2).astype(np.float64)) / log(y)
                Z = rho[np.minimum((u / hstep).astype(np.int64), len(rho) - 1)]
                Zs += Z.sum(0)
                Z12 += Z.T @ Z
                for D, ds in divs.items():
                    gD = np.gcd(F, D)
                    for i in range(k):
                        for e in ds:
                            m = gD[:, i] == e
                            key = (i, D, e)
                            cc = cls_cnt.setdefault(key, [0, 0])
                            cc[0] += int(m.sum()); cc[1] += int(sm[m, i].sum())
                n_done += n
            ns = cfg["samples"]
            for (i, j), D in dets.items():
                M1, M2, J = S1[i] / ns, S1[j] / ns, S12[i, j] / ns
                # J_gcd = sum_delta phi(D/delta)/D * sbar_1(delta) sbar_2(delta)
                Jg = 0.0
                for e in divs[D]:
                    c1, c2 = cls_cnt[(i, D, e)], cls_cnt[(j, D, e)]
                    s1 = c1[1] / c1[0] if c1[0] else 0.0
                    s2 = c2[1] / c2[0] if c2[0] else 0.0
                    phi = sum(1 for x in range(1, D // e + 1) if gcd(x, D // e) == 1)
                    Jg += phi / D * s1 * s2
                EZ1, EZ2, EZ12 = Zs[i] / ns, Zs[j] / ns, Z12[i, j] / ns
                arch = EZ12 / (EZ1 * EZ2) - 1.0
                rec = {"field": fi, "N": N, "g": g, "P1": keep[i], "P2": keep[j], "Delta": D,
                       "M1": M1, "M2": M2, "J": J, "J_gcd": Jg, "arch": arch, "J_full": Jg * (1 + arch),
                       "samples": ns}
                fh.write(json.dumps(rec) + "\n")
            fh.flush()
            print("E1b field", fi, "points", k, "pairs", len(dets), flush=True)


if __name__ == "__main__":
    which, out = sys.argv[1], sys.argv[2]
    {"E1a": run_e1a, "E1b": run_e1b}[which](out)
