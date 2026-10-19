"""HARNESS: exact version of E1 (Addendum A of notes/hypD_prereg.md).

For small N, enumerate all (u, v) in S x S, form R = u/v mod N, intersect with the squares Q_N
(membership via Legendre symbols mod the known factors: harness-only), and compare
|good R| / |Q_N| with the coprime-pair model |S|^2/phi(N) and the Theorem S bound.
G modes: "random" (each prime of J kept with prob g), "qr" (keep only l in Q_N; then every
u/v is a square: the model predicts a factor 2^omega excess).
Usage: .venv/bin/python -m src.hypD.e1exact CONFIG.json OUT.jsonl [workers]
"""
from __future__ import annotations

import json
import random
import sys
from math import isqrt
from multiprocessing import Pool

import gmpy2

from ..nfsfruit.instances import composite
from ..nfsfruit.smooth import primes_upto
from . import split as sp


def S_window_list(G, k, Y):
    Gs = sorted(G)
    out = []

    def rec(start, depth, prod):
        if depth == k:
            if Y // 2 <= prod <= Y:
                out.append(prod)
            return
        for i in range(start, len(Gs)):
            q = prod * Gs[i]
            if q * (Gs[i] ** (k - depth - 1)) > Y:
                break
            rec(i + 1, depth + 1, q)
    rec(0, 0, 1)
    return out


def job(args):
    cls, bits, seed, k, mode, g = args
    N, fs = composite(cls, bits, seed)
    rng = random.Random(f"hypD-e1x|{cls}|{bits}|{seed}|{k}|{mode}|{g}")
    Y = isqrt(N // 2)
    lo = 1 << ((Y.bit_length() - 1) // k)
    J = [p for p in primes_upto(2 * lo) if p > lo and N % p]
    if mode == "qr":
        G = [p for p in J if all(gmpy2.legendre(p, q) == 1 for q in fs)]
    else:
        G = [p for p in J if rng.random() < g]
    S = S_window_list(G, k, Y)
    inv = {v: pow(v, -1, N) for v in S}
    good = set()
    for u in S:
        for v in S:
            good.add(u * inv[v] % N)
    nsq = sum(1 for R in good if all(gmpy2.legendre(R, q) == 1 for q in fs))
    phiN = 1
    for p in fs:
        phiN *= p - 1
    QN = phiN // (2 ** len(fs))
    exact = nsq / QN
    return dict(exp="e1x", cls=cls, bits=bits, seed=seed, k=k, mode=mode, g=g, N=N,
                omega=len(fs), S=len(S), G=len(G), exact=exact, model=len(S) ** 2 / phiN,
                bound=sp.theorem_s_bound(len(S), Y, N, len(fs)))


def main():
    cfg = json.load(open(sys.argv[1]))
    jobs = [(c, b, s, k, mo, g) for c in cfg["classes"] for b in cfg["bits"]
            for s in range(cfg["seeds"]) for k in cfg["k"] for (mo, g) in cfg["modes"]]
    with Pool(int(sys.argv[3]) if len(sys.argv) > 3 else 4) as pool, open(sys.argv[2], "w") as out:
        for rec in pool.imap_unordered(job, jobs):
            out.write(json.dumps(rec) + "\n")
            out.flush()


if __name__ == "__main__":
    main()
