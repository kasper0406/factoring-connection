"""E3 (DIAGNOSTIC census over primes): Legendre-symbol biases of compressed pure values (Theorem PT-J).

For an N-independent integer c given by a straight-line program, a Jacobi test (c|N) = (c|p)(c|q) on the population
N = 1 (mod 4), eps(p) = eps(q), carries information on eps iff beta(c) * beta(-c) != 0, where beta(c) = E_X (c|X)
over the prime law (notes/leakmap_results.md, Theorem PT-J).  This module estimates z1 = sqrt(P) beta(c),
z2 = sqrt(P) beta(-c) over all primes X in [2^23, 2^24) for the pre-registered menu (notes/leakmap_prereg.md E3).

Usage: .venv/bin/python -m src.leakmap.jacobi [--bits 24] [--out experiments/raw/leakmap_e3.jsonl]
"""
from __future__ import annotations

import argparse
import json
import math
import random
from typing import Callable, Dict, List, Tuple

import numpy as np


def primes_in(lo: int, hi: int) -> np.ndarray:
    sieve = np.ones(hi, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(hi ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.flatnonzero(sieve[lo:]) + lo


def mulmod(a: np.ndarray, b: np.ndarray, X: np.ndarray) -> np.ndarray:
    return (a * b) % X          # X < 2^24 so products < 2^48 fit in int64


def powmod_scalar_exp(base: np.ndarray, e: int, X: np.ndarray) -> np.ndarray:
    r = np.ones_like(X)
    b = base % X
    while e:
        if e & 1:
            r = mulmod(r, b, X)
        b = mulmod(b, b, X)
        e >>= 1
    return r


def legendre_vec(c: np.ndarray, X: np.ndarray) -> np.ndarray:
    """(c|X) for odd primes X (vectorised Euler criterion with per-element exponent (X-1)/2)."""
    c = c % X
    e = (X - 1) // 2
    r = np.ones_like(X)
    b = c.copy()
    for _ in range(int(e.max()).bit_length()):
        bit = (e & 1).astype(bool)
        r = np.where(bit, mulmod(r, b, X), r)
        b = mulmod(b, b, X)
        e = e >> 1
    out = np.where(r == 1, 1, -1)
    return np.where(c == 0, 0, out)


def lucas_V(P: int, M: int, X: np.ndarray) -> np.ndarray:
    v0 = np.full_like(X, 2) % X
    v1 = np.full_like(X, P) % X
    for bit in bin(M)[2:]:
        if bit == "1":
            v0, v1 = (mulmod(v0, v1, X) - P) % X, (mulmod(v1, v1, X) - 2) % X
        else:
            v0, v1 = (mulmod(v0, v0, X) - 2) % X, (mulmod(v0, v1, X) - P) % X
    return v0


def lcm_upto(k: int) -> int:
    M = 1
    for i in range(2, k + 1):
        M = M * i // math.gcd(M, i)
    return M


def menu(X: np.ndarray) -> List[Tuple[str, Callable[[], np.ndarray], str]]:
    """(name, thunk -> c mod X array, structural note)."""
    items = []
    two = np.full_like(X, 2)
    three = np.full_like(X, 3)
    for k in range(5, 61):
        def f(k=k):
            r = two % X
            for _ in range(k):
                r = mulmod(r, r, X)
            return (r + 1) % X
        items.append((f"F1:2^(2^{k})+1", f, ""))
    for k in range(5, 31):
        items.append((f"F2:2^({k}!)-1", lambda k=k: (powmod_scalar_exp(two, math.factorial(k), X) - 1) % X, ""))
    for k in range(5, 61, 5):
        items.append((f"F3:3^lcm(1..{k})+2", lambda k=k: (powmod_scalar_exp(three, lcm_upto(k), X) + 2) % X, ""))
    for k in range(5, 31):
        items.append((f"F4:V_{k}!(3)-2", lambda k=k: (lucas_V(3, math.factorial(k), X) - 2) % X,
                      "= 5*U^2 (structural: symbol (5|X))"))
        items.append((f"F4:V_{k}!(3)+2", lambda k=k: (lucas_V(3, math.factorial(k), X) + 2) % X,
                      "EXCLUDED: = V_{k!/2}^2, a square"))
    for B in (10, 30, 100, 300, 1000):
        def f(B=B):
            r = np.ones_like(X)
            for x in range(1, B + 1):
                r = mulmod(r, np.full_like(X, x * x + 1) % X, X)
            return r
        items.append((f"F5:prod_(x<={B})(x^2+1)", f, ""))
    for B in (10, 100, 1000, 10000):
        def f(B=B, plus=0):
            r = np.ones_like(X)
            for x in range(2, B + 1):
                r = mulmod(r, np.full_like(X, x) % X, X)
            return (r + plus) % X
        items.append((f"F6:{B}!", f, ""))
        items.append((f"F6:{B}!+1", lambda B=B, f=f: f(B, 1), ""))
    rng = random.Random("leakmapE3")
    for j in range(300):
        prog = []
        for _ in range(40):
            op = rng.choice(("sq", "mul", "add", "inc"))
            prog.append((op, rng.random()))

        def f(prog=prog):
            vals = [np.full_like(X, 2), np.full_like(X, 3)]
            for op, u in prog:
                a = vals[-1]
                b = vals[int(u * len(vals))]
                if op == "sq":
                    v = mulmod(a, a, X)
                elif op == "mul":
                    v = mulmod(a, b, X)
                elif op == "add":
                    v = (a + b) % X
                else:
                    v = (a + 1) % X
                vals.append(v)
            return vals[-1]
        items.append((f"F7:rand{j}", f, ""))
    return items


def run(bits: int, out: str) -> Dict:
    X = primes_in(1 << (bits - 1), 1 << bits).astype(np.int64)
    P = len(X)
    eps = np.where(X % 4 == 1, 1, -1)
    rng = np.random.default_rng(12345)
    probe = rng.choice(P, 64, replace=False)
    recs = []
    hits = 0
    maxprod = 0.0
    with open(out, "w") as f:
        for name, thunk, note in menu(X):
            c = thunk()
            sym = legendre_vec(c, X)
            zeros = int((sym == 0).sum())
            b1 = float(sym.mean())
            b2 = float((sym * eps).mean())
            z1, z2 = b1 * math.sqrt(P), b2 * math.sqrt(P)
            # structural square check on 64 probe primes (for c, -c, 2c, -2c)
            sq = {}
            for lab, mult in (("c", 1), ("-c", -1), ("2c", 2), ("-2c", -2)):
                s = legendre_vec(mult * c[probe], X[probe])
                sq[lab] = bool(np.all(s >= 0))
            excluded = note.startswith("EXCLUDED") or any(sq.values())
            hit = (not excluded) and min(abs(z1), abs(z2)) >= 4.5
            hits += int(hit)
            if not excluded:
                maxprod = max(maxprod, abs(z1 * z2))
            rec = dict(name=name, primes=P, zeros=zeros, beta=b1, beta_minus=b2, z1=z1, z2=z2,
                       square_probe=sq, excluded=excluded, note=note, hit=hit)
            recs.append(rec)
            f.write(json.dumps(rec) + "\n")
    summary = dict(bits=bits, primes=P, eps_mean=float(eps.mean()), items=len(recs), hits=hits,
                   max_abs_z1z2_nonexcluded=maxprod,
                   max_min_abs_z=max(min(abs(r["z1"]), abs(r["z2"])) for r in recs if not r["excluded"]))
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, default=24)
    ap.add_argument("--out", default="experiments/raw/leakmap_e3.jsonl")
    a = ap.parse_args(argv)
    s = run(a.bits, a.out)
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
