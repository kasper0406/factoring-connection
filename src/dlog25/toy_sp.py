"""E1 (notes/dlog25_prereg.md): exhaustive check of Lemma S_p'' (split counting for T uniform on F_p^*).

Lemma S_p''.  p prime; G a set of primes < p; S the products of exactly k distinct primes of G lying in [Y/2, Y],
2Y^2 < p.  For T in F_p^* let E(T) be the event that the Gauss-reduced shortest nonzero vector of
Lambda_T = {(u, v) in Z^2 : u = T v (mod p)} is +-(u, v) with u, v in S and gcd(u, v) = 1.  Then
   #{T : E(T)} = #{(u, v) in S^2 : gcd(u, v) = 1} =: P'        (exact identity, proved)
and P' >= |S|^2 (1 - 2k^3/|G|) when |S| >= C(|G|, k)/(2k)    (so Pr[E] >= |S|^2/(2(p - 1)) if |G| >= 4k^3).

DIAGNOSTIC/HARNESS code (uses only public data; no secrets are involved in this lemma).
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time
from math import comb, gcd

from sympy import isprime, primerange


def gauss_reduce(p: int, T: int):
    """Shortest nonzero vector of {(u, v): u = T v mod p} (Lagrange-Gauss reduction, exact integers)."""
    b1, b2 = (p, 0), (T % p, 1)
    n1, n2 = b1[0] ** 2 + b1[1] ** 2, b2[0] ** 2 + b2[1] ** 2
    if n1 < n2:
        b1, b2, n1, n2 = b2, b1, n2, n1
    while True:
        # make b2 the shorter; reduce b1 by b2
        dot = b1[0] * b2[0] + b1[1] * b2[1]
        # nearest integer to dot/n2
        m = (2 * dot + n2) // (2 * n2)
        b1 = (b1[0] - m * b2[0], b1[1] - m * b2[1])
        n1 = b1[0] ** 2 + b1[1] ** 2
        if n1 >= n2:
            return b2
        b1, b2, n1, n2 = b2, b1, n2, n1


def build_S(G, k, Y):
    S = []
    for c in itertools.combinations(G, k):
        x = math.prod(c)
        if Y // 2 <= x <= Y:     # window [Y/2, Y] (integer endpoints, Y/2 rounded down: still 2Y^2 < p)
            S.append(x)
    return sorted(S)


def run_case(p: int, k: int, Glo: int, Ghi: int, Y: int):
    assert isprime(p) and 2 * Y * Y < p
    G = list(primerange(Glo + 1, Ghi + 1))
    S = build_S(G, k, Y)
    Sset = set(S)
    Pprime = sum(1 for u in S for v in S if gcd(u, v) == 1)
    hits = 0
    t0 = time.time()
    for T in range(1, p):
        u, v = gauss_reduce(p, T)
        if v < 0 or (v == 0 and u < 0):
            u, v = -u, -v
        if u in Sset and v in Sset and gcd(u, v) == 1:
            hits += 1
    lb_hyp = len(S) >= comb(len(G), k) / (2 * k)
    rec = {"p": p, "k": k, "G": [Glo, Ghi], "absG": len(G), "Y": Y, "absS": len(S), "Pprime": Pprime,
           "hits": hits, "identity_ok": hits == Pprime,
           "Sbig_hyp": lb_hyp, "G_ge_4k3": len(G) >= 4 * k ** 3,
           "Pprime_ge_bound": Pprime >= len(S) ** 2 * (1 - 2 * k ** 3 / len(G)) if G else True,
           "pr_ratio": (hits / (p - 1)) / (len(S) ** 2 / (2 * (p - 1))) if S else None,
           "secs": round(time.time() - t0, 2)}
    return rec


CASES = [
    # (p, k, G-interval (lo, hi], Y)
    (200003, 1, 128, 256, 256),
    (200003, 2, 8, 32, 300),
    (1000003, 1, 256, 512, 512),
    (1000003, 2, 16, 64, 700),
    (2000003, 3, 4, 16, 1000),
    (2000003, 2, 8, 64, 1000),
    (4000037, 3, 6, 24, 1400),
    (4000037, 2, 4, 128, 1400),
]


def main(out_path: str):
    with open(out_path, "w") as fh:
        for p, k, lo, hi, Y in CASES:
            rec = run_case(p, k, lo, hi, Y)
            fh.write(json.dumps(rec) + "\n")
            print(rec, flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
