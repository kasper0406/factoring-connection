"""Square-coset re-split obstruction (Prop RS(ii) of notes/escapeR2_results.md).  HARNESS: uses the factors of a toy
modulus to build an adversarial target set.  Not operational; never imported by operational code.

A re-split of a unit x modulo N draws s, puts R = x s^2 mod N and Gauss-reduces {(u, v): u = R v mod N}.  Lemma S''
lower-bounds the success probability for x = 1 (R uniform on the squares) for ANY target set; for x outside the
square class of the target products the pair count sum_c n_c n_{cx} can vanish.
"""
from __future__ import annotations

import random
from itertools import combinations
from math import gcd
from typing import Dict, List, Tuple

from src.common.instances import balanced_semiprime
from src.thmA import arith


def sq_class(n: int, primes: List[int]) -> Tuple[int, ...]:
    """Square class of a unit n modulo prod(primes): the vector of Legendre symbols (harness)."""
    return tuple(arith.legendre(n, p) for p in primes)


def mulc(a: Tuple[int, ...], b: Tuple[int, ...]) -> Tuple[int, ...]:
    return tuple(x * y for x, y in zip(a, b))


def build(seed: int, bits: int = 48, k: int = 2) -> Dict:
    inst = balanced_semiprime(bits, seed)
    N, ps = inst.N, inst.labels["distinct"]
    n = N.bit_length()
    Y = 2 ** (-(-n // (2 * k)) - 1)
    window = [p for p in arith.primes_upto(2 * Y) if p > Y and N % p]
    Yp = 2 ** ((n - 2) // 2)                 # 2 Yp^2 < N
    good = [p for p in window if sq_class(p, ps) == (1,) * len(ps)]   # adversarial: all squares mod every prime of N
    return {"N": N, "ps": ps, "k": k, "Y": Y, "Yp": Yp, "window": window, "good": good}


def products(G: List[int], k: int, lo: int, hi: int) -> List[int]:
    out = []
    for c in combinations(G, k):
        m = 1
        for x in c:
            m *= x
        if lo <= m <= hi:
            out.append(m)
    return out


def pair_count(S: List[int], ps: List[int], xcls: Tuple[int, ...]) -> int:
    """#{(u, v) in S^2 : gcd(u, v) = 1, class(u) = class(x) class(v)}  (exact enumeration over class buckets)."""
    buckets: Dict[Tuple[int, ...], List[int]] = {}
    for s in S:
        buckets.setdefault(sq_class(s, ps), []).append(s)
    total = 0
    for c, us in buckets.items():
        vs = buckets.get(mulc(c, xcls), [])  # class(v) = class(u) * class(x)^{-1} = class(u) * class(x)
        for u in us:
            for v in vs:
                if gcd(u, v) == 1:
                    total += 1
    return total


def hits(N: int, x: int, Sset: set, draws: int, rng: random.Random) -> int:
    h = 0
    for _ in range(draws):
        s = rng.randrange(2, N - 1)
        if gcd(s, N) != 1:
            continue
        R = x * s * s % N
        (u, v), _ = arith.gauss_reduce((N, 0), (R, 1))
        if abs(u) in Sset and abs(v) in Sset:
            h += 1
    return h


def nonsquare_unit(N: int, ps: List[int], rng: random.Random) -> int:
    while True:
        x = rng.randrange(2, N - 1)
        if gcd(x, N) == 1 and arith.legendre(x, ps[0]) == -1:
            return x


def run_seed(seed: int, draws: int = 2000) -> Dict:
    d = build(seed)
    N, ps, k, Yp = d["N"], d["ps"], d["k"], d["Yp"]
    rng = random.Random(1000 + seed)
    S_adv = products(d["good"], k, Yp // 2, Yp)
    x = nonsquare_unit(N, ps, rng)
    xc = sq_class(x, ps)
    one = (1,) * len(ps)
    h = 2 ** len(ps)
    S_all = products(d["window"], k, Yp // 2, Yp)
    pc_one_all = pair_count(S_all, ps, one)
    return {
        "seed": seed, "N_bits": N.bit_length(), "k": k, "Y": d["Y"], "Yp": Yp,
        "window_primes": len(d["window"]), "adv_primes": len(d["good"]), "S_adv": len(S_adv),
        "x_class": list(xc),
        "R1_pairs_nonsquare_x": pair_count(S_adv, ps, xc),
        "R1_hits_nonsquare_x": hits(N, x, set(S_adv), draws, rng),
        "R1_pairs_x_eq_1_adv": pair_count(S_adv, ps, one),
        "S_all": len(S_all), "h": h,
        "R2_pairs_x_eq_1_all": pc_one_all,
        "R2_bound": len(S_all) ** 2 / (2 * h),
        "R2_hypothesis_G_ge_4k3h": len(d["window"]) >= 4 * k ** 3 * h,
    }
