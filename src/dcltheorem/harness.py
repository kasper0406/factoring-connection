"""HARNESS / DIAGNOSTIC code for dcltheorem: uses p, q and group orders.
Never imported by algorithm.py / params.py / stage3.py.

* sample_pair(n, seed): uniform sample from S_n = {(p, q) primes: p < q <= 2p,
  2^(n-1) <= pq < 2^n}  (rejection from a box, so uniform over pairs).
* predict(N, p, q, design): the theorem's success certificate per block --
  clean block (no band prime divides q-1), covering targets (l | Bp - Aq),
  and for each covering target whether its true giant is FULL mod q, both by
  the closed form of Lemma 7 (d | (A-B)(p-1) - l^2 delta, d = ord_q 2) and
  directly.
* classify(...): STRUCTURAL_P / STRUCTURAL_Q / ACCIDENTAL for an operational hit.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

import gmpy2
from sympy import n_order

from .params import Design


def is_prime(x: int) -> bool:
    return bool(gmpy2.is_prime(x, 50))


def sample_pair(n: int, seed: int) -> Tuple[int, int, int]:
    rng = random.Random(f"dcltheorem-pair|n={n}|{seed}")
    plo = isqrt(1 << (n - 2))
    phi = isqrt((1 << n) - 1) + 1
    qlo = isqrt(1 << (n - 1))
    qhi = isqrt(1 << (n + 1)) + 1
    while True:
        p = rng.randint(plo, phi)
        q = rng.randint(qlo, qhi)
        if not (p < q <= 2 * p):
            continue
        N = p * q
        if not ((1 << (n - 1)) <= N < (1 << n)):
            continue
        if is_prime(p) and is_prime(q):
            return N, p, q


def ord2(x: int) -> int:
    return int(n_order(2, x))


def clean_blocks(q: int, design: Design) -> List[bool]:
    return [all((q - 1) % l for l in band) for band in design.bands]


def covering(p: int, q: int, design: Design, k: int) -> List[Tuple[int, int, int]]:
    """(A, B, l) with l | Bp - Aq in block k (the theorem's coverage)."""
    out = []
    for l in design.bands[k]:
        for (A, B) in design.V:
            if (B * p - A * q) % l == 0:
                out.append((A, B, l))
    return out


def full_closed_form(p: int, q: int, A: int, B: int, l: int, m: int, d: int,
                     istar: int = None) -> bool:
    """Lemma 7: the true giant is FULL iff d | (A-B)(p-1) - l^2 delta for some
    delta = i' - i* with 0 <= i' < m.  With istar=None the SUPERSET |delta| < m
    is used (the event bounded in the proof of Lemma 10)."""
    X = (A - B) * (p - 1)
    l2 = l * l
    lo, hi = (-(m - 1), m - 1) if istar is None else (-istar, m - 1 - istar)
    g = gcd(l2, d)
    if X % g:
        return False
    dd = d // g
    delta0 = (X // g) * pow(l2 // g, -1, dd) % dd if dd > 1 else 0
    # is there delta == delta0 (mod dd) in [lo, hi]?
    first = lo + ((delta0 - lo) % dd)
    return first <= hi


def predict(N: int, p: int, q: int, design: Design, ms: List[int]) -> Dict:
    """Theorem certificate per block (needs the m actually used per block)."""
    d = ord2(q)
    cl = clean_blocks(q, design)
    blocks = []
    for k in range(design.K):
        cov = covering(p, q, design, k)
        fulls = [full_closed_form(p, q, A, B, l, ms[k], d) for (A, B, l) in cov] if ms[k] else []
        blocks.append({"clean": cl[k], "n_cover": len(cov),
                       "n_full": sum(fulls), "cert": cl[k] and any(not f for f in fulls)})
    return {"ord_q2": d, "blocks": blocks, "certified": any(b["cert"] for b in blocks)}


def classify_hit(N: int, p: int, q: int, block_res: Dict) -> str:
    """For the giant that produced the factor in a recorded block."""
    hit = block_res.get("hit")
    if not hit:
        return "GCD_" + block_res.get("status", "?")
    t = block_res["targets"][hit["target"]]
    A, B, l, c, g_lo, W = t["A"], t["B"], t["l"], t["c"], t["g_lo"], t["W"]
    m = block_res["m"]
    j = hit["j"]
    lo, hi = g_lo + j * m, g_lo + j * m + m - 1
    f = block_res["factor"]
    if (B * p + A * q - c) % (l * l) == 0 and lo <= (B * p + A * q - c) // (l * l) <= hi and f == p:
        return "STRUCTURAL_P"
    if (B * q + A * p - c) % (l * l) == 0 and lo <= (B * q + A * p - c) // (l * l) <= hi and f == q:
        return "STRUCTURAL_Q"
    return "ACCIDENTAL"


def lemma_checks(N: int, p: int, q: int, design: Design, block_res: Dict) -> Dict:
    """Verify Lemmas 1-7 on every covering target of a recorded block:
    L1 integrality of g, L3 window, L2 G == gamma^{i*} (mod p), L4 closed form mod q,
    L5/L7 one-sided gcd == p  <=>  not FULL (closed form), when the block is clean."""
    k = block_res["k"]
    m = block_res["m"]
    gamma, u = block_res["gamma"], block_res["u"]
    d = ord2(q)
    clean = all((q - 1) % l for l in design.bands[k])
    bad = {"L1": 0, "L3": 0, "L2": 0, "L4": 0, "L57": 0}
    n_cov = 0
    n_full = 0
    meta = block_res["giant_meta"]
    gidx = {}
    for h, (ti, j) in enumerate(meta):
        gidx[(ti, j)] = h
    for ti, t in enumerate(block_res["targets"]):
        A, B, l, r, c, g_lo, W = t["A"], t["B"], t["l"], t["r"], t["c"], t["g_lo"], t["W"]
        if p % l != r:
            continue
        n_cov += 1
        num = B * p + A * q - c
        if num % (l * l):
            bad["L1"] += 1
            continue
        g = num // (l * l)
        if not (g_lo <= g < g_lo + W):
            bad["L3"] += 1
            continue
        istar_full = g - g_lo
        j = istar_full // m
        istar = istar_full % m
        h = gidx[(ti, j)]
        # recompute the giant
        L_over = None
        G = pow(u[l], B + A * N - c, N) * pow(pow(gamma, -1, N), g_lo + j * m, N) % N
        if G % p != pow(gamma, istar, p):
            bad["L2"] += 1
        rhs = pow(gamma, istar, q) * pow(u[l], (A - B) * (p - 1), q) % q
        if G % q != rhs:
            bad["L4"] += 1
        full_cf = full_closed_form(p, q, A, B, l, m, d, istar)
        n_full += full_closed_form(p, q, A, B, l, m, d)
        gg = block_res["giant_gcds"][h]
        if clean:
            if (gg == p) != (not full_cf):
                bad["L57"] += 1
        else:
            # without cleanliness only "not full => gcd p" direction is claimed
            # via the general condition; record direct fullness instead
            pass
    return {"k": k, "clean": clean, "n_cover": n_cov, "n_full": n_full, "bad": bad}
