"""HARNESS / DIAGNOSTIC code for `unbalanced`: uses p, q and group orders.
Never imported by algorithm.py / params.py.

* sample_pair_ratio(n, j, seed): uniform sample from
      S_{n,j} = {(p, q) primes : 2^j < q/p <= 2^{j+1}, 2^{n-1} <= pq < 2^n}
  (rejection from a box of integer pairs, so uniform over prime pairs).
* covering(...): instances / band primes with l | Bp - Aq AND p in the
  instance's ratio range (the theorem's coverage, Lemma 2 window included).
* predict(...): clean blocks, covering targets, FULL (closed form, Lemma 7).
* lemma_checks(...): Lemmas 1-7 on every covering target of a recorded block.
* classify_hit(...): STRUCTURAL_P / STRUCTURAL_Q / ACCIDENTAL.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Dict, List, Tuple

import gmpy2

from src.dcltheorem.harness import full_closed_form, ord2

from .algorithm import ratio_range
from .params import WideDesign


def is_prime(x: int) -> bool:
    return bool(gmpy2.is_prime(x, 50))


def pair_box(n: int, j: int) -> Tuple[int, int, int, int]:
    """p in [plo, phi], q in [qlo, qhi] contains S_{n,j}."""
    plo = isqrt((1 << (n - 1)) >> (j + 1)) - 1 if n - 1 >= j + 1 else 2
    phi = isqrt((1 << n) >> j) + 1
    qlo = isqrt((1 << (n - 1)) << j)
    qhi = isqrt((1 << n) << (j + 1)) + 1
    return max(2, plo), phi, qlo, qhi


def sample_pair_ratio(n: int, j: int, seed: int) -> Tuple[int, int, int]:
    rng = random.Random(f"unbalanced-pair|n={n}|j={j}|{seed}")
    plo, phi, qlo, qhi = pair_box(n, j)
    lo, hi = 1 << (n - 1), 1 << n
    while True:
        p = rng.randint(plo, phi)
        q = rng.randint(qlo, qhi)
        N = p * q
        if not (lo <= N < hi):
            continue
        if not ((q << 0) > (p << j) and q <= (p << (j + 1))):
            continue
        if is_prime(p) and is_prime(q):
            return N, p, q


def in_range(N: int, p: int, t) -> bool:
    x_lo, x_hi = ratio_range(N, t)
    return x_lo <= p <= x_hi


def clean_blocks(q: int, design: WideDesign) -> List[bool]:
    return [all((q - 1) % l for l in band) for band in design.bands]


def covering(N: int, p: int, q: int, design: WideDesign, k: int) -> List[Tuple]:
    out = []
    insts = [t for t in design.insts if in_range(N, p, t)]
    for l in design.bands[k]:
        for t in insts:
            if (t.B * p - t.A * q) % l == 0:
                out.append((t.A, t.B, l, t))
    return out


def predict(N: int, p: int, q: int, design: WideDesign, ms: List[int], d: int = None) -> Dict:
    d = d or ord2(q)
    cl = clean_blocks(q, design)
    blocks = []
    for k in range(design.K):
        cov = covering(N, p, q, design, k)
        fulls = [full_closed_form(p, q, A, B, l, ms[k], d) for (A, B, l, _) in cov] if ms[k] else []
        blocks.append({"clean": cl[k], "n_cover": len(cov), "n_full": sum(fulls),
                       "cert": cl[k] and any(not f for f in fulls)})
    return {"ord_q2": d, "blocks": blocks, "certified": any(b["cert"] for b in blocks)}


def classify_hit(N: int, p: int, q: int, blk: Dict) -> str:
    hit = blk.get("hit")
    if not hit:
        return "GCD_" + blk.get("status", "?")
    t = blk["targets"][hit["target"]]
    A, B, l, c, g_lo = t["A"], t["B"], t["l"], t["c"], t["g_lo"]
    m, j = blk["m"], hit["j"]
    lo, hi = g_lo + j * m, g_lo + j * m + m - 1
    f = blk["factor"]
    l2 = l * l
    if (B * p + A * q - c) % l2 == 0 and lo <= (B * p + A * q - c) // l2 <= hi and f == p:
        return "STRUCTURAL_P"
    if (B * q + A * p - c) % l2 == 0 and lo <= (B * q + A * p - c) // l2 <= hi and f == q:
        return "STRUCTURAL_Q"
    return "ACCIDENTAL"


def lemma_checks(N: int, p: int, q: int, design: WideDesign, blk: Dict, d: int = None) -> Dict:
    """On every covering target (p == r mod l AND p in the instance's ratio range):
    L1 integrality of g; L2 window contains g; L3 giant == gamma^{i*} (mod p);
    L4 giant == gamma^{i*} u^{(A-B)(p-1)} (mod q); L57 (clean block) gcd == p <=> not FULL."""
    k, m = blk["k"], blk["m"]
    gamma, u = blk["gamma"], blk["u"]
    d = d or ord2(q)
    clean = all((q - 1) % l for l in design.bands[k])
    bad = {"L1": 0, "L2": 0, "L3": 0, "L4": 0, "L57": 0}
    n_cov = n_full = 0
    gidx = {tj: h for h, tj in enumerate(blk["giant_meta"])}
    ginv = pow(gamma, -1, N)
    by_band = {}
    for ti, t in enumerate(blk["targets"]):
        A, B, l, r, c, g_lo, W = t["A"], t["B"], t["l"], t["r"], t["c"], t["g_lo"], t["W"]
        if p % l != r:
            continue
        a1, b1, a2, b2 = t["rho"]
        # ratio range membership: a1/b1 < N/p^2 <= a2/b2
        if not (a1 * p * p < N * b1 and a2 * p * p >= N * b2):
            continue
        n_cov += 1
        by_band[t["j"]] = by_band.get(t["j"], 0) + 1
        num = B * p + A * q - c
        if num % (l * l):
            bad["L1"] += 1
            continue
        g = num // (l * l)
        if not (g_lo <= g < g_lo + W):
            bad["L2"] += 1
            continue
        jj, istar = (g - g_lo) // m, (g - g_lo) % m
        h = gidx[(ti, jj)]
        G = pow(u[l], B + A * N - c, N) * pow(ginv, g_lo + jj * m, N) % N
        if G % p != pow(gamma, istar, p):
            bad["L3"] += 1
        if G % q != pow(gamma, istar, q) * pow(u[l], (A - B) * (p - 1), q) % q:
            bad["L4"] += 1
        full_cf = full_closed_form(p, q, A, B, l, m, d, istar)
        n_full += full_closed_form(p, q, A, B, l, m, d)
        if clean and (blk["giant_gcds"][h] == p) != (not full_cf):
            bad["L57"] += 1
    return {"k": k, "clean": clean, "n_cover": n_cov, "n_full": n_full, "bad": bad,
            "cover_by_band": by_band}
