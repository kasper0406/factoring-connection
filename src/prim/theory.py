"""Numerical companions of the theorems in notes/prim_results.md.  DIAGNOSTIC (uses the primes).

* d2_bound      : Theorem D2 (deterministic two-layer nodes, one additive shift) zero-set bound at a prime p.
* theorem_r_*   : exact evaluation of Theorem R (fresh-coin decoupling) on small fresh-base PRIM programs.
* collision_*   : the chain collision bound Pr_x[chain = w] <= prod_i gcd(e_i, p-1) / (p-1) (Lemma CB).
"""
from __future__ import annotations

from collections import defaultdict
from itertools import product
from math import gcd, isqrt, sqrt
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
import sympy

from .prog import vpow


def d2_bound(p: int, y: int, a: int, alpha: int) -> Tuple[float, int]:
    """Theorem D2: |Z| <= min_{d | p-1} gcd(alpha,d)*tau(d)*(n/d + d*(n/o)*(sqrt(p)+2)) + n/o,  o = ord_p(y^a).
    Returns (bound, best d)."""
    n = p - 1
    o = int(sympy.n_order(pow(y, a, p), p))
    best, bd = float("inf"), 1
    for d in sympy.divisors(n):
        val = gcd(alpha, d) * sympy.divisor_count(d) * (n / d + d * (n / o) * (sqrt(p) + 2)) + n / o
        if val < best:
            best, bd = val, d
    return float(best), int(bd)


# ----------------------------------------------------------------------------- Theorem R, exact
def chain_value(x: np.ndarray, p: int, exps: Sequence[int], shifts: Sequence[int]) -> np.ndarray:
    """((x^{e1} + c1)^{e2} + c2)...^{e_d} mod p, vectorised over x (exps already reduced mod p-1, >= 1)."""
    v = vpow(x % p, np.full(x.shape, exps[0]), p)
    for e, c in zip(exps[1:], shifts):
        v = vpow((v + c) % p, np.full(x.shape, e), p)
    return v


def local_law(p: int, N: int, probes) -> Dict[tuple, float]:
    """Exact law of the (non-adaptive) transcript at p over independent uniform fresh bases x_v in F_p
    (one coin per probe).  probes: list of (exp_fns, shifts, w) with exp_fns callables N -> exponent."""
    xs = np.arange(p, dtype=np.int64)
    marg = []
    for exp_fns, shifts, w in probes:
        exps = [max(1, f(N)) % (p - 1) or (p - 1) for f in exp_fns]
        z = chain_value(xs, p, exps, shifts) == (w % p)
        marg.append(float(z.mean()))
    law = {}
    for bits in product((0, 1), repeat=len(probes)):
        pr = 1.0
        for b_, m in zip(bits, marg):
            pr *= m if b_ else 1 - m
        law[bits] = pr
    return law


def theorem_r_check(primes: Sequence[int], probes, target: Callable[[int, int], int]):
    """Exact Pr[correct or split] (Bayes labels), s, and the Theorem R right side B + 4s + 2Pr[D_modal & m_p = m_q]
    with kappa0 = all-zero transcript, over uniform distinct ordered pairs."""
    pairs = [(p, q) for p in primes for q in primes if p != q]
    w = 1.0 / len(pairs)
    s = 0.0
    dm = 0.0
    joint = defaultdict(lambda: [0.0, 0.0])   # transcript -> mass by target value (non-split)
    Fmass = [0.0, 0.0]
    ref = None
    for p, q in pairs:
        N = p * q
        Lp, Lq = local_law(p, N, probes), local_law(q, N, probes)
        if ref is None:
            ref = tuple(0 for _ in probes)
        F = target(p, q)
        Fmass[F] += w
        same = 0.0
        for t, a in Lp.items():
            b_ = Lq.get(t, 0.0)
            same += a * b_
            joint[t][F] += w * a * b_
        s += w * (1 - same)
        mp = max(Lp, key=Lp.get)
        mq = max(Lq, key=Lq.get)
        if mp == mq and (mp != ref or mq != ref):
            dm += w
    correct_ns = sum(max(v) for v in joint.values())
    lhs = s + correct_ns
    B = max(Fmass)
    return {"lhs": lhs, "B": B, "s": s, "Pr_Dmodal": dm, "rhs": B + 4 * s + 2 * dm, "holds": lhs <= B + 4 * s + 2 * dm + 1e-12}


def collision_check(p: int, exps: Sequence[int], shifts: Sequence[int], w: int) -> Tuple[float, float]:
    """(exact Pr_x[chain = w], collision bound prod gcd(e_i, p-1)/(p-1))."""
    xs = np.arange(1, p, dtype=np.int64)
    pr = float((chain_value(xs, p, exps, shifts) == w % p).mean())
    G = 1
    for e in exps:
        G *= gcd(e, p - 1)
    return pr, G / (p - 1)
