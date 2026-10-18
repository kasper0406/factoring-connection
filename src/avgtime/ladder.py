"""OPERATIONAL: the complete ladder algorithm LAD_n (notes/avgtime_results.md, section 2).

Sees only N and n-only public parameters.  Stages:
  0. small-order detector with T = 2^{floor(3n/10)}            (detector.py)
  1..k. DCL-det_{n,H_i} for the n-only ladder H_1 < H_2 < ... (src/dcltheorem)
  k+1. complete deterministic fallback (Harvey--Hittmeir 2022 in the theorem).
Whenever a proper divisor d is found, the output is the complete factorisation:
{d, N/d} if both are prime, else the fallback's factorisation of N.

In the theorem the fallback is HH 2022 and primality is AKS (deterministic,
polynomial).  In this reference implementation (small sizes, correctness
tests only) the fallback is Lehman's deterministic N^{1/3} method and
primality uses gmpy2.is_prime (BPSW-type; disclosed per AGENT_RULES rule 2:
used only to certify the final output, never inside a search).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Callable, Dict, List, Optional, Sequence

import gmpy2

from src.common.metrics import Metrics
from src.dcltheorem.algorithm import dcl_det
from src.dcltheorem.params import Design, icbrt, pow2_frac

from .detector import small_order_detector


def is_prime(x: int) -> bool:
    return x >= 2 and bool(gmpy2.is_prime(x, 50))


def _iroot(x: int, k: int) -> int:
    lo, hi = 0, 1 << (x.bit_length() // k + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= x:
            lo = mid
        else:
            hi = mid - 1
    return lo


def lehman_split(N: int) -> Optional[int]:
    """Lehman 1974: a proper divisor of N, or None iff N is 1 or prime. Deterministic."""
    if N < 4:
        return None
    c = _iroot(N, 3)
    for a in range(2, c + 2):
        if a < N and N % a == 0:
            return a
    for k in range(1, c + 2):
        four_kN = 4 * k * N
        r = isqrt(four_kN)
        xmax = _iroot(N // (4096 * k ** 3), 6) if N >= 4096 * k ** 3 else 0
        for a in range(r, r + xmax + 2):
            b2 = a * a - four_kN
            if b2 < 0:
                continue
            b = isqrt(b2)
            if b * b == b2:
                d = gcd(a + b, N)
                if 1 < d < N:
                    return d
    return None


def complete_factorisation(N: int, split: Callable[[int], Optional[int]] = lehman_split) -> List[int]:
    if N < 2:
        return []
    if is_prime(N):
        return [N]
    d = split(N)
    if d is None:            # cannot happen for composite N with a correct split
        raise RuntimeError("fallback failed on composite N")
    return sorted(complete_factorisation(d, split) + complete_factorisation(N // d, split))


def stage_design(n: int, H: int) -> Design:
    """DCL-det_{n,H}: Y = ceil(2^{n/6} / floor(H^{1/3})) (D' convention), K = 4, l = 3 mod 4.
    (At toy sizes Y > H^2 is enforced by max; the theorem has it for n >= n_0.)"""
    Y = -(-pow2_frac(n, 1, 6) // max(1, icbrt(H)))
    return Design(n=n, H=H, Y=max(3, H * H + 1, Y), K=4, mod4=True, label=f"ladder H={H}")


def ladder_factor(N: int, n: int, Hs: Sequence[int], T: Optional[int] = None,
                  split: Callable[[int], Optional[int]] = lehman_split,
                  met: Optional[Metrics] = None) -> Dict:
    """Run LAD on N.  Returns {'factors': sorted prime list, 'stage': name}."""
    met = met or Metrics()
    if N < 2:
        return {"factors": [], "stage": "trivial"}
    T = (1 << (3 * n // 10)) if T is None else T

    def finish(d: int, stage: str) -> Dict:
        e = N // d
        if is_prime(d) and is_prime(e):
            return {"factors": sorted([d, e]), "stage": stage}
        return {"factors": complete_factorisation(N, split), "stage": stage + "+fallback"}

    d = small_order_detector(N, T, met)
    if d is not None:
        return finish(d, "detector")
    for i, H in enumerate(Hs, start=1):
        r = dcl_det(N, stage_design(n, H), met)
        if r["factor"] is not None:
            return finish(r["factor"], f"dcl{i}")
    return {"factors": complete_factorisation(N, split), "stage": "fallback"}
