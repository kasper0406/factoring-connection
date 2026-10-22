"""Route 3: multi-sample (Cohn–Heninger-type) joint certificates over several windows.

Operational pieces see only window positions, h, P and y. Nothing uses a factor of N.

- `inheritance_scan`: exhaustive check of Lemma IH (same-N windows): if d > K*h divides an element of
  [k_i N, k_i N + h] and of [k_j N, k_j N + h] (k_i, k_j <= K), then the offsets are proportional,
  k_j r_i = k_i r_j.
- `joint_tuples_indep`: Monte Carlo of certified tuples (d, z_1..z_m) with d | P, ln d >= sigma, d | z_i,
  z_i in independently placed windows [x_i, x_i + h].
- `predicted_tuples`: E[#tuples] = sum_{d | P, ln d >= sigma} prod_i E[#multiples of d in W_i]
  = sum_d ((h+1)/d)^m for uniformly placed windows (exact up to boundary terms O(d/X)).
"""
from __future__ import annotations

import math
from typing import Dict, List, Sequence

import numpy as np

from .core import divisors_from_factored, is_smooth_int, primes_upto, smooth_mask_values


def inheritance_scan(N: int, ks: Sequence[int], h: int) -> Dict[str, int]:
    """Exhaustive over all pairs k_i < k_j in ks and all offsets r_i, r_j in [0, h].

    Returns counts: pairs with gcd > K*h ('big'), of which proportional ('inherited') and violations.
    Also counts pairs with h < gcd <= K*h ('mid') and how many of those are proportional.
    """
    K = max(ks)
    if K * N + h >= (1 << 62):
        raise ValueError("too large for int64")
    r = np.arange(h + 1, dtype=np.int64)
    out = {"big": 0, "big_inherited": 0, "violations": 0, "mid": 0, "mid_inherited": 0}
    for a in range(len(ks)):
        for b in range(a + 1, len(ks)):
            ki, kj = ks[a], ks[b]
            zi = ki * N + r
            zj = kj * N + r
            g = np.gcd.outer(zi, zj)
            prop = np.equal.outer(kj * r, ki * r)
            big = g > K * h
            mid = (g > h) & ~big
            out["big"] += int(big.sum())
            out["big_inherited"] += int((big & prop).sum())
            out["violations"] += int((big & ~prop).sum())
            out["mid"] += int(mid.sum())
            out["mid_inherited"] += int((mid & prop).sum())
    return out


def _divisors_of(g: int, primes: Sequence[int]) -> List[int]:
    f = {}
    for p in primes:
        if g % p == 0:
            e = 0
            while g % p == 0:
                g //= p
                e += 1
            f[p] = e
        if g == 1:
            break
    assert g == 1, "g must divide P"
    return [int(d) for d in divisors_from_factored(f)]


def predicted_tuples(P_divs: np.ndarray, h: int, m: int, sigma: float) -> float:
    s = 0.0
    lim = math.exp(sigma)
    for d in P_divs:
        d = int(d)
        if d >= lim:
            s += ((h + 1) / d) ** m
    return s


def predicted_smooth_first(P_divs: np.ndarray, h: int, m: int, sigma: float, X: int, dens) -> float:
    """Model: E[#tuples with z_1 y-smooth] = sum_d ((h+1)/d)^m * dens(1.5 X / d) (cofactor of size ~1.5X/d
    smooth like a random integer of its size; d | P is smooth)."""
    s = 0.0
    lim = math.exp(sigma)
    for d in P_divs:
        d = int(d)
        if d >= lim:
            c = max(1.5 * X / d, 1.0)
            s += ((h + 1) / d) ** m * (1.0 if c <= dens.y else float(dens(c)))
    return s


def joint_tuples_indep(X: int, h: int, P: int, B: int, m: int, sigma: float, y: int, n_calls: int,
                       rng: np.random.Generator) -> Dict[str, float]:
    """Monte Carlo over independent uniform window positions x_i in [X, 2X).

    A certified tuple is (d, z_1, ..., z_m) with d | P, ln d >= sigma, d | z_i, z_i in [x_i, x_i + h].
    Counts tuples, tuples with z_1 y-smooth, tuples with all z_i y-smooth, and calls with >= 1 tuple.
    """
    primes = primes_upto(B)
    lim = math.exp(sigma)
    tot = smooth1 = smooth_all = calls_hit = 0
    smooth_z1_distinct = 0
    r = np.arange(h + 1, dtype=np.int64)
    for _ in range(n_calls):
        xs = [int(v) for v in rng.integers(X, 2 * X, size=m)]
        z1 = xs[0] + r
        g = np.gcd(z1, np.int64(P))
        idx = np.nonzero(g >= lim)[0]
        call_tuples = 0
        z1_smooth_seen = set()
        for i in idx:
            zz = int(z1[i])
            for d in _divisors_of(int(g[i]), primes):
                if d < lim:
                    continue
                # multiples of d in the other windows
                others: List[List[int]] = []
                ok = True
                for x in xs[1:]:
                    first = -(-x // d) * d
                    mult = list(range(first, x + h + 1, d))
                    if not mult:
                        ok = False
                        break
                    others.append(mult)
                if not ok:
                    continue
                cnt = 1
                for o in others:
                    cnt *= len(o)
                call_tuples += cnt
                s1 = is_smooth_int(zz, y)
                if s1:
                    smooth1 += cnt
                    z1_smooth_seen.add(zz)
                    # all coordinates smooth: count combinations
                    good = 1
                    for o in others:
                        good *= sum(1 for w in o if is_smooth_int(w, y))
                    smooth_all += good
        tot += call_tuples
        calls_hit += call_tuples > 0
        smooth_z1_distinct += len(z1_smooth_seen)
    return {"tuples": tot, "tuples_z1_smooth": smooth1, "tuples_all_smooth": smooth_all,
            "calls_with_tuple": calls_hit, "distinct_smooth_z1": smooth_z1_distinct, "calls": n_calls}


def divisor_trick_yield(X: int, h: int, y: int, n_calls: int, rng: np.random.Generator) -> Dict[str, float]:
    """Baseline: d = 2^floor(log2(h+1)), z = d*ceil(x/d) in [x, x+h]; success iff z y-smooth."""
    d = 1 << int(math.floor(math.log2(h + 1)))
    xs = rng.integers(X, 2 * X, size=n_calls)
    z = ((xs + d - 1) // d) * d
    ok = smooth_mask_values(z.astype(np.uint64), y)
    return {"d": d, "successes": int(ok.sum()), "calls": n_calls}
