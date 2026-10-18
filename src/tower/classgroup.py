"""Exact 2-Sylow structure of Cl(D), D < 0, D = 0, 1 (mod 4), primitive positive definite forms.  DIAGNOSTIC.

Method (no GRH, no analytic shortcuts):
  1. enumerate ALL reduced primitive forms (a, b, c) of discriminant D (numpy-vectorized over b for each a);
     h = their number (= |Cl(D)| for the form class group of the order of discriminant D);
  2. write h = 2^e m with m odd.  The map x -> x^m sends Cl(D) onto its 2-Sylow subgroup S (|S| = 2^e);
  3. draw uniformly random reduced forms f (seeded), adjoin f^m to a subgroup T <= S built by closure, until
     |T| = 2^e (then T = S; the stopping test is exact, so the result is exact whatever the random draws);
  4. count |S[2^j]| for all j and return the 2^j-ranks r_{2^j} = log2|S[2^j]| - log2|S[2^(j-1)]|.
Composition and reduction are copied from src/repair/verify_r4.py (imported, not edited); this module is validated
against that brute-force code in tests/test_tower.py.
"""
from __future__ import annotations

import random
from math import isqrt
from typing import Dict, List, Tuple

import numpy as np

from src.repair.verify_r4 import compose, identity_form, reduce_form

Form = Tuple[int, int, int]


def reduced_forms_np(D: int) -> np.ndarray:
    """All reduced primitive forms of discriminant D < 0 as an (h, 3) int64 array.

    Reduced: |b| <= a <= c, and b >= 0 if |b| = a or a = c.  Same set as verify_r4.reduced_forms."""
    assert D < 0 and D % 4 in (0, 1)
    assert -D < 2 ** 40
    out = []
    amax = isqrt(-D // 3)
    par = D % 2
    for a in range(1, amax + 1):
        # b in (-a, a], b = D (mod 2)
        lo = -a + 1
        if (lo - par) % 2:
            lo += 1
        b = np.arange(lo, a + 1, 2, dtype=np.int64)
        t = b * b - D
        m = t % (4 * a) == 0
        if not m.any():
            continue
        b = b[m]
        c = t[m] // (4 * a)
        keep = c >= a
        keep &= ~((c == a) & (b < 0))
        b, c = b[keep], c[keep]
        if b.size == 0:
            continue
        g = np.gcd(np.gcd(np.int64(a), np.abs(b)), c)
        keep = g == 1
        b, c = b[keep], c[keep]
        if b.size:
            out.append(np.stack([np.full(b.size, a, dtype=np.int64), b, c], axis=1))
    return np.concatenate(out) if out else np.zeros((0, 3), dtype=np.int64)


def form_pow(f: Form, n: int, one: Form) -> Form:
    r, x = one, f
    while n:
        if n & 1:
            r = compose(r, x)
        n >>= 1
        if n:
            x = compose(x, x)
    return r


def two_sylow(D: int, seed: int = 0, max_draws: int = 10000, return_elems: bool = False) -> Dict[str, object]:
    """Exact 2-Sylow data of Cl(D): h, e = v2(h), counts |S[2^j]|, ranks r_{2^j} (j = 1..), abelian type.
    With return_elems, also the 2-Sylow subgroup itself as a set of reduced forms (key "elems")."""
    forms = reduced_forms_np(D)
    h = int(forms.shape[0])
    e = (h & -h).bit_length() - 1
    m = h >> e
    one = identity_form(D)
    T = {one}
    rng = random.Random(f"tower128|{D}|{seed}")
    draws = 0
    while len(T) < (1 << e):
        draws += 1
        if draws > max_draws:
            raise RuntimeError(f"2-Sylow closure did not finish for D={D}")
        a, b, c = (int(x) for x in forms[rng.randrange(h)])
        g = form_pow((a, b, c), m, one)
        if g in T:
            continue
        # closure: T <- T * <g>, with g^k the first power landing in T
        powers = [one]
        x = g
        while x not in T:
            powers.append(x)
            x = compose(x, g)
        T = {compose(t, y) for t in T for y in powers}
    # orders in S
    counts = []  # counts[j] = |S[2^j]|, j = 0..e
    elems = list(T)
    order_exp = {}
    for x in elems:
        k, y = 0, x
        while y != one:
            y = compose(y, y)
            k += 1
        order_exp[x] = k
    for j in range(0, e + 1):
        counts.append(sum(1 for x in elems if order_exp[x] <= j))
    logs = [c.bit_length() - 1 for c in counts]
    assert all(c == 1 << l for c, l in zip(counts, logs))
    ranks = [logs[j] - logs[j - 1] for j in range(1, e + 1)]  # ranks[0] = r2, ranks[1] = r4, ...
    # abelian type: number of cyclic factors of order exactly 2^j is r_{2^j} - r_{2^(j+1)}
    typ = []
    for j in range(len(ranks)):
        nxt = ranks[j + 1] if j + 1 < len(ranks) else 0
        typ += [1 << (j + 1)] * (ranks[j] - nxt)
    res = {"D": D, "h": h, "e": e, "counts": counts, "ranks": ranks, "type": sorted(typ, reverse=True),
           "draws": draws}
    if return_elems:
        res["elems"] = T
    return res


def rank(info: Dict[str, object], j: int) -> int:
    """r_{2^j} (j = 1: 2-rank, 2: 4-rank, 3: 8-rank, 4: 16-rank)."""
    rk = info["ranks"]
    return rk[j - 1] if j - 1 < len(rk) else 0
