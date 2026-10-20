"""DIAGNOSTIC checkers for the systeq barriers (they may use the factors of toy N; never
imported by operational code).  Exact arithmetic.

Prop F  (additive power maps): if z -> z^s is additive on a Z/N-algebra A, then (p-1) | (s-1)
        for every prime p dividing char A.  Checked by brute force on tiny algebras.
Lemma CF (constant field): a systematic product prod_{c in S}(T - c) vanishing on the image
        of Z in A has |S| >= p for each p | char A; over F_p[sqrt D] (D a non-residue) the
        Wilson-type identity prod_{c=0}^{p-1}(z - c) = sigma(z) - z holds.
Lemma SUB (substitution Frobenius in char 0): with sigma_q(f)(t) = f(t^q),
        sigma_q(f) g - f sigma_q(g) = (ad - bc)(t^q - t) for f = at+b, g = ct+d: the
        homography trick yields ONE relation, not q^3.
Zsigmondy: m^n - 1 (m >= 2, n >= B >= 7) is never B-smooth.
"""
from __future__ import annotations

import random
from itertools import product
from typing import List, Tuple

from sympy import Poly, symbols

from .intanalog import is_smooth


# --- Prop F ------------------------------------------------------------------------------
def _alg_mul(a, b, F, N):
    """Multiply in (Z/N)[t]/(F), F monic given low->high; elements are coefficient tuples."""
    d = len(F) - 1
    r = [0] * (2 * d - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    for k in range(len(r) - 1, d - 1, -1):
        c = r[k]
        if c:
            for i in range(d + 1):
                r[k - d + i] -= c * F[i]
    return tuple(v % N for v in r[:d])


def _alg_pow(a, s, F, N):
    d = len(F) - 1
    r = tuple([1] + [0] * (d - 1))
    while s:
        if s & 1:
            r = _alg_mul(r, a, F, N)
        a = _alg_mul(a, a, F, N)
        s >>= 1
    return r


def additive_exponents(N: int, F: List[int], s_max: int) -> List[int]:
    """All s in [2, s_max] with (x+y)^s = x^s + y^s for ALL x, y in A = (Z/N)[t]/(F)."""
    d = len(F) - 1
    elems = list(product(range(N), repeat=d))
    out = []
    for s in range(2, s_max + 1):
        pw = {e: _alg_pow(e, s, F, N) for e in elems}
        ok = True
        for x in elems:
            for y in elems:
                z = tuple((u + v) % N for u, v in zip(x, y))
                if pw[z] != tuple((u + v) % N for u, v in zip(pw[x], pw[y])):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            out.append(s)
    return out


# --- Lemma CF ------------------------------------------------------------------------------
def wilson_frobenius_check(p: int, D: int, trials: int = 50, seed: int = 0) -> bool:
    """In F_p[sqrt D], D a non-residue: prod_{c=0}^{p-1} (z - c) == sigma(z) - z."""
    assert pow(D % p, (p - 1) // 2, p) == p - 1
    rng = random.Random(f"systeq|wilson|{p}|{D}|{seed}")
    mul = lambda u, v: ((u[0] * v[0] + D * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p)
    for _ in range(trials):
        z = (rng.randrange(p), rng.randrange(p))
        pr = (1, 0)
        for c in range(p):
            pr = mul(pr, ((z[0] - c) % p, z[1]))
        sig = (z[0], (-z[1]) % p)
        if pr != ((sig[0] - z[0]) % p, (sig[1] - z[1]) % p):
            return False
    return True


def systematic_degree_lower_bound(S: List[int], p: int) -> bool:
    """prod_{c in S}(T - c) vanishes on F_p  =>  S covers all residues mod p (so |S| >= p)."""
    return {c % p for c in S} == set(range(p))


# --- Lemma SUB -----------------------------------------------------------------------------
def substitution_homography(q: int):
    a, b, c, d, t = symbols("a b c d t")
    f, g = a * t + b, c * t + d
    sf, sg = a * t ** q + b, c * t ** q + d
    lhs = Poly(sf * g - f * sg, t)
    rhs = Poly((a * d - b * c) * (t ** q - t), t)
    return lhs == rhs


# --- Zsigmondy ---------------------------------------------------------------------------
def zsigmondy_check(B: int, m_max: int, n_extra: int) -> Tuple[int, int]:
    """Count (m, n), 2 <= m <= m_max, B <= n < B + n_extra, with m^n - 1 B-smooth (expect 0)."""
    hits, tot = 0, 0
    for m in range(2, m_max + 1):
        for n in range(B, B + n_extra):
            tot += 1
            hits += is_smooth(m ** n - 1, B)
    return hits, tot
