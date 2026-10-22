"""OPERATIONAL reductions for notes/selfred_results.md §2 (status of the open bits).

Each function receives N (or an integer X computed from N) and, where stated, an oracle passed in as an
argument; oracles live in harness.py (HARNESS).  Nothing here factors an integer: only gcd, exact division,
integer square roots, CRT and a probable-prime test (gmpy2.is_prime) are used.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Callable, Dict, List, Optional, Tuple

import gmpy2


# ----------------------------------------------------------------------------- Prop RAD
def squarefree_decomposition_via_rad(N: int, rad: Callable[[int], int]) -> Tuple[Dict[int, int], int]:
    """Return ({k: P_k}, calls) with P_k = product of the primes of exponent exactly k in N.

    r_k := rad(N / (r_1 ... r_{k-1})) = prod_{e_p >= k} p, so P_k = r_k / r_{k+1}.  Calls <= max_p e_p + 1.
    """
    out: Dict[int, int] = {}
    M, calls, rs = N, 0, []
    while M > 1:
        r = rad(M)
        calls += 1
        rs.append(r)
        M //= r
    rs.append(1)
    for k in range(1, len(rs)):
        Pk = rs[k - 1] // rs[k]
        if Pk > 1:
            out[k] = Pk
    return out, calls


# ----------------------------------------------------------------------------- Lemma SIG
def valuation_split(X: int, N: int) -> Tuple[str, int, int]:
    """For N = pq (p != q) and any nonzero integer X:

    G = gcd(X, N^inf) = p^a q^b.  If a != b this returns ('split', d, 0) with 1 < d < N, d | N.
    Otherwise it returns ('balanced', a, Y) with X = N^a * Y and gcd(Y, N) = 1.
    """
    X = abs(X)
    a = 0
    while True:
        g = gcd(X, N)
        if g == 1:
            return ("balanced", a, X)
        if g != N:
            return ("split", g, 0)
        X //= N
        a += 1


# ----------------------------------------------------------------------------- phi(N) mod l, l an input
def _is_small_prime(l: int) -> bool:
    """Deterministic trial division (l is O(n))."""
    if l < 2:
        return False
    d = 2
    while d * d <= l:
        if l % d == 0:
            return False
        d += 1
    return True


def small_primes_product_exceeding(bound: int, N: int) -> List[int]:
    ls, M, l = [], 1, 3
    while M <= bound:
        if _is_small_prime(l) and N % l:
            ls.append(l)
            M *= l
        l += 2
    return ls


def factor_via_phi_mod_ell(N: int, phi_mod: Callable[[int, int], int],
                           balanced: bool = True) -> Tuple[Optional[int], int]:
    """phi(N) mod l for input primes l  ==>  s = p + q exactly (CRT)  ==>  p.

    PROMISE: balanced N = pq (p, q < 2^ceil(n/2)), so s < 2^(ceil(n/2)+1); pass balanced=False to use the general
    bound s <= N + 1.  Calls = number of odd primes up to B with theta(B) ~ (n/2) ln 2, i.e. ~ 0.347 n / ln n
    asymptotically (finite n: 12, 26, 75, 131 calls at n = 80, 256, 1024, 2048).  Returns (factor, calls).
    """
    bound = (1 << ((N.bit_length() + 1) // 2 + 2)) if balanced else N + 2
    ls = small_primes_product_exceeding(bound, N)
    s, M = 0, 1
    for l in ls:
        r = (N + 1 - phi_mod(N, l)) % l          # s mod l
        # CRT: s = s0 mod M, s = r mod l
        t = ((r - s) * pow(M, -1, l)) % l
        s, M = s + M * t, M * l
    disc = s * s - 4 * N
    if disc < 0:
        return None, len(ls)
    d = isqrt(disc)
    if d * d != disc:
        return None, len(ls)
    p = (s - d) // 2
    return (p if 1 < p < N and N % p == 0 else None), len(ls)


# ----------------------------------------------------------------------------- Euclid over Z/N
def poly_gcd_mod_N(f: List[int], g: List[int], N: int) -> Tuple[str, object]:
    """Monic Euclid in (Z/N)[x].  Returns ('split', d) when a leading coefficient is a zero divisor, else
    ('gcd', h).  Used for the discriminant row: if p | disc f but q does not, deg gcd(f, f') differs mod p and
    mod q, so this must return 'split' (Lemma EUC)."""
    def tr(a):
        a = [c % N for c in a]
        while a and a[-1] == 0:
            a.pop()
        return a
    a, b = tr(f), tr(g)
    while b:
        lc = b[-1]
        d = gcd(lc, N)
        if d != 1:
            return ("split", d)
        inv = pow(lc, -1, N)
        b = [(c * inv) % N for c in b]
        while len(a) >= len(b):
            c = a[-1]
            k = len(a) - len(b)
            for i, bi in enumerate(b):
                a[i + k] = (a[i + k] - c * bi) % N
            a = tr(a)
        a, b = b, a
    return ("gcd", a)
