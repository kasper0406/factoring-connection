"""abcsunit core: exact enumeration of smooth numbers, S-unit congruences mod N, triply-smooth
relations (S'-unit equations), F2 rank / fruitfulness, and N-free smooth abc triples.

Labelling (AGENT_RULES 2):
* OPERATIONAL routines see only N and public parameters (y, X): smooth_numbers, congruent_buckets,
  relations_from_buckets, f2_dependencies, square_from_dependency.
* DIAGNOSTIC routines (suffix ``_diag``) use the factors p, q (subgroup order |H|) or sympy factorint
  (radicals of non-smooth multipliers k).  They are never imported by operational code.

Exact integer arithmetic throughout; numpy int64 is used only as an exact integer container
(all values are < 2**62 and every product is range-checked before it is formed).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

INT64_SAFE = 1 << 62


class TooMany(RuntimeError):
    """Raised when an enumeration would exceed its registered cap."""


# ----------------------------------------------------------------------------------------------
# primes and smooth numbers (operational)
# ----------------------------------------------------------------------------------------------

def primes_upto(y: int) -> List[int]:
    if y < 2:
        return []
    sieve = bytearray([1]) * (y + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(y) + 1):
        if sieve[i]:
            sieve[i * i:: i] = bytearray(len(sieve[i * i:: i]))
    return [i for i in range(y + 1) if sieve[i]]


def smooth_numbers(X: int, primes: Sequence[int], cap: int = 30_000_000) -> np.ndarray:
    """All positive integers <= X whose prime factors lie in `primes`, sorted (exact).

    Built prime by prime: vals <- vals ∪ vals*p ∪ vals*p^2 ∪ ... (each product formed only
    for vals <= X // p, so no int64 overflow).  Raises TooMany if the count exceeds `cap`.
    """
    if X >= INT64_SAFE:
        raise ValueError("X must be < 2**62")
    vals = np.array([1], dtype=np.int64)
    for p in sorted(primes):
        pieces = [vals]
        cur = vals
        while True:
            cur = cur[cur <= X // p] * p
            if cur.size == 0:
                break
            pieces.append(cur)
        vals = np.concatenate(pieces)
        if vals.size > cap:
            raise TooMany(f"Psi > cap={cap} after prime {p}")
    vals.sort()
    return vals


def smooth_count(X: int, primes: Sequence[int]) -> int:
    """Psi(X, primes) by exact recursion (no storage; for cap checks and tests)."""
    ps = sorted(primes, reverse=True)

    def rec(x: int, i: int) -> int:
        if i == len(ps):
            return 1
        p = ps[i]
        tot = 0
        while x >= 1:
            tot += rec(x, i + 1)
            x //= p
        return tot

    return rec(X, 0)


def exponent_vector(n: int, primes: Sequence[int]) -> Optional[List[int]]:
    """Exponents of n over `primes` (trial division), or None if n is not smooth."""
    if n <= 0:
        raise ValueError("n must be positive")
    e = []
    for p in primes:
        k = 0
        while n % p == 0:
            n //= p
            k += 1
        e.append(k)
    return e if n == 1 else None


def is_smooth(n: int, primes: Sequence[int]) -> bool:
    if n <= 0:
        return False
    for p in primes:
        while n % p == 0:
            n //= p
    return n == 1


# ----------------------------------------------------------------------------------------------
# S-unit congruences u ≡ v (mod N)  (operational)
# ----------------------------------------------------------------------------------------------

def residue_class_sizes(vals: np.ndarray, N: int) -> np.ndarray:
    """Sizes c_r of the nonempty residue classes r = u mod N among `vals` (exact)."""
    if N >= INT64_SAFE:
        raise ValueError("N too large for the int64 container")
    _, counts = np.unique(vals % N, return_counts=True)
    return counts


def congruent_pair_count(vals: np.ndarray, N: int) -> Tuple[int, int, int]:
    """(unordered pairs u < v with u ≡ v mod N, number of occupied classes, max class size)."""
    c = residue_class_sizes(vals, N).astype(object)
    pairs = int(sum(int(x) * (int(x) - 1) // 2 for x in c if x > 1))
    return pairs, int(len(c)), int(max(c)) if len(c) else 0


def congruent_buckets(vals: np.ndarray, N: int) -> List[List[int]]:
    """The residue classes mod N containing >= 2 smooth numbers, as sorted lists of Python ints."""
    r = vals % N
    order = np.lexsort((vals, r))
    rs = r[order]
    vs = vals[order]
    out: List[List[int]] = []
    if rs.size == 0:
        return out
    brk = np.flatnonzero(np.diff(rs)) + 1
    starts = np.concatenate(([0], brk))
    ends = np.concatenate((brk, [rs.size]))
    big = np.flatnonzero(ends - starts >= 2)
    for j in big:
        out.append([int(t) for t in vs[starts[j]:ends[j]]])
    return out


def relations_from_buckets(buckets: Iterable[List[int]], N: int, coprime_only: bool = True,
                           consecutive_only: bool = False) -> List[Tuple[int, int, int]]:
    """Doubly-smooth relations (u, v, k) with u < v, v - u = kN.

    consecutive_only=True keeps only neighbouring members of a class (the others are
    multiplicatively implied up to gcds); default lists all pairs.
    """
    rel = []
    for b in buckets:
        if consecutive_only:
            pairs = zip(b[:-1], b[1:])
        else:
            pairs = ((b[i], b[j]) for i in range(len(b)) for j in range(i + 1, len(b)))
        for u, v in pairs:
            if coprime_only and gcd(u, v) != 1:
                continue
            k, rem = divmod(v - u, N)
            assert rem == 0 and k > 0
            rel.append((u, v, k))
    return rel


def f2_dependencies(rows: Sequence[int], ncols: int) -> Tuple[int, List[int]]:
    """Gaussian elimination over F2 on bit-rows.  Returns (rank, list of dependency masks), where a
    mask is an int whose set bits index the rows that XOR to zero."""
    pivots: Dict[int, Tuple[int, int]] = {}  # pivot column -> (row bits, combination mask)
    deps: List[int] = []
    for i, r in enumerate(rows):
        comb = 1 << i
        while r:
            col = r.bit_length() - 1
            if col in pivots:
                pr, pc = pivots[col]
                r ^= pr
                comb ^= pc
            else:
                pivots[col] = (r, comb)
                break
        if r == 0:
            deps.append(comb)
    return len(pivots), deps


def relation_row(u: int, v: int, primes: Sequence[int]) -> Tuple[int, List[int]]:
    """Parity row of u*v (the relation u ≡ v gives u*v ≡ v^2), and the exponent vector of u*v."""
    eu = exponent_vector(u, primes)
    ev = exponent_vector(v, primes)
    assert eu is not None and ev is not None
    e = [a + b for a, b in zip(eu, ev)]
    bits = 0
    for j, x in enumerate(e):
        if x & 1:
            bits |= 1 << j
    return bits, e


def square_from_dependency(mask: int, rels: Sequence[Tuple[int, int, int]], exps: Sequence[List[int]],
                           primes: Sequence[int], N: int) -> int:
    """For a dependency (set of relations with prod(u_i v_i) = W^2), return gcd(W - prod v_i, N).

    prod u_i ≡ prod v_i (mod N) gives prod(u_i v_i) ≡ (prod v_i)^2, i.e. W^2 ≡ V^2 (mod N)."""
    tot = [0] * len(primes)
    V = 1
    i = 0
    m = mask
    while m:
        if m & 1:
            for j, x in enumerate(exps[i]):
                tot[j] += x
            V = V * rels[i][1] % N
        m >>= 1
        i += 1
    W = 1
    for p, t in zip(primes, tot):
        assert t % 2 == 0
        W = W * pow(p, t // 2, N) % N
    return gcd((W - V) % N, N)


# ----------------------------------------------------------------------------------------------
# N-free smooth abc triples (operational; no N involved)
# ----------------------------------------------------------------------------------------------

def primitive_smooth_triples(vals: np.ndarray, H: int) -> List[Tuple[int, int, int]]:
    """All (a, b, c) with a <= b, a + b = c <= H, gcd(a, b) = 1, a, b, c in `vals` (sorted)."""
    vals = vals[vals <= H]
    out = []
    for c in vals[1:]:
        c = int(c)
        a = vals[vals <= c // 2]
        b = c - a
        idx = np.searchsorted(vals, b)
        idx[idx >= vals.size] = vals.size - 1
        hit = vals[idx] == b
        for aa in a[hit]:
            aa = int(aa)
            if gcd(aa, c) == 1:
                out.append((aa, c - aa, c))
    return out


# ----------------------------------------------------------------------------------------------
# radicals and abc quality
# ----------------------------------------------------------------------------------------------

def rad_smooth(n: int, primes: Sequence[int]) -> int:
    e = exponent_vector(n, primes)
    assert e is not None
    r = 1
    for p, x in zip(primes, e):
        if x:
            r *= p
    return r


def rad_diag(n: int) -> int:
    """Radical of an arbitrary positive integer (DIAGNOSTIC: uses sympy factorint)."""
    from sympy import factorint
    r = 1
    for p in factorint(n):
        r *= p
    return r


def powerful_part_ratio(k: int, radk: int) -> int:
    """k / rad(k) (exact)."""
    assert k % radk == 0
    return k // radk


# ----------------------------------------------------------------------------------------------
# subgroup order |H| = |<S>| in (Z/N)^*  (DIAGNOSTIC: uses p, q)
# ----------------------------------------------------------------------------------------------

def _primitive_root(p: int) -> int:
    from sympy import primitive_root
    return int(primitive_root(p))


def subgroup_order_diag(primes: Sequence[int], p: int, q: int) -> int:
    """|<primes>| inside (Z/pq)^* ≅ Z/(p-1) x Z/(q-1).  DIAGNOSTIC (uses the factors).

    H ≅ L / M with M = (p-1)Z x (q-1)Z and L = M + span{(dlog_g l, dlog_h l)}; |H| = det M / det L,
    det L = gcd of all 2x2 minors of the generating set."""
    from sympy.ntheory import discrete_log
    g, h = _primitive_root(p), _primitive_root(q)
    vecs = [(p - 1, 0), (0, q - 1)]
    for l in primes:
        if l % p == 0 or l % q == 0:
            raise ValueError("factor-base prime divides N")
        vecs.append((int(discrete_log(p, l % p, g)), int(discrete_log(q, l % q, h))))
    d = 0
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            d = gcd(d, abs(vecs[i][0] * vecs[j][1] - vecs[i][1] * vecs[j][0]))
    M = (p - 1) * (q - 1)
    assert M % d == 0
    return M // d


def subgroup_order_bruteforce_diag(primes: Sequence[int], N: int) -> int:
    """|<primes>| in (Z/N)^* by closure (DIAGNOSTIC; tiny N only, for tests)."""
    seen = {1}
    frontier = [1]
    while frontier:
        nxt = []
        for x in frontier:
            for l in primes:
                z = x * l % N
                if z not in seen:
                    seen.add(z)
                    nxt.append(z)
        frontier = nxt
    return len(seen)
