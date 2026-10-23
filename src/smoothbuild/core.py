"""smoothbuild core helpers (agent `smoothbuild`, phase 24).

Sections
--------
* OPERATIONAL (sees only N, public config and a seed): ``primes_upto``, ``trial_factor``, ``gauss_split``,
  ``hidden_root_class_factor`` (the Dixon-style reduction of Prop EQ1: an average-case smooth-class constructor
  gives a congruence of squares; instantiated with the trivial Gauss-split constructor).
* THEORY / EXACT CHECKS (no N): ``burnside_rank`` (Galois rank r_G), ``root_count`` (rho_F(p)),
  ``atom_product_count`` and ``psi_count`` (Lemma AT), ``left_kernel_rank`` (Lemma KI), ``inv_square_sum_fixed``
  (exact fixed-point sum of 1/p^2).

Nothing here calls a factoring routine on N. ``trial_factor`` only divides small pieces by factor-base primes.
"""
from __future__ import annotations

import math
import random
from itertools import product as iproduct
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------------------------------------------
# primes and smoothness (operational)


def primes_upto(n: int) -> List[int]:
    """All primes <= n (sieve of Eratosthenes)."""
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]


def trial_factor(z: int, primes: Sequence[int]) -> Optional[Dict[int, int]]:
    """Factor |z| over ``primes``; return {p: e} if |z| is smooth over them (z != 0), else None."""
    z = abs(z)
    if z == 0:
        return None
    out: Dict[int, int] = {}
    for p in primes:
        if p * p > z:
            break
        if z % p == 0:
            e = 0
            while z % p == 0:
                z //= p
                e += 1
            out[p] = e
    if z > 1:
        if z <= primes[-1]:
            out[z] = out.get(z, 0) + 1
        else:
            return None
    return out


def gauss_split(t: int, N: int) -> Tuple[int, int]:
    """Rational reconstruction: (u, v) with u = t*v (mod N), |u| < sqrt(N), |v| <= sqrt(N) (approximately).

    Extended Euclid on (N, t): r_i = s_i * t (mod N); stop at the first r_i with r_i^2 < N.
    """
    t %= N
    r0, r1 = N, t
    s0, s1 = 0, 1
    while r1 * r1 >= N:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    return r1, s1


# ---------------------------------------------------------------------------------------------------------------
# Prop EQ1 (operational): hidden-root reduction from smooth-class construction to a congruence of squares


def _gf2_dependency(rows: List[int]) -> Optional[int]:
    """Return a nonempty subset (bitmask over row indices) whose rows XOR to zero, or None."""
    pivots: Dict[int, Tuple[int, int]] = {}
    for i, r in enumerate(rows):
        comb = 1 << i
        while r:
            lead = r.bit_length() - 1
            if lead in pivots:
                pr, pc = pivots[lead]
                r ^= pr
                comb ^= pc
            else:
                pivots[lead] = (r, comb)
                break
        if r == 0:
            return comb
    return None


def hidden_root_class_factor(N: int, y: int, seed: int,
                             constructor: Optional[Callable[[int, int], Tuple[int, int]]] = None,
                             max_targets: int = 10 ** 6) -> Dict[str, object]:
    """OPERATIONAL. Factor N from a smooth-class constructor applied to random squares (Prop EQ1).

    The constructor sees only (N, t) with t = r^2 mod N and returns (u, v) with u = t*v (mod N). A relation is kept
    when u and v are y-smooth. Since the constructor never sees r, each r is a uniformly random square root of t
    given everything the constructor saw; a dependency gives a nontrivial square root with probability >= 1/2 for
    N with >= 2 distinct odd prime factors. Every gcd is checked; any nonunit met is returned as a factor.
    Default constructor: ``gauss_split`` (so this is Dixon's method with rational reconstruction).
    """
    if constructor is None:
        constructor = lambda NN, tt: gauss_split(tt, NN)  # noqa: E731
    rng = random.Random(f"smoothbuild-eq1|{N}|{y}|{seed}")
    P = primes_upto(y)
    idx = {p: i + 1 for i, p in enumerate(P)}  # column 0 = sign
    K = len(P) + 1
    rels: List[Tuple[int, Dict[int, int], int]] = []  # (r, exponent dict incl. -1 key, row parity)
    counters = {"targets": 0, "splits": 0, "tests": 0, "gcd": 0, "relations": 0, "dependencies": 0}
    for _ in range(max_targets):
        counters["targets"] += 1
        r = rng.randrange(1, N)  # full unit range, so r is uniform among the roots of t (codex adversarial #6)
        g = math.gcd(r, N)
        counters["gcd"] += 1
        if g != 1:
            return {"factor": g, "kind": "ACCIDENTAL", **counters}
        t = r * r % N
        u, v = constructor(N, t)
        counters["splits"] += 1
        if u == 0 or v == 0:
            continue
        counters["tests"] += 1
        fu = trial_factor(u, P)
        if fu is None:
            continue
        counters["tests"] += 1
        fv = trial_factor(v, P)
        if fv is None:
            continue
        # relation: r^2 * v = u (mod N)  ->  r^2 = u * v^{-1}; exponent vector e = e(u) - e(v)
        e: Dict[int, int] = dict(fu)
        for p, k in fv.items():
            e[p] = e.get(p, 0) - k
        sgn = (1 if u < 0 else 0) ^ (1 if v < 0 else 0)
        row = sgn
        for p, k in e.items():
            if k % 2:
                row |= 1 << idx[p]
        rels.append((r, e, sgn))
        counters["relations"] += 1
        if len(rels) < K + 1 and len(rels) % 8:
            continue
        rows = [0] * len(rels)
        for i, (_, ee, sg) in enumerate(rels):
            rr = sg
            for p, k in ee.items():
                if k % 2:
                    rr |= 1 << idx[p]
            rows[i] = rr
        # try dependencies, removing one member each time to get fresh ones
        work = list(range(len(rels)))
        for _attempt in range(8):
            dep = _gf2_dependency([rows[i] for i in work])
            if dep is None:
                break
            counters["dependencies"] += 1
            S = [work[j] for j in range(len(work)) if dep >> j & 1]
            X = 1
            tot: Dict[int, int] = {}
            neg = 0
            for i in S:
                ri, ee, sg = rels[i]
                X = X * ri % N
                neg ^= sg
                for p, k in ee.items():
                    tot[p] = tot.get(p, 0) + k
            assert neg == 0 and all(k % 2 == 0 for k in tot.values())
            Y = 1
            for p, k in tot.items():
                if k >= 0:
                    Y = Y * pow(p, k // 2, N) % N
                else:
                    gp = math.gcd(p, N)
                    counters["gcd"] += 1
                    if gp != 1:
                        return {"factor": gp, "kind": "ACCIDENTAL", **counters}
                    Y = Y * pow(pow(p, -1, N), (-k) // 2, N) % N
            g = math.gcd((X - Y) % N, N)
            counters["gcd"] += 1
            if 1 < g < N:
                return {"factor": g, "kind": "STRUCTURAL", **counters}
            work.remove(S[0])
    return {"factor": None, "kind": None, **counters}


# ---------------------------------------------------------------------------------------------------------------
# theory helpers (no N)


def burnside_rank(gens: Sequence[Sequence[int]], d: int) -> Tuple[int, int, int]:
    """For the permutation group G on {0..d-1} generated by ``gens`` (each a list image), return
    (|G|, #orbits on points, rank = #orbits on ordered pairs), via Burnside: rank = (1/|G|) sum fix(g)^2."""
    ident = tuple(range(d))
    G = {ident}
    frontier = [ident]
    gens_t = [tuple(g) for g in gens]
    while frontier:
        new = []
        for h in frontier:
            for g in gens_t:
                c = tuple(g[h[i]] for i in range(d))
                if c not in G:
                    G.add(c)
                    new.append(c)
        frontier = new
    fix = [sum(1 for i in range(d) if g[i] == i) for g in G]
    n = len(G)
    s1, s2 = sum(fix), sum(f * f for f in fix)
    assert s1 % n == 0 and s2 % n == 0
    return n, s1 // n, s2 // n


# explicit transitive groups used in E1 (generators as image lists)
GROUPS = {
    "S4": [[1, 0, 2, 3], [1, 2, 3, 0]],
    "A4": [[1, 2, 0, 3], [0, 2, 3, 1]],
    "D4": [[1, 2, 3, 0], [2, 1, 0, 3]],
    "C4": [[1, 2, 3, 0]],
    "V4": [[1, 0, 3, 2], [2, 3, 0, 1]],
    "S6": [[1, 0, 2, 3, 4, 5], [1, 2, 3, 4, 5, 0]],
    "C6": [[1, 2, 3, 4, 5, 0]],
}

# E1 families: name -> (coefficients low-to-high, Galois group name)
FAMILIES = {
    "S4:x4+x+1": ([1, 1, 0, 0, 1], "S4"),
    "A4:x4+8x+12": ([12, 8, 0, 0, 1], "A4"),
    "D4:x4-2": ([-2, 0, 0, 0, 1], "D4"),
    "C4:Phi5": ([1, 1, 1, 1, 1], "C4"),
    "V4:Phi8": ([1, 0, 0, 0, 1], "V4"),
    "S6:x6+x+1": ([1, 1, 0, 0, 0, 0, 1], "S6"),
    "C6:Phi7": ([1, 1, 1, 1, 1, 1, 1], "C6"),
    "C6:Phi9": ([1, 0, 0, 1, 0, 0, 1], "C6"),
}


def poly_eval(coeffs: Sequence[int], a: int) -> int:
    v = 0
    for c in reversed(coeffs):
        v = v * a + c
    return v


def root_list(coeffs: Sequence[int], p: int) -> List[int]:
    """Distinct roots of F mod p (python-flint nmod_poly)."""
    import flint
    f = flint.nmod_poly([c % p for c in coeffs], p)
    if f.degree() < 1:
        return []
    return sorted(int(r) for r, _m in f.roots())


def root_count(coeffs: Sequence[int], p: int) -> int:
    return len(root_list(coeffs, p))


def inv_square_sum_fixed(primes: Iterable[int], scale_bits: int = 160) -> int:
    """Exact fixed-point floor(2^scale / p^2) summed: sum_p 1/p^2 to within #primes * 2^-scale."""
    S = 1 << scale_bits
    return sum(S // (p * p) for p in primes)


# ---------------------------------------------------------------------------------------------------------------
# Lemma AT and Lemma KI (exact checks)


def atom_product_count(atoms: Sequence[int], X: int) -> int:
    """Number of distinct integers <= X that are products of the given atoms (with repetition)."""
    seen = {1}
    stack = [1]
    while stack:
        m = stack.pop()
        for a in atoms:
            z = m * a
            if z <= X and z not in seen:
                seen.add(z)
                stack.append(z)
    return len(seen)


def psi_count(X: int, y: int) -> int:
    """Psi(X, y): number of y-smooth integers in [1, X] (exact, by enumeration)."""
    return atom_product_count(primes_upto(y), X)


def left_kernel_rank(rows: Sequence[int]) -> int:
    """dim of the GF(2) left kernel of the row list (rows as int bitmasks) = #rows - rank."""
    pivots: Dict[int, int] = {}
    rank = 0
    for r in rows:
        while r:
            lead = r.bit_length() - 1
            if lead in pivots:
                r ^= pivots[lead]
            else:
                pivots[lead] = r
                rank += 1
                break
    return len(rows) - rank
