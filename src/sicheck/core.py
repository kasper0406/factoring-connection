"""sicheck core: Prop SI (smooth numbers in single intervals [kN, kN + h]) at toy scale.

Labelling (AGENT_RULES 2):
* OPERATIONAL routines see only N, the public parameters (y, h, constructor name, seed) and kN. They never use p, q,
  phi(N) or any factorisation routine: primes_upto, smooth_in_interval, exponent_vector, the constructors,
  collect_relations, Dependencies, square_root_of_unity, factor_si, beta_blum (uses only N via Jacobi symbols).
* DIAGNOSTIC routines (suffix ``_diag``) may use p, q; they live in src/sicheck/diag.py and are never imported here.

Exact integer arithmetic throughout. The interval sieve accumulates rounded integer logarithms (int32) only to select
candidates; the rounding slack is chosen so that no y-smooth integer can be missed, and every candidate is verified by
exact trial division. numpy is used only as an exact integer container.
"""
from __future__ import annotations

import math
import random
from math import gcd, isqrt
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

SCALE = 16  # integer log units per bit


# ----------------------------------------------------------------------------------------------------------------------
# primes, smoothness (operational)
# ----------------------------------------------------------------------------------------------------------------------

def primes_upto(y: int) -> List[int]:
    if y < 2:
        return []
    s = bytearray([1]) * (y + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(y) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
    return [i for i in range(y + 1) if s[i]]


def exponent_vector(n: int, primes: Sequence[int]) -> Optional[List[int]]:
    """Exponents of n over `primes` if n is smooth over them, else None (exact trial division). n >= 1."""
    if n <= 0:
        return None
    e = [0] * len(primes)
    for j, p in enumerate(primes):
        if n == 1:
            break
        while n % p == 0:
            n //= p
            e[j] += 1
    return e if n == 1 else None


def is_smooth(n: int, primes: Sequence[int]) -> bool:
    return exponent_vector(n, primes) is not None


def jacobi(a: int, n: int) -> int:
    """Jacobi symbol (a/n), n odd positive (exact)."""
    assert n > 0 and n & 1
    a %= n
    t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


def smooth_in_interval(A: int, H: int, primes: Sequence[int]) -> List[int]:
    """All integers n in [A, A + H) (A >= 2) whose prime factors lie in `primes`, sorted. Exact.

    Sieve: acc[i] = sum over prime powers p^e | A + i of w_p, w_p = round(SCALE log2 p). For a smooth n,
    acc >= SCALE log2 n - Omega(n)/2 >= SCALE log2 A - log2(A + H)/2, so the threshold below cannot miss one.
    Candidates are verified exactly.
    """
    assert A >= 2 and H >= 1
    top = A + H - 1
    acc = np.zeros(H, dtype=np.int32)
    for p in primes:
        w = int(round(SCALE * math.log2(p)))
        pe = p
        while pe <= top:
            start = (-A) % pe
            if start < H:
                acc[start::pe] += w
            pe *= p
    slack = math.ceil(0.5 * math.log2(top)) + 2
    thr = math.floor(SCALE * math.log2(A)) - slack
    idx = np.nonzero(acc >= thr)[0]
    out = []
    for i in idx.tolist():
        n = A + i
        if exponent_vector(n, primes) is not None:
            out.append(n)
    return out


# ----------------------------------------------------------------------------------------------------------------------
# constructors: functions of (X = kN, h, primes, rng). They see kN only (not N, not k separately), except the two
# explicitly labelled degenerate ones, which also see N (still operational: no p, q).
# ----------------------------------------------------------------------------------------------------------------------

def divisors_from_exps(primes: Sequence[int], e: Sequence[int], cap: int) -> List[int]:
    """All divisors <= cap of prod p^e (exact)."""
    divs = [1]
    for p, a in zip(primes, e):
        if a == 0:
            continue
        new = []
        for d in divs:
            x = d
            for _ in range(a):
                x *= p
                if x > cap:
                    break
                new.append(x)
        divs.extend(new)
    return divs


def choose_next(X: int, h: int, smooth: Sequence[int], primes, rng) -> Optional[int]:
    """Smallest y-smooth u in [X, X + h]."""
    return smooth[0] if smooth else None


def choose_last(X: int, h: int, smooth: Sequence[int], primes, rng) -> Optional[int]:
    """Largest y-smooth u in [X, X + h]."""
    return smooth[-1] if smooth else None


def choose_rand(X: int, h: int, smooth: Sequence[int], primes, rng) -> Optional[int]:
    """Uniformly random y-smooth u in [X, X + h] (rng seeded by the caller from (seed, X))."""
    return rng.choice(list(smooth)) if smooth else None


def cofactor_choice(X: int, h: int, smooth: Sequence[int], primes) -> Optional[Tuple[int, int]]:
    """(u0, u) for the greedy cofactor construction (see choose_cofactor), or None."""
    best = None
    for u in smooth:
        v = u - X
        e = exponent_vector(u, primes)
        cands = [d for d in divisors_from_exps(primes, e, h) if d > v]
        if cands:
            d = min(cands)
            if best is None or (d, u) < best:
                best = (d, u)
    return best


def choose_cofactor(X: int, h: int, smooth: Sequence[int], primes, rng) -> Optional[int]:
    """Greedy cofactor construction: the least smooth multiplier u0 <= h such that u0 * ceil(X/u0) is smooth.

    u0 * ceil(X/u0) lies in [X, X + u0); it is smooth iff it is one of `smooth` and u0 | u with u0 > u - X. So this is
    computed exactly from the smooth list: minimise u0 over pairs (u, u0), u0 | u, u - X < u0 <= h (ties: smaller u).
    Its residue is v = (-X) mod u0.
    """
    best = cofactor_choice(X, h, smooth, primes)
    return None if best is None else best[1]


CONSTRUCTORS: Dict[str, Callable] = {
    "next": choose_next,
    "last": choose_last,
    "rand": choose_rand,
    "cofactor": choose_cofactor,
}


def beta_blum(u: int, v: int, primes: Sequence[int], N: int) -> int:
    """Computable fruitfulness bit for Blum N (Prop B in the notes): ((Omega_T(u) - Omega_T(v)) / 2) mod 2,
    T = {l : J_N(l) = -1}. Uses only N. Requires J_N(u) = J_N(v), which holds for u = v (mod N)."""
    eu = exponent_vector(u, primes)
    ev = exponent_vector(v, primes)
    T = [jacobi(l, N) == -1 for l in primes]
    a = sum(x for x, t in zip(eu, T) if t)
    b = sum(x for x, t in zip(ev, T) if t)
    assert (a - b) % 2 == 0
    return ((a - b) // 2) % 2


def qs_factor_base(primes: Sequence[int], N: int) -> List[int]:
    """Degenerate factor base: primes l <= y with J_N(l) = +1 (the shape of a QS factor base for Blum N)."""
    return [l for l in primes if jacobi(l, N) == 1]


# ----------------------------------------------------------------------------------------------------------------------
# relations and linear algebra (operational)
# ----------------------------------------------------------------------------------------------------------------------

class Relation:
    __slots__ = ("k", "u", "v", "eu", "ev", "g", "key")

    def __init__(self, k: int, u: int, v: int, eu: List[int], ev: List[int]):
        self.k, self.u, self.v, self.eu, self.ev = k, u, v, eu, ev
        self.g = gcd(u, v)
        self.key = (u // self.g, v // self.g)

    def parity(self) -> int:
        bits = 0
        for j, (a, b) in enumerate(zip(self.eu, self.ev)):
            if (a + b) & 1:
                bits |= 1 << j
        return bits


def square_root_of_unity(rels: Sequence[Relation], mask: int, primes: Sequence[int], N: int) -> Tuple[int, int]:
    """For a dependency (bitmask over rels with sum of (eu + ev) even) return (Z, g): Z = W / V mod N with
    W = prod l^{E_l / 2}, V = prod v_i, so Z^2 = 1 (mod N); g = gcd(V, N) (if g > 1, Z = 0 and g is a factor)."""
    tot = [0] * len(primes)
    V = 1
    i = 0
    m = mask
    while m:
        if m & 1:
            r = rels[i]
            for j in range(len(primes)):
                tot[j] += r.eu[j] + r.ev[j]
            V = V * r.v % N
        m >>= 1
        i += 1
    W = 1
    for p, t in zip(primes, tot):
        assert t % 2 == 0
        if t:
            W = W * pow(p, t // 2, N) % N
    g = gcd(V, N)
    if g != 1:
        return 0, g
    Z = W * pow(V, -1, N) % N
    assert Z * Z % N == 1
    return Z, 1


def classify_root(Z: int, N: int) -> Tuple[bool, int]:
    """(fruitful, factor) for a square root of unity Z mod N. Fruitful iff Z not in {1, -1}."""
    if Z == 1 or Z == N - 1:
        return False, 1
    f = gcd(Z - 1, N)
    assert 1 < f < N
    return True, f


class Dependencies:
    """Incremental F2 elimination on relation parity rows (bit ints); collects a basis of the dependency space."""

    def __init__(self):
        self.pivots: Dict[int, Tuple[int, int]] = {}
        self.basis: List[int] = []
        self.n = 0

    def add(self, row: int) -> Optional[int]:
        comb = 1 << self.n
        self.n += 1
        while row:
            col = row.bit_length() - 1
            if col in self.pivots:
                pr, pc = self.pivots[col]
                row ^= pr
                comb ^= pc
            else:
                self.pivots[col] = (row, comb)
                return None
        self.basis.append(comb)
        return comb

    @property
    def rank(self) -> int:
        return len(self.pivots)


# ----------------------------------------------------------------------------------------------------------------------
# end-to-end algorithm (operational)
# ----------------------------------------------------------------------------------------------------------------------

def interval_rng(seed: int, X: int) -> random.Random:
    return random.Random(f"sicheck|{seed}|{X}")


def factor_si(N: int, y: int, h: int, constructor: str = "rand", seed: int = 0, excess: int = 20,
              k_max: int = 10 ** 6, n_random_deps: int = 64, factor_base: Optional[List[int]] = None,
              chooser: Optional[Callable] = None, trace: Optional[list] = None) -> Dict:
    """Prop SI with a brute-force (interval-sieving) constructor. Operational: sees N, y, h, seed only.

    For k = 1, 2, ...: the constructor receives X = kN and returns a smooth u in [X, X + h] (or None); v = u - X;
    keep (u, v) if v is smooth and its primitive pair (u, v)/gcd is new. Stop at #relations = #factor base + excess.
    Then F2 dependencies, square roots of unity, gcd. Returns the factor (or None) and counters.
    """
    primes = factor_base if factor_base is not None else primes_upto(y)
    g = gcd(N, math.prod(primes)) if primes else 1
    if g != 1:  # ACCIDENTAL: N has a factor <= y (AGENT_RULES 3)
        return {"factor": g, "kind": "ACCIDENTAL_small_factor"}
    choose = chooser if chooser is not None else CONSTRUCTORS[constructor]
    need = len(primes) + excess
    rels: List[Relation] = []
    seen = set()
    stats = {"intervals": 0, "empty": 0, "v_smooth": 0, "duplicates": 0, "nonprimitive": 0, "sieved": 0,
             "smooth_u_total": 0}
    deps = Dependencies()
    k = 0
    while len(rels) < need and k < k_max:
        k += 1
        X = k * N
        sm = smooth_in_interval(X, h + 1, primes)
        stats["intervals"] += 1
        stats["sieved"] += h + 1
        stats["smooth_u_total"] += len(sm)
        u = choose(X, h, sm, primes, interval_rng(seed, X))
        if u is None:
            stats["empty"] += 1
            continue
        v = u - X
        if trace is not None:
            trace.append((k, v))
        if v == 0:
            continue
        ev = exponent_vector(v, primes)
        if ev is None:
            continue
        stats["v_smooth"] += 1
        r = Relation(k, u, v, exponent_vector(u, primes), ev)
        if r.g > 1:
            stats["nonprimitive"] += 1
        if r.key in seen:
            stats["duplicates"] += 1
            continue
        seen.add(r.key)
        rels.append(r)
        deps.add(r.parity())
    stats.update({"k_last": k, "relations": len(rels), "rank": deps.rank, "dim_deps": len(deps.basis),
                  "fb_size": len(primes)})
    # basis roots, then random combinations
    rng = random.Random(f"sicheck-deps|{seed}|{N}")
    factor = None
    phi_basis = []
    for mask in deps.basis:
        Z, gg = square_root_of_unity(rels, mask, primes, N)
        if gg != 1:
            factor = factor or gg
            phi_basis.append(None)
            continue
        fr, f = classify_root(Z, N)
        phi_basis.append(int(fr))
        if fr and factor is None:
            factor = f
    fruitful = 0
    tried = 0
    for _ in range(n_random_deps if deps.basis else 0):
        mask = 0
        while mask == 0:
            for b in deps.basis:
                if rng.getrandbits(1):
                    mask ^= b
        Z, gg = square_root_of_unity(rels, mask, primes, N)
        tried += 1
        if gg == 1:
            fr, _ = classify_root(Z, N)
            fruitful += int(fr)
    stats.update({"phi_basis": phi_basis, "phi_nonzero": any(x == 1 for x in phi_basis),
                  "random_deps": tried, "random_fruitful": fruitful})
    return {"factor": factor, "kind": "STRUCTURAL" if factor else None, "stats": stats, "relations": rels,
            "deps": deps, "primes": primes}
