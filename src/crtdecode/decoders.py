"""Operational toy decoders for crtdecode E1/E2 (pre-registered in notes/crtdecode_prereg.md).

OPERATIONAL code: every decoder receives only (X, h, y), a public repertoire fixed independently of X, and a budget.
None of them sees the ground-truth sieve. Smoothness tests are exact trial division (`core.is_smooth`).
Costs: 1 per pointer (computing the multiples of one d in the window), 1 per test of one position.

E2's information decoder interacts with X ONLY through `JohnsonOracle.call` (idealised Theorem-J certificates);
the oracle is the harness.
"""
from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Sequence, Tuple

from .core import certified, is_smooth, prime_powers_upto, primes_upto


# ------------------------------------------------------------------------------------------------------ E1

class Result:
    __slots__ = ("cost", "found", "z", "tests", "pointers", "tested_log")

    def __init__(self):
        self.cost = 0
        self.found = False
        self.z = None
        self.tests = 0
        self.pointers = 0
        self.tested_log: List[Tuple[int, bool]] = []  # (hits at test time, smooth?)


def _multiples(X: int, h: int, d: int) -> List[int]:
    j0 = (-X) % d
    return list(range(j0, h + 1, d))


def divisor_trick(X: int, h: int, y: int, repertoire: Sequence[int], budget: int) -> Result:
    """DT / DT8: walk the repertoire; pay 1 per pointer, 1 per test of each untested multiple."""
    r = Result()
    tested = set()
    for d in repertoire:
        if r.cost + 1 > budget:
            return r
        r.cost += 1
        r.pointers += 1
        for j in _multiples(X, h, d):
            if j in tested:
                continue
            if r.cost + 1 > budget:
                return r
            r.cost += 1
            r.tests += 1
            tested.add(j)
            s = is_smooth(X + j, y)
            r.tested_log.append((1, s))
            if s:
                r.found, r.z = True, X + j
                return r
    return r


def heavy_hitter(X: int, h: int, y: int, repertoire: Sequence[int], budget: int, theta: float) -> Result:
    """PH-theta: adaptive pointer aggregation. A hit position j is tested once ln L[j] >= theta * ln h,
    L[j] = lcm of the repertoire elements hitting j. After the repertoire is exhausted, remaining hit positions
    are tested in decreasing order of L[j]."""
    r = Result()
    L: Dict[int, int] = {}
    hits: Dict[int, int] = {}
    tested = set()
    # threshold on L as an exact integer comparison: L >= ceil(h^theta) (theta rational with small denominator)
    thr = _int_pow_ceil(h, theta)

    def test(j: int) -> bool:
        r.cost += 1
        r.tests += 1
        tested.add(j)
        s = is_smooth(X + j, y)
        r.tested_log.append((hits[j], s))
        if s:
            r.found, r.z = True, X + j
        return s

    for d in repertoire:
        if r.cost + 1 > budget:
            return r
        r.cost += 1
        r.pointers += 1
        for j in _multiples(X, h, d):
            Lj = L.get(j, 1)
            L[j] = Lj * d // math.gcd(Lj, d)
            hits[j] = hits.get(j, 0) + 1
            if j not in tested and L[j] >= thr:
                if r.cost + 1 > budget:
                    return r
                if test(j):
                    return r
    for j in sorted((j for j in L if j not in tested), key=lambda j: (-L[j], j)):
        if r.cost + 1 > budget:
            return r
        if test(j):
            return r
    return r


def _int_pow_ceil(h: int, theta: float) -> int:
    """ceil(h^theta) exactly for theta in {1, 1.25, 1.5, ...} with h a power of two (as in the configs)."""
    k = h.bit_length() - 1
    if 1 << k != h:
        raise ValueError("h must be a power of two")
    num = theta * 4
    if abs(num - round(num)) > 1e-12:
        raise ValueError("theta must be a multiple of 1/4")
    num = int(round(num))
    # h^theta = 2^(k*num/4); ceil of 2^(a/4) for integer a
    a = k * num
    if a % 4 == 0:
        return 1 << (a // 4)
    # ceil(2^(a/4)) = least integer m with m^4 >= 2^a
    lo, hi = 1 << (a // 4), 1 << (a // 4 + 1)
    while lo < hi:
        mid = (lo + hi) // 2
        if mid ** 4 >= (1 << a):
            hi = mid
        else:
            lo = mid + 1
    return lo


def random_tests(X: int, h: int, y: int, order: Sequence[int], budget: int) -> Result:
    """RT reference: test positions in a fixed random order."""
    r = Result()
    for j in order:
        if r.cost + 1 > budget:
            return r
        r.cost += 1
        r.tests += 1
        if is_smooth(X + j, y):
            r.found, r.z = True, X + j
            return r
    return r


# ------------------------------------------------------------------------------------------------------ E2

class JohnsonOracle:
    """HARNESS. Holds X; answers idealised unweighted Johnson calls (Theorem J's class).

    call(M, lo, length): region [X + lo, X + lo + length - 1]; returns [(offset, gcd(X + offset, M))] for every
    certified point. Counts calls.
    """

    def __init__(self, X: int):
        self._X = X
        self.calls = 0

    def call(self, M: int, lo: int, length: int) -> List[Tuple[int, int]]:
        self.calls += 1
        X = self._X
        out = []
        # Only points with gcd >= 2 can be certified; enumerate candidates through divisors of M is expensive,
        # so scan the multiples of each prime-power part when M is a prime power, else of M itself.
        if length == M:
            j = (-(X + lo)) % M
            out.append((lo + j, M))  # the unique multiple of M in a region of length M; gcd = M, certified
            return out
        q = _pb(M)
        if q is not None:
            m = q
            seen = {}
            while M % m == 0:
                j0 = (-(X + lo)) % m
                for j in range(j0, length, m):
                    seen[lo + j] = m
                m *= q
            for off, g in sorted(seen.items()):
                if certified(g, M, length):
                    out.append((off, g))
            return out
        # general M: only the full-gcd points are needed by the tree (final certification call)
        j0 = (-(X + lo)) % M
        for j in range(j0, length, M):
            if certified(M, M, length):
                out.append((lo + j, M))
        return out


def _pb(m: int) -> Optional[int]:
    for p in primes_upto(1 << 12):
        if m % p == 0:
            while m % p == 0:
                m //= p
            return p if m == 1 else None
    return None


def information_decoder(oracle: JohnsonOracle, X_public_size_bound: int, h: int, y: int) -> Dict[str, object]:
    """Decision tree of idealised Johnson calls (Prop DT). Uses X only through oracle answers.

    X_public_size_bound: an upper bound on X + h (public: the bit size), used only to cap exponents e_q.
    Returns depth, the reconstructed smooth offsets, the leaf cost (incidences) and the output.
    """
    V = h + 1
    primes = primes_upto(y)
    # Step 1: residues mod q^{a_q} <= h + 1, grouped into products M <= h + 1 (greedy, descending).
    qa = []
    for q in primes:
        a, m = 0, 1
        while m * q <= V:
            m *= q
            a += 1
        qa.append((q, a, m))
    qa.sort(key=lambda t: -t[2])
    groups: List[List[Tuple[int, int, int]]] = []
    cur, prod = [], 1
    for t in qa:
        if prod * t[2] <= V:
            cur.append(t)
            prod *= t[2]
        else:
            groups.append(cur)
            cur, prod = [t], t[2]
    if cur:
        groups.append(cur)
    res: Dict[int, int] = {}  # q -> (-X) mod q^{a_q}, as offset residue r with X + r == 0 mod q^{a_q}
    for g in groups:
        M = 1
        for (_q, _a, m) in g:
            M *= m
        ans = oracle.call(M, 0, M)
        (off, gg), = ans
        assert gg == M
        for (q, _a, m) in g:
            res[q] = off % m
    # Step 2: higher powers on the full window.
    high: Dict[int, Dict[int, int]] = {q: {} for q in primes}  # q -> {offset: v_q} for v_q > a_q
    for (q, a, m) in qa:
        e = a + 1
        mq = m * q
        while mq <= X_public_size_bound:
            ans = oracle.call(mq, 0, V)
            got = [off for (off, g) in ans if g == mq]
            if not got:
                break
            for off in got:
                high[q][off] = e
            e += 1
            mq *= q
    # Step 3 (leaf computation): rebuild ln of the smooth part at every offset from residues; count incidences.
    logpart = [0] * V  # exact: store the smooth part as an integer product
    part = [1] * V
    inc = 0
    for (q, a, m) in qa:
        r = res[q]
        # offsets with v_q >= b for b <= a: j == r mod q^b
        for b in range(1, a + 1):
            mb = q ** b
            for j in range(r % mb, V, mb):
                part[j] *= q
                inc += 1
        for off, v in high[q].items():
            for _ in range(v - a):
                part[off] *= q
                inc += 1
    return {"part": part, "incidences": inc, "groups": len(groups), "calls_before_output": oracle.calls}


def information_decode_and_certify(X: int, h: int, y: int, bits: int) -> Dict[str, object]:
    """Run the tree on X in [2^(bits-1), 2^bits); output the least smooth offset; certify it with one more call.

    Exact leaf rule: part[j] is the full y-smooth part of X + j (all valuations are known), so the cofactor
    c = (X + j)/part[j] is 1 or has all prime factors > y >= 2, i.e. c >= 3. Since 2^(bits-1) <= X + j <= 2^bits + h
    < 3 * 2^(bits-1), X + j is y-smooth iff part[j] >= 2^(bits-1). The final call (modulus part[j], region = the
    window) certifies it: gcd = modulus >= h + 1.
    """
    assert (1 << (bits - 1)) <= X < (1 << bits) and h < (1 << (bits - 1))
    orc = JohnsonOracle(X)
    out = information_decoder(orc, (1 << bits) + h, h, y)
    part = out["part"]
    lo = 1 << (bits - 1)
    smooth_offsets = [j for j, pj in enumerate(part) if pj >= lo]
    output = None
    if smooth_offsets:
        j = smooth_offsets[0]
        ans = orc.call(part[j], 0, h + 1)
        if any(off == j for (off, _g) in ans):
            output = j
    return {"depth": orc.calls, "smooth_offsets": smooth_offsets, "output": output,
            "incidences": out["incidences"], "groups": out["groups"]}
