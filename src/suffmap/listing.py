"""Search-to-listing reduction and toy exact checks (agent suffmap; notes/suffmap_results.md, Lemma SL).

Lemma SL [PROVED]. Let S be a set of integers and O a search oracle that, given an interval [a, b], returns some
element of S in [a, b] or reports none. Then all K elements of S in [a, b] can be listed with at most 2K + 1 calls
(and with at most 2K + J calls for a union of J disjoint intervals).

Proof: recurse on [a, z-1] and [z+1, b] after each hit z, never calling on empty intervals. Each element is
returned exactly once (the sub-intervals are disjoint and exclude z), so the recursion tree has K internal (hit)
nodes, each with at most two children, and the leaves are misses; a tree with K internal nodes of out-degree <= 2
has <= K + 1 leaves.

Everything here is HARNESS / DIAGNOSTIC: the toy oracle decides smoothness by exact trial division and is not an
operational factoring routine. It never sees any factor of N.
"""
from __future__ import annotations

from typing import Callable, List, Optional, Tuple

Oracle = Callable[[int, int], Optional[int]]


def is_smooth(n: int, y: int) -> bool:
    """Exact y-smoothness of |n| (n != 0) by trial division (harness)."""
    n = abs(n)
    if n == 0:
        return False
    d = 2
    while d <= y and d * d <= n:
        while n % d == 0:
            n //= d
        d += 1 if d == 2 else 2
    return n <= y


def make_smooth_search_oracle(y: int, pick: str = "first") -> Tuple[Oracle, dict]:
    """A search oracle for y-smooth integers; counts calls. pick in {'first', 'last', 'middle'} varies which hit
    is returned, to show the call bound does not depend on the oracle's choice."""
    stats = {"calls": 0}

    def oracle(a: int, b: int) -> Optional[int]:
        stats["calls"] += 1
        hits = [n for n in range(a, b + 1) if n != 0 and is_smooth(n, y)]
        if not hits:
            return None
        if pick == "first":
            return hits[0]
        if pick == "last":
            return hits[-1]
        return hits[len(hits) // 2]

    return oracle, stats


def list_by_search(a: int, b: int, oracle: Oracle) -> List[int]:
    """Lemma SL: list every element of S in [a, b] using the search oracle (iterative, no empty calls)."""
    out: List[int] = []
    stack = [(a, b)]
    while stack:
        lo, hi = stack.pop()
        if lo > hi:
            continue
        z = oracle(lo, hi)
        if z is None:
            continue
        out.append(z)
        stack.append((lo, z - 1))
        stack.append((z + 1, hi))
    return sorted(out)


def check_lemma_sl(a: int, b: int, y: int, pick: str = "first") -> dict:
    oracle, stats = make_smooth_search_oracle(y, pick)
    got = list_by_search(a, b, oracle)
    truth = [n for n in range(a, b + 1) if n != 0 and is_smooth(n, y)]
    return {"K": len(truth), "calls": stats["calls"], "bound": 2 * len(truth) + 1, "complete": got == truth}


def rational_side_rows(m: int, E: int, y: int, bmax: int) -> dict:
    """Toy NFS rational side: list y-smooth values of a - b m, |a| <= E, 1 <= b <= bmax, row by row (Lemma SL with
    J = bmax intervals). Returns outputs, oracle calls and the box size; calls <= 2K + bmax (proved)."""
    oracle, stats = make_smooth_search_oracle(y)
    K = 0
    for b in range(1, bmax + 1):
        # values a - b m for |a| <= E form the interval [-E - b m, E - b m]
        K += len(list_by_search(-E - b * m, E - b * m, oracle))
    return {"K": K, "calls": stats["calls"], "bound": 2 * K + bmax, "box": (2 * E + 1) * bmax}


class AdversarialPromiseOracle:
    """Prop PC (a): a promise constructor valid only on intervals of length >= h0. Adversary: return an already
    revealed element if the query contains one; otherwise reveal the element of S nearest the midpoint of the
    maximal revealed-free gap containing the query. If every gap of S is <= h0/4, all revealed points are pairwise
    >= h0/4 apart, so at most 4H/h0 + 1 distinct outputs can ever be obtained from a range of length H."""

    def __init__(self, S: List[int], lo: int, hi: int, h0: int):
        self.S, self.lo, self.hi, self.h0 = sorted(S), lo, hi, h0
        self.revealed: List[int] = []
        self.calls = 0

    def __call__(self, a: int, b: int) -> int:
        if b - a + 1 < self.h0:
            raise ValueError("outside the promise")
        self.calls += 1
        for r in self.revealed:
            if a <= r <= b:
                return r
        left = max([r for r in self.revealed if r < a], default=self.lo - 1)
        right = min([r for r in self.revealed if r > b], default=self.hi + 1)
        mid2 = left + right  # twice the gap midpoint (exact)
        cands = [s for s in self.S if a <= s <= b]
        z = min(cands, key=lambda s: (abs(2 * s - mid2), s))
        self.revealed.append(z)
        return z


def tiling_outputs(lo: int, hi: int, h0: int, oracle: Callable[[int, int], int]) -> List[int]:
    """Prop PC (b): disjoint tiles of length h0 give >= floor(H/h0) distinct outputs from an honest constructor."""
    out = []
    a = lo
    while a + h0 - 1 <= hi:
        out.append(oracle(a, a + h0 - 1))
        a += h0
    return out


def congruence_max_lower_bound(N: int) -> bool:
    """Elementary: u != v, u = v (mod N) implies max(|u|, |v|) >= N/2 (|u| + |v| >= |u - v| >= N). Exhaustive
    check over |u|, |v| < N for small N."""
    for u in range(-N + 1, N):
        for v in range(-N + 1, N):
            if u != v and (u - v) % N == 0 and 2 * max(abs(u), abs(v)) < N:
                return False
    return True


if __name__ == "__main__":  # pragma: no cover
    for pick in ("first", "last", "middle"):
        print(pick, check_lemma_sl(10 ** 6, 10 ** 6 + 2000, 50, pick))
    print(rational_side_rows(m=10 ** 5 + 3, E=200, y=60, bmax=20))
