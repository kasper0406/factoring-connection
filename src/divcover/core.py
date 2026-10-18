"""divcover core: generalized arithmetic progressions (GAPs) and divisor covers.

Diagnostic / combinatorial code for the Umans-Wang divisor-cover question. Nothing
here factors an RSA-type N; the only factoring is of GAP elements (search objects).

A GAP is Q = {b + sum_k x_k c_k : 0 <= x_k < L_k}.  Its *box size* is prod L_k (the
cost that Umans-Wang Lemma 5.3/5.4 pay), its *height* is max |element|.
Covers:
  * full  : every d in [1, n] divides some NONZERO element (UW Def. 3.1 on S - T),
  * prime : every prime p <= n divides some nonzero element (what UW Thm 5.2 uses;
            see notes/divcover_results.md, Prop. R1).
"""
from __future__ import annotations

import itertools
import math
import os
import subprocess
import tempfile
from dataclasses import dataclass
from functools import lru_cache

import sympy


@dataclass(frozen=True)
class GAP:
    b: int
    c: tuple
    L: tuple

    @property
    def rank(self) -> int:
        return len(self.c)

    @property
    def box(self) -> int:
        return math.prod(self.L)

    def elements(self):
        """Iterate over the box (with multiplicity)."""
        for x in itertools.product(*[range(l) for l in self.L]):
            yield self.b + sum(xi * ci for xi, ci in zip(x, self.c))

    def height(self) -> int:
        return max(abs(e) for e in self.elements())

    def residues(self, d: int) -> set:
        """Exact set Q mod d, computed as an iterated sumset mod d (cost O(d * rank))."""
        cur = {self.b % d}
        for ck, lk in zip(self.c, self.L):
            step = ck % d
            nxt = set()
            for s in cur:
                v = s
                for _ in range(min(lk, d)):
                    nxt.add(v)
                    v = (v + step) % d
            cur = nxt
        return cur


def nonzero_multiple_exists(Q: GAP, d: int) -> bool:
    """Does d divide some NONZERO element of Q? Exact; zero elements are excluded."""
    for e in Q.elements():
        if e != 0 and e % d == 0:
            return True
    return False


def covers(Q: GAP, n: int, mode: str = "full") -> bool:
    """Exact cover test by enumerating elements (use for small boxes / tests)."""
    els = [e for e in set(Q.elements()) if e != 0]
    targets = sympy.primerange(2, n + 1) if mode == "prime" else range(2, n + 1)
    for d in targets:
        if not any(e % d == 0 for e in els):
            return False
    return True


def uncovered(Q: GAP, n: int, mode: str = "prime") -> list:
    """List of targets not covered (element enumeration; zero elements excluded)."""
    els = [e for e in set(Q.elements()) if e != 0]
    targets = sympy.primerange(2, n + 1) if mode == "prime" else range(2, n + 1)
    return [d for d in targets if not any(e % d == 0 for e in els)]


def brute_min_height(n: int, shape: tuple, Hmax: int, mode: str = "prime"):
    """Reference brute force over positive GAPs (b>=1, c_k>=1) of the given shape.
    Returns (H*, GAP) or (None, None). Only for tiny parameters (tests)."""
    r = len(shape)
    w = [l - 1 for l in shape]
    best, arg = None, None
    for c in itertools.product(range(1, Hmax + 1), repeat=r):
        base = sum(wi * ci for wi, ci in zip(w, c))
        if base + 1 > Hmax:
            continue
        for b in range(1, Hmax - base + 1):
            h = b + base
            if best is not None and h >= best:
                break
            Q = GAP(b, tuple(c), tuple(shape))
            if covers(Q, n, mode):
                best, arg = h, Q
                break
    return best, arg


# ---------------------------------------------------------------- C search wrapper
_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gapsearch.c")


@lru_cache(maxsize=1)
def _binary() -> str:
    out_dir = os.path.join(tempfile.gettempdir(), "divcover_build")
    os.makedirs(out_dir, exist_ok=True)
    exe = os.path.join(out_dir, "gapsearch")
    if (not os.path.exists(exe)) or os.path.getmtime(exe) < os.path.getmtime(_SRC):
        subprocess.run(["gcc", "-O2", "-o", exe, _SRC], check=True)
    return exe


def c_min_height(n: int, shape: tuple, Hmax: int, mode: str = "prime", timeout=None):
    """Exact minimal height via the C search. Returns dict(best, b, c, nparams)."""
    m = {"prime": 0, "full": 1}[mode]
    args = [_binary(), str(n), str(m), str(Hmax)] + [str(l) for l in shape]
    res = subprocess.run(args, capture_output=True, text=True, check=True, timeout=timeout)
    kv = dict(tok.split("=") for tok in res.stdout.split())
    best = int(kv["best"])
    cs = tuple(int(v) for v in kv["c"].split(","))[: len(shape)]
    return {"best": best or None, "b": int(kv["b"]), "c": cs, "nparams": int(kv["nparams"])}


# ---------------------------------------------------------------- random model
def residue_cover_prob(p: int, shape: tuple) -> float:
    """P_p = Pr_{b,c uniform mod p}[0 in Q mod p] for a GAP of the given shape (exact)."""
    r = len(shape)
    hit = 0
    tot = 0
    for c in itertools.product(range(p), repeat=r):
        Q0 = GAP(0, c, shape)
        res = Q0.residues(p)  # set of <x,c> mod p
        hit += len(res)  # number of b with -b in res  == |res|
        tot += p
    return hit / tot


def count_gaps(H: float, shape: tuple) -> float:
    """#{b>=1, c_k>=1 : b + sum (L_k-1) c_k <= H} ~ H^{r+1} / ((r+1)! prod (L_k-1))."""
    r = len(shape)
    return H ** (r + 1) / (math.factorial(r + 1) * math.prod(l - 1 for l in shape))


def random_model_log_height(n: int, shapes, probs_cache=None) -> float:
    """Natural-log height at which the expected number of prime-covering GAPs
    (summed over the given shapes) equals 1, assuming independence across primes."""
    probs_cache = {} if probs_cache is None else probs_cache
    terms = []
    for shape in shapes:
        lp = 0.0
        for p in sympy.primerange(2, n + 1):
            key = (p, tuple(shape))
            if key not in probs_cache:
                probs_cache[key] = residue_cover_prob_np(p, tuple(shape))
            P = probs_cache[key]
            if P <= 0:
                lp = -math.inf
                break
            lp += math.log(P)
        r = len(shape)
        # log E = (r+1) log H - log((r+1)! prod(L-1)) + lp
        const = -math.log(math.factorial(r + 1) * math.prod(l - 1 for l in shape)) + lp
        terms.append((r + 1, const))
    # solve sum_i exp(a_i h + const_i) = 1 for h = log H by bisection
    def logE(h):
        vals = [a * h + cst for a, cst in terms if cst > -math.inf]
        mx = max(vals)
        return mx + math.log(sum(math.exp(v - mx) for v in vals))
    lo, hi = 0.0, 1.0
    while logE(hi) < 0:
        hi *= 2
    for _ in range(100):
        mid = (lo + hi) / 2
        if logE(mid) < 0:
            lo = mid
        else:
            hi = mid
    return hi


def shapes_of_rank(r: int, Lbox: int, exact: bool = True):
    """All non-increasing shapes (L_1 >= ... >= L_r >= 2) with prod == Lbox (exact) or <= Lbox."""
    out = []

    def rec(prefix, remaining_rank, maxside, prod):
        if remaining_rank == 0:
            if (prod == Lbox) or (not exact and prod <= Lbox):
                out.append(tuple(prefix))
            return
        for l in range(min(maxside, Lbox // prod), 1, -1):
            rec(prefix + [l], remaining_rank - 1, l, prod * l)

    rec([], r, Lbox, 1)
    return out


def residue_cover_prob_np(p: int, shape: tuple) -> float:
    """Vectorized exact P_p (same value as residue_cover_prob)."""
    import numpy as np
    r = len(shape)
    grids = np.meshgrid(*[np.arange(p, dtype=np.int64)] * r, indexing="ij")
    C = np.stack([g.ravel() for g in grids], axis=1)  # (p^r, r)
    X = np.array(list(itertools.product(*[range(l) for l in shape])), dtype=np.int64)  # (box, r)
    total = 0
    chunk = max(1, 4_000_000 // max(1, len(X)))
    for s in range(0, len(C), chunk):
        V = (C[s:s + chunk] @ X.T) % p  # (chunk, box)
        V.sort(axis=1)
        distinct = 1 + (np.diff(V, axis=1) != 0).sum(axis=1)
        total += int(distinct.sum())
    return total / (p ** r * p)


# ---------------------------------------------------------------- incidence lemma (Prop. HS-r)
def incidence_params(elements, P):
    """For a finite set of nonzero integers and a set P of primes, return
    (v, Delta, mu, pairs_covered, blocks_proper) where blocks B_e = {p in P : p | e},
    Delta = max_p #{e : p | e}, mu = max_{e != e'} |B_e cap B_e'|."""
    els = sorted({e for e in elements if e != 0})
    P = sorted(set(P))
    blocks = [frozenset(p for p in P if e % p == 0) for e in els]
    v = len(P)
    Delta = max((sum(1 for B in blocks if p in B) for p in P), default=0)
    mu = 0
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            mu = max(mu, len(blocks[i] & blocks[j]))
    pairs_ok = all(any(p in B and q in B for B in blocks)
                   for p, q in itertools.combinations(P, 2))
    proper = all(len(B) < v for B in blocks)
    return v, Delta, mu, pairs_ok, proper


def hs_r_bound_holds(elements, P) -> bool:
    """Prop. HS-r: if every pair of P lies in a block and all blocks are proper, then
    v <= 1 + Delta*(mu*Delta - 1). Returns True when the hypotheses fail (vacuous)."""
    v, D, mu, pairs_ok, proper = incidence_params(elements, P)
    if v < 2 or not pairs_ok or not proper:
        return True
    return v <= 1 + D * (mu * D - 1)
