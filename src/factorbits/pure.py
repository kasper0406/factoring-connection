"""Pure generic-ring decision trees (Theorem PT, notes/factorbits_results.md §3).  DIAGNOSTIC.

A *pure* test is an integer c given by a straight-line program whose constants and shape do not
depend on N.  On input N the tree computes c mod N, branches on [c == 0 mod N], and (instrumented)
reports a split when 1 < gcd(c mod N, N) < N.  For a prime x the test answer is [x | c].

Partition lemma (proved in the note): the tree, run on N = xy (x != y primes), either splits N or
ends in the leaf reached by the "diagonal" path of x, which equals that of y.  So the non-split
behaviour is a partition of the primes into leaf classes S_l, and N splits iff x, y lie in
different classes.

Tests are objects with
    .mod(m)  -> residue of c modulo m   (SLP evaluated mod m; used both mod a prime and mod N)
Trees are nested tuples: ("leaf", label) or ("node", test, child_if_nonzero, child_if_zero).
Classes are diagonal transcripts (leaves of the unrolled tree), never output labels.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt, prod
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

import gmpy2


# ----------------------------------------------------------------------------- pure tests (SLPs)
@dataclass(frozen=True)
class PowTest:
    """c = a^M - 1 (Pollard p-1 shape)."""
    a: int
    M: int

    def mod(self, m: int) -> int:
        return (pow(self.a, self.M, m) - 1) % m


@dataclass(frozen=True)
class LucasTest:
    """c = V_M(P) - 2 with V_0 = 2, V_1 = P, V_{k+1} = P V_k - V_{k-1}  (Williams p+1 shape)."""
    P: int
    M: int

    def mod(self, m: int) -> int:
        # Montgomery-style ladder on (V_k, V_{k+1})
        v0, v1 = 2 % m, self.P % m
        for bit in bin(self.M)[2:]:
            if bit == "1":
                v0, v1 = (v0 * v1 - self.P) % m, (v1 * v1 - 2) % m
            else:
                v0, v1 = (v0 * v0 - 2) % m, (v0 * v1 - self.P) % m
        return (v0 - 2) % m


@dataclass(frozen=True)
class ECMTest:
    """Montgomery curve B y^2 = x^3 + A x^2 + x with A = 4*a24 - 2 and integer constants
    (a24, x0) independent of N; c = Z-coordinate of M * P0 from the x-only ladder in projective
    coordinates (no inversions, so c is an integer polynomial expression in a24, x0)."""
    a24: int
    x0: int
    M: int

    def mod(self, m: int) -> int:
        a24 = self.a24 % m
        XP, ZP = self.x0 % m, 1

        def dbl(X, Z):
            s, d = (X + Z) * (X + Z) % m, (X - Z) * (X - Z) % m
            t = (s - d) % m
            return s * d % m, t * (d + a24 * t) % m

        def add(X1, Z1, X2, Z2):
            u = (X1 - Z1) * (X2 + Z2) % m
            v = (X1 + Z1) * (X2 - Z2) % m
            return ZP * (u + v) * (u + v) % m, XP * (u - v) * (u - v) % m

        R0, R1 = (1, 0), (XP, ZP)
        for bit in bin(self.M)[2:]:
            if bit == "0":
                R1 = add(*R0, *R1)
                R0 = dbl(*R0)
            else:
                R0 = add(*R0, *R1)
                R1 = dbl(*R1)
        return R0[1] % m


@dataclass(frozen=True)
class ProductTest:
    """c = prod_{x in [lo, hi]} f(x) with f(x) = x (interval divisor test, Lipton style) or
    f(x) = x^2 + 1 (vanishes mod p iff sqrt(-1) mod p has a representative in [lo, hi])."""
    lo: int
    hi: int
    kind: str = "x"

    def mod(self, m: int) -> int:
        r = 1
        for x in range(self.lo, self.hi + 1):
            r = r * (x if self.kind == "x" else x * x + 1) % m
        return r


# ----------------------------------------------------------------------------- trees
def leaf(label) -> tuple:
    return ("leaf", label)


def node(test, if_nonzero, if_zero) -> tuple:
    return ("node", test, if_nonzero, if_zero)


def nonadaptive(tests: Sequence) -> tuple:
    """Full binary tree over a fixed test list; leaf label = answer tuple."""
    def build(i: int, ans: tuple) -> tuple:
        if i == len(tests):
            return leaf(ans)
        return node(tests[i], build(i + 1, ans + (0,)), build(i + 1, ans + (1,)))
    return build(0, ())


def run_on_prime(tree: tuple, x: int) -> tuple:
    """Diagonal path of x: the transcript of answers [x | c_v] (1 = zero branch).  Classes are
    transcripts (= leaves of the unrolled tree), not output labels: two leaves may share a label."""
    t, tr = tree, ()
    while t[0] == "node":
        _, test, nz, z = t
        a = 1 if test.mod(x) == 0 else 0
        tr += (a,)
        t = z if a else nz
    return tr


def leaf_label(tree: tuple, transcript: tuple):
    t = tree
    for a in transcript:
        t = t[3] if a else t[2]
    return t[1]


def run_on_N(tree: tuple, N: int):
    """Execute the instrumented tree mod N. Returns ("leaf", transcript) or ("split", divisor)."""
    t, tr = tree, ()
    while t[0] == "node":
        _, test, nz, z = t
        r = test.mod(N)
        g = gcd(r, N)
        if 1 < g < N:
            return ("split", g)
        a = 1 if r == 0 else 0
        tr += (a,)
        t = z if a else nz
    return ("leaf", tr)


def classes(tree: tuple, primes: Iterable[int]) -> Dict[int, tuple]:
    return {x: run_on_prime(tree, x) for x in primes}


# ----------------------------------------------------------------------------- exact quantities
def class_dist(cls: Dict[int, object], primes: Sequence[int]) -> Dict[object, float]:
    """Class probabilities under the uniform law on `primes` (exact rationals as floats for
    reporting; the identities are checked with integer counts in the tests)."""
    out: Dict[object, int] = {}
    for x in primes:
        out[cls[x]] = out.get(cls[x], 0) + 1
    return {k: v / len(primes) for k, v in out.items()}


def split_prob_product(cls: Dict[int, object], P: Sequence[int], Q: Sequence[int]) -> float:
    """P(class(X) != class(Y)) for X ~ U(P), Y ~ U(Q) independent: 1 - sum_l pi_l rho_l."""
    a = class_dist(cls, P)
    b = class_dist(cls, Q)
    return 1.0 - sum(a[k] * b.get(k, 0.0) for k in a)


def pt_bound_check(cls: Dict[int, object], primes: Sequence[int],
                   target: Callable[[int, int], int]) -> Dict[str, float]:
    """iid uniform X, Y on `primes` (ordered pairs, x != y not imposed; x == y has measure 1/|P|
    and is counted as non-split).  Computes the Bayes-optimal success of the instrumented tree for
    predicting target(x, y) (a split counts as success: it factors N), the baseline max_b P(F=b),
    P_split, eta = 1 - max_l pi_l, and the Theorem PT bound baseline + P_split + eta^2."""
    n = len(primes)
    by_class: Dict[object, List[int]] = {}
    for x in primes:
        by_class.setdefault(cls[x], []).append(x)
    tot = n * n
    split = tot - sum(len(v) ** 2 for v in by_class.values())
    # baseline
    cnt: Dict[int, int] = {}
    for x in primes:
        for y in primes:
            b = target(x, y)
            cnt[b] = cnt.get(b, 0) + 1
    baseline = max(cnt.values()) / tot
    # Bayes success on non-split pairs, per class
    good = 0
    for S in by_class.values():
        c: Dict[int, int] = {}
        for x in S:
            for y in S:
                b = target(x, y)
                c[b] = c.get(b, 0) + 1
        good += max(c.values())
    succ = (good + split) / tot
    pi_max = max(len(v) for v in by_class.values()) / n
    eta = 1.0 - pi_max
    ps = split / tot
    return {"success": succ, "baseline": baseline, "p_split": ps, "eta": eta,
            "bound": baseline + ps + eta ** 2, "n_classes": len(by_class)}


def eps_cross_check(cls: Dict[int, object], P1: Sequence[int], P3: Sequence[int]) -> Dict[str, float]:
    """Population A: eps = +1 w.p. 1/2 (X, Y iid U(P1)), eps = -1 w.p. 1/2 (X, Y iid U(P3)).
    Population B: X ~ U(P1), Y ~ U(P3).  Bayes success (split = success) on A, and the bound
    1/2 + P_split(A)/2 + sqrt(2 P_split(B))/2."""
    a = class_dist(cls, P1)
    b = class_dist(cls, P3)
    keys = set(a) | set(b)
    sA1 = 1 - sum(v * v for v in a.values())
    sA3 = 1 - sum(v * v for v in b.values())
    pA = (sA1 + sA3) / 2
    pB = 1 - sum(a.get(k, 0) * b.get(k, 0) for k in keys)
    good = 0.5 * sum(max(a.get(k, 0) ** 2, b.get(k, 0) ** 2) for k in keys)
    succ = good + pA
    return {"success": succ, "p_split_A": pA, "p_split_B": pB,
            "bound": 0.5 + pA / 2 + (2 * pB) ** 0.5 / 2}


# ----------------------------------------------------------------------------- prime ranges
def primes_top2(bits: int) -> List[int]:
    """All primes with exactly `bits` bits and the top two bits set (the balanced_semiprime range)."""
    lo, hi = 3 << (bits - 2), (1 << bits) - 1
    out, x = [], int(gmpy2.next_prime(lo - 1))
    while x <= hi:
        out.append(x)
        x = int(gmpy2.next_prime(x))
    return out


def smooth_exponent(B: int, bound_bits: int) -> int:
    """M_B = prod_{l <= B prime} l^{floor(log_l 2^bound_bits)}."""
    M, l = 1, 2
    while l <= B:
        e = 1
        while l ** (e + 1) < (1 << bound_bits):
            e += 1
        M *= l ** e
        l = int(gmpy2.next_prime(l))
    return M
