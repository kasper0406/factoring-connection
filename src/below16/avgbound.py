"""HARNESS: exact pair counts for local covering generators and the explicit
bounds of Theorems T1 / T2 (notes/below16_results.md, Section 2).

This module uses p and q (it enumerates the population S_n).  It is never
imported by operational code.

Population.  S_n = {(p, q) : p < q <= 2p primes, 2^{n-1} <= pq < 2^n}.  Then
p in [P/2, P), q in (P/sqrt2, sqrt2 P) with P = 2^{n/2}.

T1 generator (residue-local, arbitrary finite windows)
------------------------------------------------------
Fixed integers A, B (not both 0), M >= 1, K with M | K, D >= 1, and a side
('p' or 'q').  With x the class variable (x = p for side 'p', x = q for side
'q') and y the partner, the generator covers (p, q) iff

    x mod M  in R(kappa)   and   floor((B x + A y) / D)  in S(kappa),
    kappa = pq mod K.

R and S are arbitrary maps from keys to finite sets (the design may choose
them adversarially; validity of the class lock is irrelevant for the count).

T2 generator (size-local, interval windows)
-------------------------------------------
Side 'p', A, B >= 1, M | K, D >= 1, Delta = 2^j >= P*K (boxes partition
[2^{n-1}, 2^n) exactly).  Key (kappa, b) = (N mod K, (N - 2^{n-1}) // Delta);
covers iff p mod M in R(key) and floor((Bp+Aq)/D) in [U, U+W) (key-dependent).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, FrozenSet, Iterable, List, Optional, Tuple

import numpy as np
import gmpy2


# --------------------------------------------------------------------------
# population
# --------------------------------------------------------------------------

def primes_upto(m: int) -> np.ndarray:
    s = np.ones(m + 1, dtype=bool)
    s[:2] = False
    for i in range(2, math.isqrt(m) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def population(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """All (p, q) in S_n as int64 arrays (exhaustive; n <= 36 sensible)."""
    lo, hi = 1 << (n - 1), 1 << n
    pr = primes_upto(1 << ((n + 3) // 2))
    ps, qs = [], []
    for i, p in enumerate(pr):
        p = int(p)
        if p * p >= hi:
            break
        if 2 * p * p < lo:
            continue
        qlo = max(p + 1, -(-lo // p))
        qhi = min(2 * p, (hi - 1) // p)
        if qlo > qhi:
            continue
        a = np.searchsorted(pr, qlo, side="left")
        b = np.searchsorted(pr, qhi, side="right")
        qq = pr[a:b]
        ps.append(np.full(len(qq), p, dtype=np.int64))
        qs.append(qq)
    return np.concatenate(ps), np.concatenate(qs)


def class_var_range(n: int, side: str) -> Tuple[int, int, int]:
    """(xmin, xmax, partner_len_max): integer range of the class variable over
    S_n and an upper bound for the length of the partner's range."""
    lo, hi = 1 << (n - 1), 1 << n
    if side == "p":
        xmin = math.isqrt((lo - 1) // 2) + 1          # 2p^2 >= lo
        xmax = math.isqrt(hi - 1)                       # p^2 < hi
        return xmin, xmax, xmax                        # q in (p, 2p]: length <= p
    xmin = math.isqrt(lo - 1) + 1                       # q^2 > pq >= lo
    xmax = math.isqrt(2 * (hi - 1))                     # q^2 <= 2pq < 2hi
    return xmin, xmax, xmax // 2 + 1                   # p in [q/2, q)


def count_in_class(lo: int, hi: int, r: int, m: int) -> int:
    """#{x in [lo, hi] : x == r (mod m)}."""
    if hi < lo:
        return 0
    return (hi - r) // m - (lo - 1 - r) // m


def n_count_in_class(lo: int, hi: int, r: int, m: int) -> int:
    """#{N in [lo, hi) : N == r (mod m)}."""
    return count_in_class(lo, hi - 1, r, m)


# --------------------------------------------------------------------------
# T1
# --------------------------------------------------------------------------

@dataclass
class GenT1:
    A: int
    B: int
    M: int
    K: int
    D: int
    side: str = "p"
    R: Dict[int, FrozenSet[int]] = field(default_factory=dict)
    S: Dict[int, FrozenSet[int]] = field(default_factory=dict)

    def check(self, n: int) -> None:
        assert self.K % self.M == 0 and self.D >= 1 and (self.A, self.B) != (0, 0)
        xmin, _, _ = class_var_range(n, "p")
        assert self.K < xmin, "need K < P/2 so that both primes are units mod K"


def count_T1(g: GenT1, p: np.ndarray, q: np.ndarray) -> int:
    x, y = (p, q) if g.side == "p" else (q, p)
    kap = (p * q) % g.K
    cls = x % g.M
    val = (g.B * x + g.A * y) // g.D
    keys = np.array(sorted(set(g.R) & set(g.S)), dtype=np.int64)
    if len(keys) == 0:
        return 0
    # R membership via table
    rt = np.zeros((g.K, g.M), dtype=bool)
    for k, rs in g.R.items():
        for r in rs:
            rt[k, r % g.M] = True
    ok = rt[kap, cls]
    # S membership via combined codes
    vmin = int(val.min()) - 1
    span = int(val.max()) - vmin + 2
    codes = []
    for k, ss in g.S.items():
        for s in ss:
            if vmin <= s < vmin + span:
                codes.append(k * span + (s - vmin))
    codes = np.array(sorted(codes), dtype=np.int64)
    c = kap * span + (val - vmin)
    ok &= np.isin(c, codes, assume_unique=False)
    return int(ok.sum())


def bound_T1_core(g: GenT1, n: int) -> int:
    """B1 of Lemma 1 (exact integer): sum over keys of #class-variable integers
    times the per-(x, key) partner bound."""
    xmin, xmax, plen = class_var_range(n, g.side)
    tot = 0
    for k in set(g.R) & set(g.S):
        if math.gcd(k, g.K) != 1:
            continue
        s = len(g.S[k])
        if s == 0 or len(g.R[k]) == 0:
            continue
        ypart = plen // g.K + 1
        if g.A != 0:
            nx = sum(count_in_class(xmin, xmax, r, g.M) for r in g.R[k])
            phi = min(s * (g.D // (abs(g.A) * g.K) + 1), ypart)
            tot += nx * phi
        else:
            # value depends on x only: count x exactly
            nx = 0
            for r in g.R[k]:
                for v in g.S[k]:
                    # x == r (M) and floor(Bx/D) == v
                    if g.B > 0:
                        lo = -(-v * g.D // g.B)
                        hi = -(-(v + 1) * g.D // g.B) - 1
                    else:
                        lo = (v + 1) * g.D // g.B + 1
                        hi = v * g.D // g.B
                    nx += count_in_class(max(lo, xmin), min(hi, xmax), r, g.M)
            tot += nx * ypart
    return tot


def expect_over_N(n: int, K: int, f: Dict[int, float]) -> Fraction:
    """E_N[f(N mod K)] for N uniform on [2^{n-1}, 2^n) (exact rational)."""
    lo, hi = 1 << (n - 1), 1 << n
    s = Fraction(0)
    for k, v in f.items():
        s += Fraction(v) * n_count_in_class(lo, hi, k, K)
    return s / (hi - lo)


def bound_T1_closed(g: GenT1, n: int) -> float:
    """Closed form CF1 = 2 c_side (P/M) sqrt(2P(D/|A|+K)) E_N[|R| sqrt|S|]
    (A != 0) or 4 c_side (P/M) sqrt(P(D/|B|+M)) E_N[|R| sqrt|S|] (A == 0).
    Float only as the final comparison value (sqrt)."""
    P = 2.0 ** (n / 2)
    c_side = 1.0 if g.side == "p" else 1.25
    f = {k: len(g.R[k]) * math.sqrt(len(g.S[k])) for k in set(g.R) & set(g.S)}
    E = float(expect_over_N(n, g.K, f))
    if g.A != 0:
        return 2 * c_side * (P / g.M) * math.sqrt(2 * P * (g.D / abs(g.A) + g.K)) * E
    return 4 * c_side * (P / g.M) * math.sqrt(P * (g.D / abs(g.B) + g.M)) * E


def t1_conditions(g: GenT1, n: int) -> bool:
    P = 2.0 ** (n / 2)
    return g.M <= P / 2 and g.K <= P / 2 and g.K <= 2 ** (n - 2)


# --------------------------------------------------------------------------
# T2
# --------------------------------------------------------------------------

@dataclass
class GenT2:
    A: int
    B: int
    M: int
    K: int
    D: int
    Delta: int
    R: Dict[Tuple[int, int], FrozenSet[int]] = field(default_factory=dict)
    UW: Dict[Tuple[int, int], Tuple[int, int]] = field(default_factory=dict)


def t2_delta(n: int, K: int) -> int:
    """Smallest power of two >= P*K (P = 2^{n/2}); it divides 2^{n-1}."""
    j = math.ceil(n / 2 + math.log2(K) - 1e-12)
    while (1 << j) < math.isqrt(1 << n) * K:
        j += 1
    return 1 << j


def keys_T2(g: GenT2, n: int, p: np.ndarray, q: np.ndarray):
    N = p * q
    return N % g.K, (N - (1 << (n - 1))) // g.Delta


def count_T2(g: GenT2, n: int, p: np.ndarray, q: np.ndarray) -> int:
    kap, b = keys_T2(g, n, p, q)
    val = (g.B * p + g.A * q) // g.D
    tot = 0
    # group by key present in design
    want = set(g.R) & set(g.UW)
    if not want:
        return 0
    code = kap * (int(b.max()) + 1) + b
    nb = int(b.max()) + 1
    for (k, bb) in want:
        m = code == k * nb + bb
        if not m.any():
            continue
        U, W = g.UW[(k, bb)]
        rs = np.array(sorted(g.R[(k, bb)]), dtype=np.int64)
        sel = m & np.isin(p % g.M, rs) & (val >= U) & (val < U + W)
        tot += int(sel.sum())
    return tot


def cplus_T2(g: GenT2, n: int, key: Tuple[int, int]) -> int:
    """|C^+(key)|: integers p (class var) in the p-range with p mod M in R and
    some integer q == kappa/p (mod K) in the box with floor((Bp+Aq)/D) in window."""
    k, bb = key
    lo_n = (1 << (n - 1)) + bb * g.Delta
    hi_n = min(lo_n + g.Delta, 1 << n)
    U, W = g.UW[key]
    xmin, xmax, _ = class_var_range(n, "p")
    cnt = 0
    for r in g.R[key]:
        start = xmin + ((r - xmin) % g.M)
        for x in range(start, xmax + 1, g.M):
            if math.gcd(x, g.K) != 1:
                continue
            a = (k * pow(x, -1, g.K)) % g.K
            # q range: q in (x, 2x], lo_n <= xq < hi_n, D*U <= Bx + Aq < D*(U+W)
            qlo = max(x + 1, -(-lo_n // x), -(-(g.D * U - g.B * x) // g.A))
            qhi = min(2 * x, (hi_n - 1) // x, (g.D * (U + W) - g.B * x - 1) // g.A)
            if qhi < qlo:
                continue
            if count_in_class(qlo, qhi, a, g.K) > 0:
                cnt += 1
    return cnt


def bound_T2_core(g: GenT2, n: int) -> int:
    """B2 = sum_keys |C^+(key)| * (floor((Delta/xmin)/K) + 2)."""
    xmin, _, _ = class_var_range(n, "p")
    per = (g.Delta // xmin + 1) // g.K + 1
    return sum(cplus_T2(g, n, key) * per for key in set(g.R) & set(g.UW))


def cplus_closed(g: GenT2, n: int, key) -> float:
    P = 2.0 ** (n / 2)
    U, W = g.UW[key]
    I = g.D * W + 2 * g.A * g.Delta / P
    return len(g.R[key]) * (2 * math.sqrt(2 * I * P / g.A) / g.M + 2)


def bound_T2_closed(g: GenT2, n: int) -> float:
    """CF2 = (3 * 2^n / P) * E_N[ |R| (2 sqrt(2 I P / A)/M + 2) ]."""
    P = 2.0 ** (n / 2)
    lo, hi = 1 << (n - 1), 1 << n
    tot = 0.0
    for key in set(g.R) & set(g.UW):
        k, bb = key
        a = lo + bb * g.Delta
        bnd = min(a + g.Delta, hi)
        w = n_count_in_class(a, bnd, k, g.K) / (hi - lo)
        tot += w * cplus_closed(g, n, key)
    return 3 * (2 ** n) / P * tot


# --------------------------------------------------------------------------
# designs
# --------------------------------------------------------------------------

def sqrt_mod_prime(a: int, l: int) -> Optional[int]:
    a %= l
    if a == 0:
        return 0
    if pow(a, (l - 1) // 2, l) != 1:
        return None
    # Tonelli-Shanks
    qq, s = l - 1, 0
    while qq % 2 == 0:
        qq //= 2
        s += 1
    z = 2
    while pow(z, (l - 1) // 2, l) != l - 1:
        z += 1
    m, c, t, r = s, pow(z, qq, l), pow(a, qq, l), pow(a, (qq + 1) // 2, l)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % l
            i += 1
        b = pow(c, 1 << (m - i - 1), l)
        m, c, t, r = i, b * b % l, t * b * b % l, r * b % l
    return r


def dcl_generators(n: int, H: int, Y: int, K_sq: bool = False) -> List[GenT1]:
    """DCL as T1 generators: (A,B) coprime <= H, A != B, primes l in (Y, 2Y],
    R(kappa) = roots of B r^2 == A kappa (mod l), n-only whole-range window
    (a superset of DCL-det's window).  K = l (or l^2 if K_sq)."""
    lo, hi = 1 << (n - 1), 1 << n
    xmin, xmax, _ = class_var_range(n, "p")
    gens = []
    for l in primes_upto(2 * Y):
        l = int(l)
        if l <= Y or l <= 2:
            continue
        for A in range(1, H + 1):
            for B in range(1, H + 1):
                if A == B or math.gcd(A, B) != 1 or l <= max(A, B):
                    continue
                K = l * l if K_sq else l
                D = l * l
                # n-only window: values floor((Bp+Aq)/l^2) over p in [xmin,xmax], q in (p,2p]
                Fmin = 2 * math.isqrt(A * B * lo)        # AM-GM lower bound (floor)
                Fmax = max(B * xmin + A * (2 * xmin + 2), B * xmax + A * (hi // xmax + 1),
                           B * xmin + A * (hi // xmin + 1), B * xmax + A * (2 * xmax))
                S = frozenset(range(Fmin // D - 1, Fmax // D + 2))
                R, SS = {}, {}
                for kap in range(K):
                    if math.gcd(kap, K) != 1:
                        continue
                    r = sqrt_mod_prime(A * kap * pow(B, -1, l), l)
                    if r is None or r == 0:
                        continue
                    R[kap] = frozenset({r, (l - r) % l})
                    SS[kap] = S
                gens.append(GenT1(A, B, l, K, D, "p", R, SS))
    return gens
