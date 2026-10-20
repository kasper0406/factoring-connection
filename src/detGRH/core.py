"""detGRH core (OPERATIONAL: sees only N and public parameters; never uses factors).

Rational index calculus modulo N with a deterministic relation supply, and two ways to extract a factor:
  * the F2 route (square dependencies, Lemma H's functional Phi), and
  * the Z route (full-rank kernel lattice -> multiple of the exponent of G_B -> Miller's test, Prop Z).

Group: Gamma = Z^n with coordinates FB = [-1, 2, 3, 5, ..., primes <= B]; ev(e) = prod q^{e_q} mod N.
A *kernel relation* is e in Z^n with ev(e) = 1.  -1 is treated as a free generator; the relation 2*e_{-1}
is always added, so Gamma/Lambda_c is finite iff Lambda_c has full rank.

Every inversion mod N is gcd-checked; a nonunit found anywhere is returned as a factor exit
(rule 3 of notes/AGENT_RULES.md).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

import flint


# ----------------------------------------------------------------------------- primes, factor base

def primes_upto(B: int) -> List[int]:
    if B < 2:
        return []
    s = bytearray([1]) * (B + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(B) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
    return [i for i in range(B + 1) if s[i]]


def grh_bound(N: int) -> int:
    """ceil(2 (ln N)^2): Bach's ERH bound for generators of (Z/N)^* (Math. Comp. 55 (1990)).

    Parameter selection only (float log); rounded up generously by the caller."""
    return math.ceil(2.0 * math.log(N) ** 2)


def fb_bound(N: int) -> int:
    """B(n): least multiple of 64 that is >= 2 ln^2 N (pre-registered)."""
    b = grh_bound(N)
    return ((b + 63) // 64) * 64


class FactorExit(Exception):
    def __init__(self, divisor: int, where: str):
        super().__init__(f"factor {divisor} found at {where}")
        self.divisor = divisor
        self.where = where


@dataclass
class FactorBase:
    N: int
    B: int
    primes: List[int]                 # rational primes <= B
    index: Dict[int, int] = field(default_factory=dict)   # prime -> column (column 0 is -1)
    prod: int = 1

    @property
    def n(self) -> int:
        return 1 + len(self.primes)

    @classmethod
    def make(cls, N: int, B: Optional[int] = None) -> "FactorBase":
        B = fb_bound(N) if B is None else B
        ps = primes_upto(B)
        for p in ps:                       # trial division = gcd check of every FB element
            g = gcd(p, N)
            if 1 < g < N:
                raise FactorExit(g, "factor base")
        fb = cls(N=N, B=B, primes=ps)
        fb.index = {p: i + 1 for i, p in enumerate(ps)}
        P = 1
        for p in ps:
            P *= p
        fb.prod = P
        return fb

    def vector(self, v: int) -> Optional[List[int]]:
        """Exponent vector of a nonzero integer v over FB (sign in column 0), or None if not B-smooth."""
        if v == 0:
            return None
        e = [0] * self.n
        if v < 0:
            e[0] = 1
            v = -v
        # fast smoothness filter: strip gcds with the primorial
        w = v
        g = gcd(w, self.prod)
        while g > 1:
            while w % g == 0:
                w //= g
            g = gcd(w, g)
        if w != 1:
            return None
        for p in self.primes:
            if v % p == 0:
                c = 0
                while v % p == 0:
                    v //= p
                    c += 1
                e[self.index[p]] = c
                if v == 1:
                    break
        return e

    def ev(self, e: Sequence[int]) -> int:
        """prod q^{e_q} mod N (negative exponents via gcd-checked inverses)."""
        N = self.N
        x = 1
        if e[0] % 2:
            x = N - 1
        for p, i in self.index.items():
            k = e[i]
            if k:
                x = x * pow(p, k, N) % N      # p is a unit: checked in make()
        return x


# ----------------------------------------------------------------------------- incremental rank mod a prime

RANK_PRIME = (1 << 61) - 1   # fixed, deterministic


class ModRank:
    """Incremental row echelon basis over F_P (P = 2^61 - 1). rank mod P <= rank over Q."""

    def __init__(self, n: int, P: int = RANK_PRIME):
        self.n, self.P = n, P
        self.piv: Dict[int, List[int]] = {}

    def add(self, v: Sequence[int]) -> bool:
        P = self.P
        w = [x % P for x in v]
        for j in range(self.n):
            if w[j]:
                if j in self.piv:
                    r = self.piv[j]
                    c = w[j]
                    w = [(a - c * b) % P for a, b in zip(w, r)]
                else:
                    inv = pow(w[j], -1, P)
                    self.piv[j] = [a * inv % P for a in w]
                    return True
        return False

    @property
    def rank(self) -> int:
        return len(self.piv)


# ----------------------------------------------------------------------------- deterministic relation supply

@dataclass
class Relations:
    fb: FactorBase
    kernel: List[List[int]] = field(default_factory=list)      # kernel relations (ev = 1)
    kernel_src: List[str] = field(default_factory=list)
    squares: List[Tuple[int, int, List[int]]] = field(default_factory=list)  # (r, R = r^2 - N, e(R))
    qrank_modP: int = 0
    tested: Dict[str, int] = field(default_factory=dict)


def collect(N: int, fb: FactorBase, K: int, X: int, T: int, surplus: int = 10) -> Relations:
    """Deterministic enumeration (pre-registered order):
    1. the relation 2*e_{-1};
    2. F_seed: r = ceil(sqrt N) + t, t < T, with r and r^2 - N both B-smooth (square relation with root r,
       kernel relation 2 e(r) - e(r^2 - N));
    3. F_lin: k = 1..K, x = 1..X, x + kN B-smooth -> kernel relation e(x) - e(x + kN).
    Stops once rank mod P = n and `surplus` further kernel relations exist, or at the caps."""
    rel = Relations(fb=fb)
    n = fb.n
    mr = ModRank(n)
    e_m1 = [0] * n
    e_m1[0] = 2
    rel.kernel.append(e_m1)
    rel.kernel_src.append("sign")
    mr.add(e_m1)
    extra = 0
    m = isqrt(N)
    if m * m < N:
        m += 1
    if m * m == N:
        raise FactorExit(m, "perfect square")
    cnt_seed = 0
    for t in range(T):
        r = m + t
        cnt_seed += 1
        er = fb.vector(r)
        if er is None:
            continue
        R = r * r - N
        eR = fb.vector(R)
        if eR is None:
            continue
        rel.squares.append((r, R, eR))
        kv = [2 * a - b for a, b in zip(er, eR)]
        if any(kv):
            rel.kernel.append(kv)
            rel.kernel_src.append("seed")
            if not mr.add(kv) and mr.rank == n:
                extra += 1
    cnt_lin = 0
    done = mr.rank == n and extra >= surplus
    for k in range(1, K + 1):
        if done:
            break
        for x in range(1, X + 1):
            cnt_lin += 1
            ex = fb.vector(x)
            if ex is None:
                continue
            ey = fb.vector(x + k * N)
            if ey is None:
                continue
            kv = [a - b for a, b in zip(ex, ey)]
            rel.kernel.append(kv)
            rel.kernel_src.append("lin")
            added = mr.add(kv)
            if not added and mr.rank == n:
                extra += 1
            if mr.rank == n and extra >= surplus:
                done = True
                break
    rel.qrank_modP = mr.rank
    rel.tested = {"seed": cnt_seed, "lin": cnt_lin}
    return rel


# ----------------------------------------------------------------------------- Z route (Prop Z)

def lattice_index(rows: List[List[int]], n: int) -> Tuple[int, int]:
    """(rank over Q, |Z^n / L|) for the row lattice L; index 0 if rank < n. Exact (flint HNF)."""
    if not rows:
        return 0, 0
    M = flint.fmpz_mat(rows)
    H = M.hnf()
    H_rows = [r for r in H.tolist() if any(int(x) for x in r)]
    rk = len(H_rows)
    if rk < n:
        return rk, 0
    D = 1
    for i in range(n):
        D *= abs(int(H_rows[i][i]))
    return rk, D


def miller(N: int, M: int, a_max: int) -> Tuple[Optional[int], Optional[int]]:
    """Miller's test with a multiple M of the exponent: for bases a = 2..a_max, write M = 2^s t,
    find the last x != 1 in a^t, a^{2t}, ...; if x != -1 then gcd(x - 1, N) is proper.
    Returns (factor, base) or (None, None).  Also gcd-checks every base."""
    if M <= 0:
        return None, None
    s, t = 0, M
    while t % 2 == 0:
        t //= 2
        s += 1
    for a in range(2, a_max + 1):
        g = gcd(a, N)
        if 1 < g < N:
            return g, a
        x = pow(a, t, N)
        if x == 1:
            continue
        for _ in range(s + 1):
            y = x * x % N
            if y == 1:
                break
            x = y
        if x * x % N != 1:
            # M is not a multiple of ord(a): cannot happen under the hypotheses; skip base
            continue
        if x != N - 1:
            g = gcd(x - 1, N)
            if 1 < g < N:
                return g, a
    return None, None


# ----------------------------------------------------------------------------- F2 route (Lemma H)

class F2Basis:
    """Reduced row echelon basis over F2 of a list of integer rows, tracking which input rows make up
    each basis row (bitmask over input row indices).  Rows are stored as int bitmasks over columns."""

    def __init__(self, rows: List[List[int]], n: int):
        self.n = n
        self.rows = rows
        self.piv: Dict[int, Tuple[int, int]] = {}     # pivot col -> (row bits, combo bits)
        self.pivmask = 0
        self.deps: List[int] = []                      # combos of input rows that vanish mod 2
        # Invariant (fully reduced RREF): each pivot row contains exactly one pivot column, its own.
        for i, r in enumerate(rows):
            bits = bits_of(r)
            combo = 1 << i
            t = bits & self.pivmask
            while t:
                j = (t & -t).bit_length() - 1
                pb, pc = self.piv[j]
                bits ^= pb
                combo ^= pc
                t = bits & self.pivmask
            if bits:
                j = (bits & -bits).bit_length() - 1    # a non-pivot column
                for jj, (pb, pc) in list(self.piv.items()):
                    if (pb >> j) & 1:
                        self.piv[jj] = (pb ^ bits, pc ^ combo)
                self.piv[j] = (bits, combo)
                self.pivmask |= 1 << j
            else:
                self.deps.append(combo)

    def express(self, target_bits: int) -> Optional[int]:
        """Unique combo (over input rows) of the RREF basis whose sum is target mod 2; None if not in span."""
        combo = 0
        t = target_bits
        p = target_bits & self.pivmask
        while p:
            j = (p & -p).bit_length() - 1
            pb, pc = self.piv[j]
            t ^= pb
            combo ^= pc
            p &= p - 1
        if t:
            return None
        return combo

    def combo_sum(self, combo: int) -> List[int]:
        s = [0] * self.n
        i = 0
        while combo:
            if combo & 1:
                s = [a + b for a, b in zip(s, self.rows[i])]
            combo >>= 1
            i += 1
        return s


def bits_of(v: Sequence[int]) -> int:
    b = 0
    for j, a in enumerate(v):
        if a & 1:
            b |= 1 << j
    return b


def half_ev(fb: FactorBase, v: Sequence[int]) -> int:
    """ev(v/2) for an even integer vector v."""
    assert all(a % 2 == 0 for a in v)
    return fb.ev([a // 2 for a in v])


def f2_fruitful(fb: FactorBase, rows: List[List[int]]) -> Tuple[int, int, Optional[int]]:
    """Test a basis of the F2 dependency space of kernel relations: (#deps, #fruitful basis vectors, a factor)."""
    N = fb.N
    basis = F2Basis(rows, fb.n)
    nfr, fac = 0, None
    for combo in basis.deps:
        s = basis.combo_sum(combo)
        y = half_ev(fb, s)          # y^2 = ev(s) = 1
        if y not in (1, N - 1):
            g = gcd(y - 1, N)
            if 1 < g < N:
                nfr += 1
                fac = fac or g
    return len(basis.deps), nfr, fac


def canonical_root(fb: FactorBase, basis: F2Basis, eR: Sequence[int]) -> Optional[int]:
    """y(R) := ev((e(R) + lambda)/2) with lambda the canonical completion from the RREF basis of the kernel
    relations (a deterministic function of the certificate e(R) and the fixed basis); None if no completion."""
    combo = basis.express(bits_of(eR))
    if combo is None:
        return None
    lam = basis.combo_sum(combo)
    v = [a + b for a, b in zip(eR, lam)]
    return half_ev(fb, v)
