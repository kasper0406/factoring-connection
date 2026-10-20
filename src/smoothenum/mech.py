"""Enumeration mechanisms for (a, b) pairs with smooth norms, with explicit cost counters (agent smoothenum).

Cost unit: one *examined point* (a lattice/box point whose norm is sieved or tested) or one *lattice setup*
(a modulus M whose lattice is reduced and scanned, even if it holds no box point). Sieve updates are counted
separately. Memory: rows are processed one at a time (O(row length) bytes); nothing materialises the box.

Mechanisms
  box_line_sieve  : plain line sieve of the box |a| <= A, 1 <= b <= Bm (rational, or both sides)
  lattice_points  : exact enumeration of {(a, b) in box : a == R b (mod M)} via a reduced basis
  lattice_sample  : uniform sample of such points (rejection in lattice coordinates) for yield-per-point
  bottom_up_hit   : given an integer s, the unique (a, b) with G(a, b) = s and |a| <= |Y0|/2 (if in box)
Floating point appears only in the approximate log-sieve threshold (as in every real siever); every reported
relation is confirmed by the exact SideTest.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from math import gcd, log2
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np

from src.smoothenum.nfsdata import FactorBase, PolyPair, SideTest

SCALE = 2.0  # log units per bit in the uint8 sieve array


@dataclass
class Counters:
    points: int = 0          # examined points (sieved or tested)
    setups: int = 0          # lattices set up / scanned
    loop_steps: int = 0      # j-steps of lattice scans (bookkeeping overhead)
    updates: int = 0         # sieve updates
    tests: int = 0           # exact side tests performed
    survivors: int = 0
    found: int = 0

    def cost(self) -> int:
        return self.points + self.setups

    def as_dict(self) -> Dict[str, int]:
        d = dict(self.__dict__)
        d["cost"] = self.cost()
        return d


# --------------------------------------------------------------------------- line sieve (box)

def _sieve_row(acc: np.ndarray, fb: FactorBase, lp: np.ndarray, b: int, a_lo: int, cnt: Counters) -> None:
    """acc[i] accumulates log p over p | norm(a_lo + i, b), i.e. a == r b (mod p)."""
    W = acc.shape[0]
    starts = (fb.r * b - a_lo) % fb.p
    for p, s, l in zip(fb.p.tolist(), starts.tolist(), lp.tolist()):
        if s < W:
            sl = acc[s::p]
            sl += l
            cnt.updates += sl.shape[0]


def box_line_sieve(P: PolyPair, A: int, b_lo: int, b_hi: int, rfb: FactorBase, rtest: SideTest,
                   afb: Optional[FactorBase] = None, atest: Optional[SideTest] = None,
                   slack_bits: float = 4.0, a_window: Optional[Tuple[int, int]] = None,
                   ) -> Tuple[List[Tuple[int, int]], Counters]:
    """Line-sieve rows b_lo <= b < b_hi over a in [-A, A) (or a_window = [a_lo, a_hi)). Rational side only if
    afb is None, else both sides. Threshold: accumulated log >= log2|norm| - (mfb + slack_bits).
    Survivors with gcd(a, b) = 1 get the exact tests. One row in memory at a time."""
    cnt = Counters()
    a_lo, a_hi = a_window if a_window is not None else (-A, A)
    W = a_hi - a_lo
    rlp = np.rint(np.log2(rfb.p.astype(np.float64)) * SCALE).astype(np.uint8)
    alp = None if afb is None else np.rint(np.log2(afb.p.astype(np.float64)) * SCALE).astype(np.uint8)
    av = np.arange(a_lo, a_hi, dtype=np.float64)
    out = []
    racc = np.zeros(W, np.uint8)
    aacc = np.zeros(W, np.uint8) if afb is not None else None
    for b in range(b_lo, b_hi):
        racc[:] = 0
        _sieve_row(racc, rfb, rlp, b, a_lo, cnt)
        gl = np.log2(np.abs(P.Y1 * av + P.Y0 * float(b)) + 1.0)
        mask = racc.astype(np.float64) >= (gl - (rtest.mfb + slack_bits)) * SCALE
        if afb is not None:
            aacc[:] = 0
            _sieve_row(aacc, afb, alp, b, a_lo, cnt)
            # homogeneous F(a, b) = sum c_i a^i b^(d-i); float only for the threshold
            fv = sum(float(c) * av ** i * float(b) ** (P.d - i) for i, c in enumerate(P.f))
            fl = np.log2(np.abs(fv) + 1.0)
            mask &= aacc.astype(np.float64) >= (fl - (atest.mfb + slack_bits)) * SCALE
        cnt.points += W
        idx = np.nonzero(mask)[0]
        cnt.survivors += idx.shape[0]
        for i in idx.tolist():
            a = a_lo + i
            if gcd(a, b) != 1:
                continue
            cnt.tests += 1
            if not rtest.ok(P.G(a, b)):
                continue
            if afb is not None:
                cnt.tests += 1
                if not atest.ok(P.F(a, b)):
                    continue
            out.append((a, b))
    cnt.found = len(out)
    return out, cnt


# --------------------------------------------------------------------------- lattices {a == R b mod M}

def reduced_basis(R: int, M: int, A: int, Bm: int) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Lagrange-Gauss reduction of {(M, 0), (R, 1)} in the skewed norm a^2 Bm^2 + b^2 A^2 (exact ints)."""
    def nrm(v):
        return v[0] * v[0] * Bm * Bm + v[1] * v[1] * A * A

    def dot(u, v):
        return u[0] * v[0] * Bm * Bm + u[1] * v[1] * A * A
    u, v = (M, 0), (R % M, 1)
    if nrm(u) < nrm(v):
        u, v = v, u
    while True:
        # make v the shorter
        nv = nrm(v)
        if nv == 0:
            break
        k = (2 * dot(u, v) + nv) // (2 * nv)  # nearest integer to dot/nv
        u = (u[0] - k * v[0], u[1] - k * v[1])
        if nrm(u) >= nv:
            break
        u, v = v, u
    return v, u


def _ceil_div(x: int, y: int) -> int:
    return -((-x) // y)


def _i_range(ca: int, cb: int, ja: int, jb: int, A: int, Bm: int) -> Tuple[int, int]:
    """Integers i with |i ca + ja| <= A and 1 <= i cb + jb <= Bm."""
    lo, hi = -(10 ** 30), 10 ** 30
    for c, k, L, H in ((ca, ja, -A, A), (cb, jb, 1, Bm)):
        if c == 0:
            if not (L <= k <= H):
                return 1, 0
            continue
        if c > 0:
            lo = max(lo, _ceil_div(L - k, c))
            hi = min(hi, (H - k) // c)
        else:
            lo = max(lo, _ceil_div(H - k, c))
            hi = min(hi, (L - k) // c)
    return lo, hi


def _j_bounds(u, v, A: int, Bm: int) -> Tuple[int, int]:
    D = u[0] * v[1] - u[1] * v[0]
    js = []
    for a, b in ((-A, 1), (-A, Bm), (A, 1), (A, Bm)):
        # point = i u + j v  =>  j = (u_a b - u_b a) / D
        js.append((u[0] * b - u[1] * a) / D)
    return int(np.floor(min(js))) - 1, int(np.ceil(max(js))) + 1


def lattice_points(R: int, M: int, A: int, Bm: int, cnt: Counters) -> Iterator[Tuple[int, int]]:
    """All (a, b) in the box with a == R b (mod M); counts one setup + j-steps + points."""
    cnt.setups += 1
    u, v = reduced_basis(R, M, A, Bm)
    j0, j1 = _j_bounds(u, v, A, Bm)
    for j in range(j0, j1 + 1):
        cnt.loop_steps += 1
        lo, hi = _i_range(u[0], u[1], j * v[0], j * v[1], A, Bm)
        for i in range(lo, hi + 1):
            a, b = i * u[0] + j * v[0], i * u[1] + j * v[1]
            cnt.points += 1
            yield a, b


def lattice_count(R: int, M: int, A: int, Bm: int) -> int:
    u, v = reduced_basis(R, M, A, Bm)
    j0, j1 = _j_bounds(u, v, A, Bm)
    tot = 0
    for j in range(j0, j1 + 1):
        lo, hi = _i_range(u[0], u[1], j * v[0], j * v[1], A, Bm)
        tot += max(0, hi - lo + 1)
    return tot


def lattice_sample(R: int, M: int, A: int, Bm: int, k: int, rng: random.Random) -> List[Tuple[int, int]]:
    """k uniform samples (with replacement) from the box points of the lattice, by choosing a uniform j-row
    weighted by its point count (exact counts), then a uniform i. Returns [] if the lattice has no box point."""
    u, v = reduced_basis(R, M, A, Bm)
    j0, j1 = _j_bounds(u, v, A, Bm)
    rows = []
    for j in range(j0, j1 + 1):
        lo, hi = _i_range(u[0], u[1], j * v[0], j * v[1], A, Bm)
        if hi >= lo:
            rows.append((j, lo, hi))
    if not rows:
        return []
    w = [hi - lo + 1 for _, lo, hi in rows]
    out = []
    for j, lo, hi in rng.choices(rows, weights=w, k=k):
        i = rng.randint(lo, hi)
        out.append((i * u[0] + j * v[0], i * u[1] + j * v[1]))
    return out


def crt_root(roots: List[Tuple[int, int]]) -> Tuple[int, int]:
    """Combine conditions a == r_i b (mod p_i) into a == R b (mod M)."""
    R, M = 0, 1
    for r, p in roots:
        # R' == R mod M, R' == r mod p
        t = ((r - R) * pow(M, -1, p)) % p
        R, M = R + M * t, M * p
    return R % M, M


# --------------------------------------------------------------------------- bottom-up

def bottom_up_hit(P: PolyPair, s: int, A: int, Bm: int) -> Optional[Tuple[int, int]]:
    """Solve Y1 a + Y0 b = s. All solutions are (a0 + k Y0', b0 - k Y1') with Y' = Y/gcd; when |Y0| > 2A there is
    at most one with |a| <= A (the centred residue a == s Y1^{-1} mod |Y0|). Returns it if it lies in the box."""
    g = gcd(P.Y0, P.Y1)
    if s % g:
        return None
    Y0, Y1, s = P.Y0 // g, P.Y1 // g, s // g
    m0 = abs(Y0)
    a = (s * pow(Y1, -1, m0)) % m0
    if a > m0 // 2:
        a -= m0
    if (s - Y1 * a) % Y0:
        return None
    b = (s - Y1 * a) // Y0
    if b < 0:
        a, b = -a, -b  # (a, b) and (-a, -b) give G = +-s; keep b > 0
    if abs(a) <= A and 1 <= b <= Bm:
        return a, b
    return None
