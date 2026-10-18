"""Continued-fraction PREFIX features of sqrt(N).  OPERATIONAL: input is N only.

Recurrence (N not a square):  P_0 = 0, Q_0 = 1, a_0 = isqrt(N);
    P_{i+1} = a_i Q_i - P_i,  Q_{i+1} = (N - P_{i+1}^2) / Q_i,  a_{i+1} = (a_0 + P_{i+1}) // Q_{i+1}.
Convergents A_{-1} = 1, A_0 = a_0, A_i = a_i A_{i-1} + A_{i-2} (same for B with B_{-1}=0, B_0=1) satisfy
    A_{i-1}^2 - N B_{i-1}^2 = (-1)^i Q_i      (i >= 1),
so (-1)^i Q_i is a square modulo N for every i >= 1 (tested in tests/test_eps.py).

One-sided signals (proved in notes/eps_prereg.md section 3; never consult factors):
  PC1  period l (least i >= 1 with Q_i = 1) is odd and l <= T  => x^2 - N y^2 = -1 solvable => eps = +1.
  PC2  some odd i <= T has Q_i = k^2                            => -k^2 square mod N, gcd(k,N)=1 => eps = +1.
  PC3  e_{-1} in F_2-span of exponent vectors of the B-smooth relations (-1)^i Q_i, i <= T  => eps = +1.
  PCQ2 some even i <= T has Q_i = 2 k^2                         => 2 square mod p and q => QR_2 = 1.
  PCQ3 e_{2} in that span                                        => QR_2 = 1.
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List

import gmpy2

T_LIST = (64, 1024)
T_MAX = 1024
FB_BOUND = 256
FB_PRIMES: List[int] = [x for x in range(2, FB_BOUND + 1) if all(x % d for d in range(2, int(x ** 0.5) + 1))]
_PRIMORIAL = 1
for _l in FB_PRIMES:
    _PRIMORIAL *= _l
# vector coordinates: bit 0 = sign (-1), bit 1+j = FB_PRIMES[j]  (bit 1 is the prime 2)
BIT_M1 = 1
BIT_2 = 1 << 1


def _is_sq(x: int) -> bool:
    return bool(gmpy2.is_square(x))


def _smooth_vector(sign_odd: bool, Q: int):
    """Exponent-parity vector of (-1)^i Q over {-1} u FB, or None if Q is not FB-smooth."""
    m = Q
    while True:
        g = gcd(m, _PRIMORIAL)
        if g == 1:
            break
        m //= g
    if m != 1:
        return None
    v = BIT_M1 if sign_odd else 0
    m = Q
    for j, l in enumerate(FB_PRIMES):
        if m % l == 0:
            e = 0
            while m % l == 0:
                m //= l
                e += 1
            if e & 1:
                v ^= 1 << (j + 1)
            if m == 1:
                break
    return v


class _Span:
    """Incremental F_2 row-echelon basis on int bitmasks (pivot = highest set bit)."""

    def __init__(self):
        self.basis: Dict[int, int] = {}

    def reduce(self, v: int) -> int:
        while v:
            h = v.bit_length() - 1
            b = self.basis.get(h)
            if b is None:
                return v
            v ^= b
        return 0

    def add(self, v: int) -> None:
        v = self.reduce(v)
        if v:
            self.basis[v.bit_length() - 1] = v

    def contains(self, v: int) -> bool:
        return self.reduce(v) == 0


def cf_features(N: int, T_max: int = T_MAX, T_list=T_LIST) -> Dict[str, int]:
    """Menu statistics + positive-control indicators from the first T_max steps of the CF of sqrt(N).

    Returns (all ints):
      cf{T}_n_a1, cf{T}_n_a2      #{1<=i<=T : a_i = 1}, #{... a_i = 2}
      cf{T}_sumbits, cf{T}_maxbits  sum / max of bit_length(a_i), 1<=i<=T
      cf{T}_n_Qeven               #{1<=i<=T : Q_i even}
      cf{T}_n_sq_even             #{even i <= T : Q_i a perfect square}
      cf{T}_n_Q1mod4_odd          #{odd i <= T : Q_i = 1 mod 4}
      cf{T}_n_2sq_odd             #{odd i <= T : Q_i = 2 k^2}
      pc1_T{T}, pc2_T{T}, pc3_T{T}, pcq2_T{T}, pcq3_T{T}   indicators in {0,1}
      first_pc1, first_pc2, first_pc3, first_pcq2, first_pcq3: first index i <= T_max at which the
      signal holds (0 if never), for frequency-vs-T curves.
      cf_nsmooth: number of FB-smooth relations found (smoothness is tested only until both PC3 and
      PCQ3 have fired, so this is an auxiliary count, not a menu statistic).
      cf_nonunit: square signals rejected because gcd(k, N) > 1 (never for balanced N);
      cf_fb_ok: 1 iff gcd(N, prod of FB primes) = 1 (PC3/PCQ3 are disabled otherwise).
    """
    a0 = int(gmpy2.isqrt(N))
    if a0 * a0 == N:
        raise ValueError("N is a perfect square")
    P, Q = 0, 1
    a = a0
    acc = {T: dict(n_a1=0, n_a2=0, sumbits=0, maxbits=0, n_Qeven=0, n_sq_even=0, n_Q1mod4_odd=0,
                   n_2sq_odd=0) for T in T_list}
    first = dict(pc1=0, pc2=0, pc3=0, pcq2=0, pcq3=0, Q1=0)
    span = _Span()
    nsmooth = 0
    nonunit = 0
    # PC3 needs every factor-base prime coprime to N (else a relation could involve p itself)
    fb_ok = gcd(N, _PRIMORIAL) == 1
    Ts = sorted(T_list)
    for i in range(1, T_max + 1):
        P = a * Q - P
        Q = (N - P * P) // Q
        a = (a0 + P) // Q
        odd = i & 1
        # per-step facts
        is_sq = _is_sq(Q)
        is_2sq = (Q & 1) == 0 and _is_sq(Q >> 1)
        ab = a.bit_length()
        if Q == 1 and first["Q1"] == 0:
            first["Q1"] = i                 # period length (any parity), if <= T_max
        if Q == 1 and odd and first["pc1"] == 0:
            first["pc1"] = i
        if odd and is_sq and first["pc2"] == 0:
            # Q = k^2; k <= (2 sqrt N)^(1/2) < p for balanced N, but check gcd(k, N) anyway (rule 3)
            if gcd(int(gmpy2.isqrt(Q)), N) == 1:
                first["pc2"] = i
            else:
                nonunit += 1
        if (not odd) and is_2sq and first["pcq2"] == 0:
            if gcd(int(gmpy2.isqrt(Q >> 1)), N) == 1:
                first["pcq2"] = i
            else:
                nonunit += 1
        if fb_ok and (first["pc3"] == 0 or first["pcq3"] == 0):
            v = _smooth_vector(bool(odd), Q)
            if v is not None:
                nsmooth += 1
                span.add(v)
                if first["pc3"] == 0 and span.contains(BIT_M1):
                    first["pc3"] = i
                if first["pcq3"] == 0 and span.contains(BIT_2):
                    first["pcq3"] = i
        for T in Ts:
            if i > T:
                continue
            d = acc[T]
            if a == 1:
                d["n_a1"] += 1
            elif a == 2:
                d["n_a2"] += 1
            d["sumbits"] += ab
            if ab > d["maxbits"]:
                d["maxbits"] = ab
            if (Q & 1) == 0:
                d["n_Qeven"] += 1
            if (not odd) and is_sq:
                d["n_sq_even"] += 1
            if odd and (Q & 3) == 1:
                d["n_Q1mod4_odd"] += 1
            if odd and is_2sq:
                d["n_2sq_odd"] += 1
    out: Dict[str, int] = {}
    for T in Ts:
        for k, v in acc[T].items():
            out[f"cf{T}_{k}"] = v
        for k in ("pc1", "pc2", "pc3", "pcq2", "pcq3"):
            out[f"{k}_T{T}"] = int(0 < first[k] <= T)
    for k, v in first.items():
        out[f"first_{k}"] = v
    out["cf_nsmooth"] = nsmooth
    out["cf_nonunit"] = nonunit
    out["cf_fb_ok"] = int(fb_ok)
    return out


def cf_sequence(N: int, T: int):
    """(a_i, P_i, Q_i) for i = 0..T (test helper; same recurrence as cf_features)."""
    a0 = int(gmpy2.isqrt(N))
    P, Q, a = 0, 1, a0
    out = [(a, P, Q)]
    for _ in range(T):
        P = a * Q - P
        Q = (N - P * P) // Q
        a = (a0 + P) // Q
        out.append((a, P, Q))
    return out
