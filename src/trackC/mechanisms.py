"""Unit-only mechanism attempts for U(K) / I(A,B) (track C, owner factC).

OPERATIONAL functions (input: N, public interval/config, seed).  No factor
information.  All inversions gcd-checked; nonunits are returned, not hidden.
Diagnostic evaluation lives in src/trackC/harness_mechanisms.py.

M1  random_shift_probe   : nonunit status of prod_{i<L}(s+i) at random shifts s.
M2  power_sum_faulhaber  : S_e(A,B) = sum_{A<=j<B} j^e mod N via Faulhaber /
                           Bernoulli numbers mod N (cost O(e^2) + O(e log B)).
M2' power_sum_direct     : same sum by direct powering (cost O(L log e)); used
                           for huge e (e = N-1, Fermat-type sums).
M3  truncated_FN_test    : Z(N,K) = [(1+X)^N == 1 mod (N, X^K)]  (= U(K)).
M4  jacobi_prefix_product: prod_{j<K} (j/N)  (codex r1 idea 2; zero iff U fails).
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Optional

from src.common.metrics import Metrics
from src.common.ring import FlintQuotient, NonUnit, inv_mod, power_of_linear

from .baselines import bgs_interval


# --------------------------------------------------------------------- M1
def random_shift_probe(N: int, L: int, t: int, seed: int, m: Optional[Metrics] = None) -> List[Dict]:
    """Evaluate the unit status of prod_{i<L}(s+i) mod N at t random shifts s."""
    m = m if m is not None else Metrics()
    rng = random.Random(f"trackC_M1|{N}|{L}|{t}|{seed}")
    out = []
    for _ in range(t):
        s = rng.randrange(1, N - L)
        r = bgs_interval(N, s, s + L, m)
        g = gcd(r["value"], N) if r["value"] is not None else r["divisor"]
        out.append({"s": s, "gcd": g, "aux_hit": r["aux_hit"]})
    return out


# --------------------------------------------------------------------- M2
def bernoulli_mod(N: int, e: int, m: Metrics) -> List[int]:
    """B_0..B_e mod N (B_1 = -1/2) by the recurrence sum_{k<=n} C(n+1,k) B_k = 0.

    Needs 1/(n+1) for n <= e: gcd-checked; raises NonUnit for a nonunit n+1
    (i.e. a prime factor of N that is <= e+1 -- trial-division-type discovery).
    """
    B = [1] + [0] * e
    binom_row = [1, 1]  # C(1, .); updated to C(n+1, .) at step n
    for n in range(1, e + 1):
        binom_row = [1] + [(binom_row[k - 1] + binom_row[k]) % N for k in range(1, n + 1)] + [1]
        s = 0
        for k in range(n):
            s = (s + binom_row[k] * B[k]) % N
        m.inc("rn_mul", n)
        m.inc("rn_add", 2 * n)
        B[n] = (-s * inv_mod(n + 1, N, m)) % N
        m.inc("rn_mul")
    return B


def _binom_row_upto(e1: int, N: int) -> List[int]:
    row = [1]
    for n in range(1, e1 + 1):
        row = [1] + [(row[k - 1] + row[k]) % N for k in range(1, n)] + [1]
    return row


def faulhaber_prefix(N: int, n: int, e: int, B: List[int], m: Metrics) -> int:
    """sum_{j=0}^{n-1} j^e mod N  =  1/(e+1) sum_k C(e+1,k) B_k n^(e+1-k)."""
    row = _binom_row_upto(e + 1, N)
    s = 0
    npow = 1
    pw = [1] * (e + 2)
    for i in range(1, e + 2):
        npow = npow * n % N
        pw[i] = npow
    for k in range(e + 1):
        s = (s + row[k] * B[k] % N * pw[e + 1 - k]) % N
    m.inc("rn_mul", 3 * (e + 1))
    return s * inv_mod(e + 1, N, m) % N


def power_sum_faulhaber(N: int, A: int, B_: int, e: int, m: Optional[Metrics] = None) -> Dict:
    m = m if m is not None else Metrics()
    try:
        Bn = bernoulli_mod(N, e, m)
        S = (faulhaber_prefix(N, B_, e, Bn, m) - faulhaber_prefix(N, A, e, Bn, m)) % N
    except NonUnit as ex:
        return {"S": None, "divisor": ex.divisor, "witness": ex.value, "kind": "denominator",
                **m.to_record()}
    m.inc("gcd")
    g = gcd(S, N)
    return {"S": S, "divisor": g if g > 1 else None, "witness": S if g > 1 else None,
            "kind": "value" if g > 1 else None, **m.to_record()}


def power_sum_direct(N: int, A: int, B_: int, e: int, m: Optional[Metrics] = None) -> Dict:
    m = m if m is not None else Metrics()
    S = 0
    for j in range(A, B_):
        S = (S + pow(j, e, N)) % N
    m.inc("rn_mul", (B_ - A) * 2 * max(1, e.bit_length()))
    m.inc("gcd")
    g = gcd(S, N)
    return {"S": S, "divisor": g if g > 1 else None, **m.to_record()}


# --------------------------------------------------------------------- M3
def truncated_FN_test(N: int, K: int, m: Optional[Metrics] = None) -> Dict:
    """Z(N,K): is (1+X)^N == 1 mod (N, X^K)?  (equals U(K) for composite N).

    If not, returns gcd(first nonzero coefficient index k, N) -- by BRIEF 2.2
    a nonzero coefficient a_k forces 1 < gcd(k,N) < N.
    """
    m = m if m is not None else Metrics()
    Q = FlintQuotient(N, [0] * K + [1], m)
    c = Q.coeffs(power_of_linear(Q, N))
    nz = [k for k in range(1, K) if c[k] % N]
    if not nz:
        return {"Z": True, "divisor": None, **m.to_record()}
    k = nz[0]
    m.inc("gcd")
    return {"Z": False, "divisor": gcd(k, N), "first_k": k, **m.to_record()}


# --------------------------------------------------------------------- M4
def jacobi(a: int, n: int) -> int:
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


def jacobi_prefix_product(N: int, K: int, m: Optional[Metrics] = None) -> Dict:
    """C_N(K) = prod_{1<=j<K} (j/N) in {-1,0,1}; cost K Jacobi symbols (odd N)."""
    m = m if m is not None else Metrics()
    c = 1
    for j in range(1, K):
        m.inc("jacobi")
        c *= jacobi(j, N)
        if c == 0:
            return {"C": 0, "witness": j, "divisor": gcd(j, N), **m.to_record()}
    return {"C": c, "witness": None, "divisor": None, **m.to_record()}
