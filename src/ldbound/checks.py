"""Exact checks of the identities and inequalities behind notes/ldbound_results.md.

OPERATIONAL (N, public parameters, seed only):
- field_w(N, w, d, g0, g1, R): f = gbe(N w) + (g1 x - g0) R.  F(g0, g1) = N w.
DIAGNOSTIC (pure arithmetic; no factorization of N is ever used):
- lemma12_w: g1^d F(P) == N w b^d (mod G(P)).
- pair_count_direct / pair_count_ramanujan: #{(F,F'): F == F' (mod m)} computed directly and via
  (1/m) sum_{m'|m} sum_{F,F'} c_{m'}(F - F')  (Ramanujan sums; = the additive Parseval identity
  Q(m) = (1/m) sum_{m'|m} sum_{(a,m')=1} |S(a/m')|^2, in exact integer form).
- bt_ratio(t, W): max_a #{w prime in (W,2W], w == a (t)} * phi(t) / #{w}  versus the
  Montgomery-Vaughan Brun-Titchmarsh bound 4 W log W' ... (see bt_bound).
- lift_chain(...): for one toy pair (P, P', g, w): Z <= (Q + Q')/2, Q <= Q0(t, Delta, r0),
  and sum over all reduced r0 of Q0 <= Q_full(t Delta).
"""
from __future__ import annotations

from math import gcd, log
from typing import Dict, List, Sequence

import sympy

from src.pcest.geometry import gbe, hom_eval


def field_w(N: int, w: int, d: int, g0: int, g1: int, R: Sequence[int]) -> List[int]:
    f = gbe(N * w, d, g0, g1)
    assert len(R) == d
    for j, r in enumerate(R):
        f[j + 1] += g1 * r
        f[j] -= g0 * r
    return f


def lemma12_w(f: Sequence[int], N: int, w: int, g0: int, g1: int, a: int, b: int) -> bool:
    d = len(f) - 1
    G = g1 * a - g0 * b
    if G == 0:
        return True
    return (pow(g1, d) * hom_eval(f, a, b) - N * w * b ** d) % abs(G) == 0


def ramanujan(m: int, n: int) -> int:
    """c_m(n) = sum_{e | gcd(m, n)} e mu(m/e)."""
    g = gcd(m, n)
    return sum(e * sympy.mobius(m // e) for e in sympy.divisors(g))


def pair_count_direct(vals: Sequence[int], m: int) -> int:
    cnt: Dict[int, int] = {}
    for v in vals:
        cnt[v % m] = cnt.get(v % m, 0) + 1
    return sum(c * c for c in cnt.values())


def pair_count_ramanujan(vals: Sequence[int], m: int) -> int:
    tot = 0
    for mp in sympy.divisors(m):
        for v in vals:
            for v2 in vals:
                tot += ramanujan(mp, v - v2)
    assert tot % m == 0
    return tot // m


def bt_ratio(t: int, W: int) -> Dict[str, float]:
    ws = list(sympy.primerange(W + 1, 2 * W + 1))
    cnt: Dict[int, int] = {}
    for w in ws:
        cnt[w % t] = cnt.get(w % t, 0) + 1
    phi = int(sympy.totient(t))
    mx = max(cnt.values())
    bound = 4 * W / (phi * log(2 * W / t))       # Montgomery-Vaughan: pi(2W; t, a) <= this
    return dict(max_class=mx, bt_bound=bound, ratio_to_uniform=mx * phi / len(ws),
                n_primes=len(ws), ok=mx <= bound)


def smooth(n: int, B: int) -> bool:
    n = abs(n)
    if n == 0:
        return False
    for p in sympy.primerange(2, B + 1):
        while n % p == 0:
            n //= p
        if n == 1:
            return True
    return n == 1


def smooth_sieve(X: int, B: int):
    """Set of B-smooth integers in [1, X] (exact; numpy sieve)."""
    import numpy as np
    v = np.arange(X + 1, dtype=np.int64)
    for p in sympy.primerange(2, B + 1):
        pk = p
        while pk <= X:
            v[::pk] //= p
            pk *= p
    return set(int(i) for i in np.nonzero(v == 1)[0] if i >= 1)


def lift_chain(N: int, w: int, d: int, g: Sequence[int], P, P2, T: int, B: int,
               Xstar: int = None, smooth_set=None) -> Dict[str, int]:
    """Exact check of steps 2-4 of the LD* proof on one pair.

    A_j(P) = #{z in [-T,T], z == j (|Delta|), Fhat(P) + G z B-smooth}; Z = sum_j A_j(P) A_{kappa^(d-1) j}(P');
    Q(P) = sum_j A_j(P)^2.  Q0 = #{(F,F') smooth in [-X*,X*]: F == F' == r0 (t), F == F' (t |Delta|)}.
    """
    g0, g1 = g
    a, b = P
    a2, b2 = P2
    Delta = a * b2 - a2 * b
    D = abs(Delta)
    fhat = gbe(N * w, d, g0, g1)
    G, G2 = g1 * a - g0 * b, g1 * a2 - g0 * b2
    F0, F02 = hom_eval(fhat, a, b), hom_eval(fhat, a2, b2)
    if Xstar is None:
        Xstar = max(abs(F0), abs(F02)) + max(abs(G), abs(G2)) * T
    if smooth_set is None:
        smooth_set = smooth_sieve(Xstar, B)
    isS = lambda x: abs(x) in smooth_set
    A = [0] * D
    A2 = [0] * D
    for z in range(-T, T + 1):
        if F0 + G * z != 0 and isS(F0 + G * z):
            A[z % D] += 1
        if F02 + G2 * z != 0 and isS(F02 + G2 * z):
            A2[z % D] += 1
    # kappa: P' == kappa P (mod Delta)
    kappa = None
    for k0 in range(D):
        if (a2 - k0 * a) % D == 0 and (b2 - k0 * b) % D == 0:
            kappa = k0
            break
    assert kappa is not None
    Z = sum(A[j] * A2[(pow(kappa, d - 1, D) * j) % D] for j in range(D))
    Q, Q2 = sum(x * x for x in A), sum(x * x for x in A2)
    t = abs(G)
    r0 = F0 % t
    vals = [F for F in range(-Xstar, Xstar + 1) if F != 0 and isS(F)]
    cls: Dict[int, int] = {}
    for F in vals:
        if F % t == r0:
            cls[F % (t * D)] = cls.get(F % (t * D), 0) + 1
    Q0 = sum(c * c for c in cls.values())
    # all reduced r0' mod t
    byres: Dict[int, int] = {}
    for F in vals:
        if gcd(F, t) == 1:
            byres[F % (t * D)] = byres.get(F % (t * D), 0) + 1
    Q0_all = sum(c * c for c in byres.values())
    Qfull = pair_count_direct(vals, t * D)
    return dict(Delta=Delta, Xstar=Xstar, Z=Z, Q=Q, Q2=Q2, Q0=Q0, Q0_all=Q0_all, Qfull=Qfull,
                r0_ok=int((pow(g1, d, t) * F0 - N * w * pow(b, d, t)) % t == 0),
                chain_ok=int(2 * Z <= Q + Q2 and Q <= Q0 and Q0 <= Q0_all <= Qfull))
