"""OPERATIONAL: elliptic-curve arithmetic over Z/N and reductions from point counts to factoring.

Inputs to every function here are N, public curve data, and the counts (or, for the Dieulefait-Urroz malleability
theorem, the factorisation of the count, which is the oracle answer of their Game 1).  Nothing here reads p, q,
calls an integer-factorisation routine, or imports harness code.  Primality is never tested.

Contents
- Split: exception carrying a nontrivial factor of N.
- Group law on short Weierstrass curves y^2 = x^3 + a x + b over Z/N (N odd squarefree, gcd(N, 6 disc) = 1):
  points are None (= O in both CRT components) or (x, y) affine in both components; any situation in which the
  two components differ in shape surfaces as a non-invertible denominator with 1 < gcd < N and raises Split.
- tw_trial / tw_factor: Theorem TW (one-curve count => factoring) via a point on the twist E^(delta).
- sr2_supersingular: the (p+1)(q+1) fallback (SR2 instance).
- coppersmith_linear: Howgrave-Graham small root of x + u modulo an unknown divisor b >= N^beta of N (flint LLL).
- du_affine_twist, du_projective_twist_repaired, du_malleability: Dieulefait-Urroz reductions (arXiv:1911.11004).
"""
from __future__ import annotations

import math
import random
from math import gcd, isqrt
from typing import Dict, Iterable, List, Optional, Tuple

import gmpy2
import flint

from src.selfred.rigidity import integer_roots

Point = Optional[Tuple[int, int]]


class Split(Exception):
    def __init__(self, factor: int, where: str = ""):
        super().__init__(f"split {factor} at {where}")
        self.factor = int(factor)
        self.where = where


class Counters:
    def __init__(self):
        self.adds = 0
        self.dbls = 0
        self.invs = 0
        self.gcds = 0
        self.lll_calls = 0
        self.max_lll_dim = 0

    def as_dict(self):
        return dict(self.__dict__)


def _inv(u: int, N: int, ctr: Optional[Counters], where: str) -> int:
    g = gcd(u % N, N)
    if ctr is not None:
        ctr.gcds += 1
        ctr.invs += 1
    if g != 1:
        if g != N:
            raise Split(g, where)
        raise ZeroDivisionError(where)
    return int(gmpy2.invert(u, N))


def ec_add(P: Point, Q: Point, a: int, N: int, ctr: Optional[Counters] = None) -> Point:
    """P + Q on y^2 = x^3 + a x + b over Z/N (b is not needed). Raises Split on a nontrivial gcd.

    Correctness (componentwise, N = prod of distinct primes, all > 3):
    - gcd(x2 - x1, N) = 1: x1 != x2 in every component, the chord formula is valid and the sum is affine everywhere.
    - gcd = N: x1 = x2 everywhere, so y1 = +-y2 in each component.  h = gcd(y1 + y2, N):
      h = N -> y1 = -y2 everywhere, sum = O everywhere; h = 1 -> y1 = y2 and 2*y1 != 0 everywhere, tangent formula.
    - any other gcd is a proper factor.
    """
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if ctr is not None:
        ctr.adds += 1
        ctr.gcds += 1
    g = gcd((x2 - x1) % N, N)
    if g == 1:
        lam = (y2 - y1) * int(gmpy2.invert(x2 - x1, N)) % N
    elif g != N:
        raise Split(g, "add:x")
    else:
        if ctr is not None:
            ctr.gcds += 1
        h = gcd((y1 + y2) % N, N)
        if h == N:
            return None
        if h != 1:
            raise Split(h, "add:y")
        lam = (3 * x1 * x1 + a) * int(gmpy2.invert(2 * y1, N)) % N
    x3 = (lam * lam - x1 - x2) % N
    y3 = (lam * (x1 - x3) - y1) % N
    return (x3, y3)


def ec_mul(k: int, P: Point, a: int, N: int, ctr: Optional[Counters] = None) -> Point:
    """[k]P by left-to-right double-and-add (k >= 0). Raises Split."""
    if k < 0:
        raise ValueError("k >= 0")
    R: Point = None
    for bit in bin(k)[2:]:
        R = ec_add(R, R, a, N, ctr)
        if ctr is not None:
            ctr.dbls += 1
        if bit == "1":
            R = ec_add(R, P, a, N, ctr)
    return R


def f_eval(x: int, a: int, b: int, N: int) -> int:
    return (x * x * x + a * x + b) % N


def tw_trial(N: int, a: int, b: int, M: int, rng: random.Random,
             ctr: Optional[Counters] = None, jacobi_filter: bool = True) -> Tuple[Optional[int], str]:
    """One trial of algorithm TW.  Returns (factor or None, outcome label).

    outcome in {"delta_gcd" (f(x0) shares a factor with N), "jacobi_ne_-1", "split" (from [M]P), "killed"
    ([M]P = O in both components, no split), "affine" (cannot happen if M is a multiple of #E(F_p): reported)}.
    """
    x0 = rng.randrange(N)
    delta = f_eval(x0, a, b, N)
    g = gcd(delta, N)
    if ctr is not None:
        ctr.gcds += 1
    if 1 < g < N:
        return g, "delta_gcd"
    if g == N:
        return None, "delta_zero"
    if jacobi_filter and gmpy2.jacobi(delta, N) != -1:
        return None, "jacobi_ne_-1"
    # (TW*: without the filter the point lies on E or E' independently at each prime; the chain's trailing
    # doublings then also perform Miller's 2-adic descent, cf. codex/ecount_1.out.)
    # P = (delta*x0, delta^2) on E^(delta): Y^2 = X^3 + a delta^2 X + b delta^3
    ad = a * delta * delta % N
    P = (delta * x0 % N, delta * delta % N)
    try:
        R = ec_mul(M, P, ad, N, ctr)
    except Split as s:
        return s.factor, "split"
    if R is None:
        return None, "killed"
    return None, "affine"


def tw_factor(N: int, a: int, b: int, M: int, rng: random.Random, max_trials: int = 64,
              ctr: Optional[Counters] = None, jacobi_filter: bool = True) -> Dict:
    """Run TW until a split or max_trials. Verifies the factor divides N (no primality test)."""
    outcomes: Dict[str, int] = {}
    for t in range(1, max_trials + 1):
        fac, out = tw_trial(N, a, b, M, rng, ctr, jacobi_filter)
        outcomes[out] = outcomes.get(out, 0) + 1
        if fac is not None:
            assert 1 < fac < N and N % fac == 0
            return {"factor": fac, "trials": t, "outcomes": outcomes}
    return {"factor": None, "trials": max_trials, "outcomes": outcomes}


def sr2_supersingular(N: int, M: int) -> Optional[int]:
    """If M = (p+1)(q+1), recover p from s = M - N - 1 = p + q (SR2 instance)."""
    s = M - N - 1
    disc = s * s - 4 * N
    if disc < 0:
        return None
    r = isqrt(disc)
    if r * r != disc or (s + r) % 2:
        return None
    for c in ((s - r) // 2, (s + r) // 2):
        if 1 < c < N and N % c == 0:
            return c
    return None


# ---------------------------------------------------------------------------------------------------------------
# Coppersmith / Howgrave-Graham for a linear polynomial x + u modulo an unknown divisor b >= 2^blog of N.
# ---------------------------------------------------------------------------------------------------------------

def _hg_xbound_bits(n_bits: int, blog: float, m: int, t: int) -> float:
    """Largest log2 X for which the Howgrave-Graham lattice (dim d = m + t) provably yields the root:
    det^(1/d) * 2^((d-1)/4) * sqrt(d) < b^m, det = X^(d(d-1)/2) * N^(m(m+1)/2)."""
    d = m + t
    rhs = m * blog - (d - 1) / 4.0 - 0.5 * math.log2(d) - 1.0
    lhs_N = m * (m + 1) / 2.0 * n_bits / d
    return 2.0 * (rhs - lhs_N) / (d - 1)


def coppersmith_linear(N: int, u: int, X: int, m: int, t: int, ctr: Optional[Counters] = None) -> List[int]:
    """Return candidate small x (|x| <= X) with gcd(x + u, N) possibly a large divisor. Caller checks gcds."""
    u %= N
    # basis polynomials in y = x (coefficient lists, low degree first), evaluated at y*X
    polys: List[List[int]] = []
    fpow = [[1]]
    for i in range(1, m + 1):
        prev = fpow[-1]
        nxt = [0] * (len(prev) + 1)
        for k, c in enumerate(prev):
            nxt[k] += c * u
            nxt[k + 1] += c
        fpow.append(nxt)
    for i in range(m + 1):
        polys.append([c * N ** (m - i) for c in fpow[i]])
    for j in range(1, t):
        polys.append([0] * j + fpow[m])
    d = len(polys)
    rows = []
    for P in polys:
        row = [0] * d
        for k, c in enumerate(P):
            row[k] = c * X ** k
        rows.append(row)
    B = flint.fmpz_mat(rows).lll()
    if ctr is not None:
        ctr.lll_calls += 1
        ctr.max_lll_dim = max(ctr.max_lll_dim, d)
    cands: List[int] = []
    for r in range(min(2, d)):
        coeffs = [int(B[r, k]) for k in range(d)]
        poly = []
        ok = True
        for k, c in enumerate(coeffs):
            if c % (X ** k):
                ok = False
                break
            poly.append(c // X ** k)
        if not ok or all(c == 0 for c in poly):
            continue
        while poly and poly[-1] == 0:
            poly.pop()
        if len(poly) < 2:
            continue
        for x in integer_roots(poly, -X, X):
            cands.append(x)
    return cands


def find_divisor_near(N: int, center: int, W: int, blog: float, ctr: Optional[Counters] = None,
                      m: int = 5, force_lattice: bool = False) -> Optional[int]:
    """Find a divisor b >= 2^blog of N with |b - center| <= W, by brute force for small W and Coppersmith chunks
    otherwise.  Deterministic.  Returns b or None."""
    n = N.bit_length()
    if W <= 4096 and not force_lattice:
        for x in range(-W, W + 1):
            g = gcd(center + x, N)
            if 1 < g < N and g.bit_length() >= int(blog):
                return g
        return None
    t = m
    xb = _hg_xbound_bits(n, blog, m, t)
    while xb < 2 and m < 30:
        m += 1
        t = m
        xb = _hg_xbound_bits(n, blog, m, t)
    X = max(2, int(2 ** xb))
    c = center - W + X
    while c - X <= center + W:
        for x in coppersmith_linear(N, c, X, m, t, ctr):
            g = gcd(c + x, N)
            if 1 < g < N:
                return g
        c += 2 * X + 1
    return None


# ---------------------------------------------------------------------------------------------------------------
# Dieulefait-Urroz reductions
# ---------------------------------------------------------------------------------------------------------------

def du_affine_twist(N: int, Eaff: int, Edaff: int, d: int) -> Optional[int]:
    """D-U Thm 4, affine counts E* = (p - a_p)(q - a_q) and E_d* for the twist by d, (d/p) = -1 promised.
    Case 2 ((d/N) = -1): E* + E_d* = 2q(p - a_p) or 2p(q - a_q); take gcd with N.
    Case 1 ((d/N) = +1, so (d/q) = -1 too): E_d* = (p + a_p)(q + a_q); the other two twist counts are the roots of
    X^2 - (4N - L) X + E* E_d* (Lemma 3, affine: the sum of the four counts is 4N); then gcd(E* + root, N)."""
    cands: List[int] = []
    if gmpy2.jacobi(d % N, N) == -1:
        cands.append(Eaff + Edaff)
    else:
        L = Eaff + Edaff
        S = 4 * N - L
        disc = S * S - 4 * Eaff * Edaff
        if disc >= 0:
            r = isqrt(disc)
            if r * r == disc and (S + r) % 2 == 0:
                for root in ((S + r) // 2, (S - r) // 2):
                    cands.append(Eaff + root)
    for c in cands:
        g = gcd(c, N)
        if 1 < g < N:
            return g
    return None


def du_projective_lemma_as_written(N: int, E: int, Ed: int, PQ: int) -> Tuple[int, int]:
    """D-U Lemma 3, case 1, projective, AS WRITTEN: it needs PQ = (p+1)(q+1) as an input.  Exposed only to show the
    dependence: without PQ the step cannot be executed (and with PQ, N is already factored)."""
    L = E + Ed
    S = 4 * PQ - L
    r = isqrt(S * S - 4 * E * Ed)
    return (S + r) // 2, (S - r) // 2


def du_projective_twist_repaired(N: int, E: int, Ed: int, d: int, ctr: Optional[Counters] = None,
                                 small_k: int = 4096) -> Optional[int]:
    """Projective case 2 ((d/N) = -1), repaired (no PQ needed).
    If (d/p) = -1: E = (P - a)A_q, Ed = (P + a)A_q with P = p + 1, so (E + Ed)/(Ed - E) = P/a; with
    r = (E + Ed)/gcd(E + Ed, Ed - E) = P/k, k = gcd(P, a) <= |a| <= 2 sqrt(p): p = k r - 1.  Search k: small k
    by gcd, then Coppersmith on y = k with f(y) = r y - 1 (monic: y - r^{-1} mod N).  Symmetric if (d/q) = -1."""
    if gmpy2.jacobi(d % N, N) != -1:
        return None
    S = E + Ed
    D = Ed - E
    if D == 0:
        return None
    r = S // gcd(S, abs(D))
    for k in range(1, small_k + 1):
        g = gcd(k * r - 1, N)
        if 1 < g < N:
            return g
    g = gcd(r, N)
    if 1 < g < N:
        return g
    rinv = int(gmpy2.invert(r, N))
    u = (-rinv) % N  # r*y - 1 = r*(y - rinv): divisor of N divides y + u
    n = N.bit_length()
    kmax = 2 * isqrt(isqrt(N) * 2) + 4  # |a| <= 2 sqrt(p) <= 2 (2N)^{1/4}
    half = (kmax + 1) // 2
    res = find_divisor_near(N, (u + half) % N, half + 1, (n - 2) / 2.0, ctr, force_lattice=small_k == 0)
    if res is None:
        return None
    return res


def divisors_in_range(fac: Dict[int, int], lo: int, hi: int, cap: int = 10 ** 6) -> Tuple[List[int], int]:
    """All divisors of prod p^e in [lo, hi] (enumerates all divisors; returns (list, total divisor count))."""
    divs = [1]
    for pr, e in fac.items():
        divs = [dv * pr ** k for dv in divs for k in range(e + 1)]
        if len(divs) > cap:
            raise OverflowError("too many divisors")
    return sorted(x for x in divs if lo <= x <= hi), len(divs)


def du_malleability(N: int, M_factorisation: Dict[int, int], ctr: Optional[Counters] = None) -> Dict:
    """D-U Thm 6: from the factorisation of M = #E(Z/N) (oracle answer), find a divisor u of M within
    2 sqrt(u) + 2 of a prime factor of N, then Coppersmith (approximation of the factor to within ~N^{1/4})."""
    n = N.bit_length()
    s = isqrt(N)
    cands, dM = divisors_in_range(M_factorisation, s // 4, 4 * s + 4)
    tried = 0
    for u in cands:
        tried += 1
        W = 2 * isqrt(u) + 3
        g = find_divisor_near(N, u, W, (n - 2) / 2.0, ctr)
        if g is not None:
            return {"factor": g, "divisors_total": dM, "divisors_in_window": len(cands), "tried": tried}
    return {"factor": None, "divisors_total": dM, "divisors_in_window": len(cands), "tried": tried}
