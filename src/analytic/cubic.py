"""Complete cubic Weyl sums mod N factor N (reduction R3; see notes/analytic_results.md Sec. 5).

    S(a, N) = sum_{x mod N} e(a x^3 / N),   e(t) = exp(2 pi i t).

For N = p q (p != q primes, p, q > 3, gcd(a, N) = 1):
  * CRT: S(a, N) = S(a qbar, p) * S(a pbar, q)  (qbar = q^{-1} mod p, pbar = p^{-1} mod q);
  * p = 2 mod 3: x -> x^3 permutes F_p, so S(b, p) = 0;
  * p = 1 mod 3: S(b, p) = 2 Re(chibar(b) G(chi)) = 2 sqrt(p) cos(theta_p - 2 pi k(b)/3), where chi
    is a cubic character mod p, chi(b) = e(k(b)/3), G(chi) = sqrt(p) e^{i theta_p}, and
    G(chi)^3 = p J(chi, chi) with J(chi, chi) = pi_p in Z[w] of norm p.  So arg(pi_p) = 3 theta_p.
Hence S(a^i b^j, N) = 4 sqrt(N) u_{alpha(i,j)} v_{beta(i,j)} with u_k = cos(theta_p - 2 pi (k+s)/3),
a rank-one 3x3 array once the cubic-class labels (alpha, beta) are guessed (48 candidates).
z_u = sum_k u_k w^k is proportional to e^{i theta_p} (up to a real scale, a cyclic shift = unit,
and complex conjugation), so z_u^3 is proportional to pi_p or its conjugate.  Then
gamma = z_u^3 z_v^3 and gamma' = z_u^3 conj(z_v^3) are real multiples of pi_p pi_q and
pi_p conj(pi_q) (up to units/conjugation), both of norm N: rescale to norm N, round in Z[w], and
gcd(gamma, gamma') in the Euclidean ring Z[w] has norm p.

Code split:
  * `weyl_oracle(a, N)` -- HARNESS oracle: brute force O(N), float64 (a disclosed exception to
    the exact-arithmetic rule: the oracle's output is an inherently real analytic quantity; its
    rounding error is <= N * 2^-45 and is far below the noise levels studied).  It sees only
    (a, N) and does not use the factorisation.
  * `cubic_decode(N, oracle, rng)` -- OPERATIONAL decoder: sees N and oracle values (as exact
    Fractions); exact rational / Z[w] arithmetic only.
Elements of Q(w), w = e^{2 pi i/3}, are pairs (A, B) meaning A + B w; w^2 = -1 - w.
"""
from __future__ import annotations

import random
from fractions import Fraction
from math import gcd, isqrt
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

Elt = Tuple[Fraction, Fraction]


# ---------------------------------------------------------------- Q(w) / Z[w] arithmetic
def qw_mul(x: Elt, y: Elt) -> Elt:
    A, B = x
    C, D = y
    return (A * C - B * D, A * D + B * C - B * D)


def qw_conj(x: Elt) -> Elt:
    A, B = x
    return (A - B, -B)


def qw_norm(x: Elt):
    A, B = x
    return A * A - A * B + B * B


def zw_divmod(x: Tuple[int, int], y: Tuple[int, int]) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Euclidean division in Z[w]: x = qy + r with N(r) < N(y)."""
    n = y[0] * y[0] - y[0] * y[1] + y[1] * y[1]
    num = qw_mul((Fraction(x[0]), Fraction(x[1])), qw_conj((Fraction(y[0]), Fraction(y[1]))))
    best = None
    fa = (num[0] / n)
    fb = (num[1] / n)
    a0 = fa.numerator // fa.denominator
    b0 = fb.numerator // fb.denominator
    for da in (0, 1):
        for db in (0, 1):
            q = (a0 + da, b0 + db)
            qy = qw_mul((Fraction(q[0]), Fraction(q[1])), (Fraction(y[0]), Fraction(y[1])))
            r = (x[0] - int(qy[0]), x[1] - int(qy[1]))
            nr = r[0] * r[0] - r[0] * r[1] + r[1] * r[1]
            if best is None or nr < best[0]:
                best = (nr, q, r)
    assert best[0] < n
    return best[1], best[2]


def zw_gcd(x: Tuple[int, int], y: Tuple[int, int]) -> Tuple[int, int]:
    while y != (0, 0):
        _, r = zw_divmod(x, y)
        x, y = y, r
    return x


def zw_norm(x: Tuple[int, int]) -> int:
    return x[0] * x[0] - x[0] * x[1] + x[1] * x[1]


def _sqrt_round(r: Fraction) -> int:
    """round(sqrt(r)) for rational r >= 0, exact."""
    s = isqrt(r.numerator // r.denominator)
    while (2 * s + 1) ** 2 <= 4 * r:   # (s + 1/2)^2 <= r  -> round up
        s += 1
    while s > 0 and (2 * s - 1) ** 2 > 4 * r:
        s -= 1
    return s


def round_to_norm(g: Elt, N: int) -> List[Tuple[int, int]]:
    """Integer points of norm N near the ray through g (rescale g to norm N, round, +-1 search)."""
    ng = qw_norm(g)
    if ng == 0:
        return []
    lam2 = Fraction(N) / ng  # scale^2
    out = []
    coords = []
    for c in g:
        s = _sqrt_round(c * c * lam2)
        coords.append(s if c >= 0 else -s)
    A0, B0 = coords
    for dA in (-1, 0, 1):
        for dB in (-1, 0, 1):
            e = (A0 + dA, B0 + dB)
            if zw_norm(e) == N:
                out.append(e)
    return out


# ---------------------------------------------------------------- oracle (harness)
def weyl_oracle(a: int, N: int) -> float:
    """S(a, N) = sum_{x mod N} cos(2 pi a x^3 / N) (the sum is real since x -> -x). N < 2^31."""
    x = np.arange(N, dtype=np.int64)
    r = (x * x) % N
    r = (r * x) % N
    r = (r * (a % N)) % N
    return float(np.cos(2.0 * np.pi * r.astype(np.float64) / N).sum())


# ---------------------------------------------------------------- decoder (operational)
def _functionals():
    """Nonzero linear functionals on F_3^2, as pairs (c1, c2): (i, j) -> c1 i + c2 j mod 3."""
    return [(c1, c2) for c1 in range(3) for c2 in range(3) if (c1, c2) != (0, 0)]


def cubic_decode(N: int, oracle: Callable[[int], Fraction], rng: random.Random,
                 attempts: int = 6, stats: Optional[Dict] = None) -> Optional[int]:
    """Factor N = pq (p = q = 1 mod 3) from values of S(., N). Sees only N and the oracle."""
    if stats is None:
        stats = {}
    stats.setdefault("oracle_calls", 0)
    stats.setdefault("candidates", 0)
    for _ in range(attempts):
        a = rng.randrange(2, N - 1)
        b = rng.randrange(2, N - 1)
        for t in (a, b):
            g = gcd(t, N)
            if 1 < g < N:
                stats["found_in"] = "gcd_accidental"
                return g
        vals: Dict[Tuple[int, int], Fraction] = {}
        for i in range(3):
            for j in range(3):
                vals[(i, j)] = oracle(pow(a, i, N) * pow(b, j, N) % N)
                stats["oracle_calls"] += 1
        if all(v == 0 for v in vals.values()):
            stats["all_zero"] = True
            return None  # some prime factor is 2 mod 3: no information from this degree
        F = _functionals()
        for (al1, al2) in F:
            for (be1, be2) in F:
                if (al1 * be2 - al2 * be1) % 3 == 0:
                    continue
                R = [[None] * 3 for _ in range(3)]
                for i in range(3):
                    for j in range(3):
                        R[(al1 * i + al2 * j) % 3][(be1 * i + be2 * j) % 3] = vals[(i, j)]
                # u from the column of largest l1 norm, v from the row of largest l1 norm
                lstar = max(range(3), key=lambda l: sum(abs(R[k][l]) for k in range(3)))
                kstar = max(range(3), key=lambda k: sum(abs(R[k][l]) for l in range(3)))
                u = [R[k][lstar] for k in range(3)]
                v = [R[kstar][l] for l in range(3)]
                zu = (u[0] - u[2], u[1] - u[2])
                zv = (v[0] - v[2], v[1] - v[2])
                zu3 = qw_mul(qw_mul(zu, zu), zu)
                zv3 = qw_mul(qw_mul(zv, zv), zv)
                g1 = qw_mul(zu3, zv3)
                g2 = qw_mul(zu3, qw_conj(zv3))
                c1s = round_to_norm(g1, N)
                c2s = round_to_norm(g2, N)
                for e1 in c1s:
                    for e2 in c2s:
                        stats["candidates"] += 1
                        d = zw_gcd(e1, e2)
                        p = zw_norm(d)
                        if 1 < p < N and N % p == 0:
                            stats["found_in"] = "structural_zw_gcd"
                            return p
    return None
