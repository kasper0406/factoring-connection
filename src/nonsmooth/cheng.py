"""OPERATIONAL: CM / anomalous-curve factoring (Cheng 2002 "4p - 1 method"; Shirase 2017).

Input: N only (plus public config and a seed). Never sees p or q.

Mechanism (KNOWN). Let -D be a class-number-one fundamental discriminant, D = 3 mod 8 squarefree
(D in {3, 11, 19, 43, 67, 163}), and E_D a curve over Q with CM by the maximal order of
Q(sqrt(-D)). (D = 7 is excluded: 4p - 1 = 3 mod 8 for odd p forces D = 3 mod 8.) If a prime p | N satisfies 4p - 1 = D s^2, then the Frobenius trace of E_D mod p is
t = +-1 (for D > 4: 4p = t^2 + D v^2 has the unique solution t = +-1 up to sign), so E_D mod p or
its quadratic twist has exactly p points: it is anomalous. For D = 3 (j = 0) one of the six twists
y^2 = x^3 + b is anomalous. Take a random x in Z/N and compute [N] x with x-only (X:Z) arithmetic,
which is valid on E and its quadratic twist simultaneously. If x lies (mod p) on the anomalous twist
(probability ~1/2), [N]P = [q]([p]P) = O mod p, so p | Z. Mod q the same happens only if the order
of the point mod q divides N (rare). No smoothness is used: the group order is p itself.

Cost: O(log N) ring operations per try. The method only applies to the thin set of p with
4p - 1 = D s^2 for small D (class polynomial of degree h(-D) in Cheng's general version).

Counters follow AGENT_RULES §4 (ring ops, gcds, max operand bits).
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from math import gcd

# j-invariants of the class-number-one orders with D = 3 mod 4 (fundamental discriminant -D).
J_INVARIANTS = {
    11: -32768,
    19: -884736,
    43: -884736000,
    67: -147197952000,
    163: -262537412640768000,
}
D_LIST = (3, 11, 19, 43, 67, 163)


@dataclass
class Counters:
    mul: int = 0
    add: int = 0
    gcd: int = 0
    tries: int = 0
    max_bits: int = 0
    log: list = field(default_factory=list)


def curve_coeffs(D, rng, N):
    """(a, b) of a short Weierstrass curve with CM by Q(sqrt(-D)), as integers mod N."""
    if D == 3:
        return 0, rng.randrange(1, N)          # y^2 = x^3 + b, random sextic twist
    j = J_INVARIANTS[D]
    k = 1728 - j
    a = 3 * j * k
    b = 2 * j * k * k
    return a % N, b % N


def _xdbl(X, Z, a, b, N, c):
    X2, Z2 = X * X % N, Z * Z % N
    XZ = X * Z % N
    t = (X2 - a * Z2) % N
    Xo = (t * t - 8 * b * XZ * Z2) % N
    Zo = 4 * Z * (X2 * X + a * X * Z2 + b * Z2 * Z) % N
    c.mul += 12
    c.add += 6
    return Xo, Zo


def _xadd(Xm, Zm, Xn, Zn, X0, Z0, a, b, N, c):
    """Differential addition: (m+n) from m, n and m-n = (X0:Z0)."""
    XX, ZZ = Xm * Xn % N, Zm * Zn % N
    t = (XX - a * ZZ) % N
    u = (Xm * Zn + Xn * Zm) % N
    Xo = Z0 * ((t * t - 4 * b * ZZ * u) % N) % N
    v = (Xm * Zn - Xn * Zm) % N
    Zo = X0 * (v * v % N) % N
    c.mul += 12
    c.add += 5
    return Xo, Zo


def ladder_simple(x, k, a, b, N, c):
    """[k](x:1): standard Montgomery ladder with invariant R1 - R0 = P."""
    if k == 0:
        return (1, 0)
    R0 = (x % N, 1)
    R1 = _xdbl(R0[0], R0[1], a, b, N, c)
    for bit in bin(k)[3:]:
        if bit == '1':
            R0 = _xadd(*R0, *R1, x, 1, a, b, N, c)
            R1 = _xdbl(*R1, a, b, N, c)
        else:
            R1 = _xadd(*R0, *R1, x, 1, a, b, N, c)
            R0 = _xdbl(*R0, a, b, N, c)
    return R0


def cheng_factor(N, seed=0, tries_per_D=24, d_list=D_LIST, tries_D3=64):
    """Try the anomalous-CM test for each D in d_list. Returns (factor or None, D, Counters).

    Every nonunit encountered is output (rule 3): gcd is taken on the curve discriminant, on the
    final Z coordinate, and on x itself.
    """
    rng = random.Random(f"cheng|{N}|{seed}")
    c = Counters()
    for g0 in (2, 3):
        if N % g0 == 0 and N != g0:
            return g0, None, c
    for D in d_list:
        for _ in range(tries_D3 if D == 3 else tries_per_D):
            c.tries += 1
            a, b = curve_coeffs(D, rng, N)
            disc = (4 * a * a * a + 27 * b * b) % N
            g = gcd(disc, N)
            c.gcd += 1
            if 1 < g < N:
                c.log.append(("disc", D))
                return g, D, c
            x = rng.randrange(2, N)
            X, Z = ladder_simple(x, N, a, b, N, c)
            c.max_bits = max(c.max_bits, N.bit_length())
            g = gcd(Z, N)
            c.gcd += 1
            if 1 < g < N:
                c.log.append(("Z", D))
                return g, D, c
    return None, None, c
