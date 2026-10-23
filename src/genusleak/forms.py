"""genusleak forms (OPERATIONAL): binary quadratic forms of negative discriminant, exact integer arithmetic.

Sees only N (through D = -4kN) and public parameters. No factorisation routine, no class number, no group order.
Primality is tested only on small public primes ell (trial division), to build prime forms.

A form (a, b, c) means a x^2 + b x y + c y^2 with b^2 - 4ac = D < 0, a > 0.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Iterator, List, Optional, Tuple

Form = Tuple[int, int, int]

__all__ = [
    "Form", "disc", "reduce_form", "is_reduced", "compose", "inverse", "power", "identity", "sqrt_mod_prime",
    "small_primes", "prime_form", "split_primes", "is_ambiguous", "ambiguous_split", "public_ambiguous",
]


def disc(f: Form) -> int:
    a, b, c = f
    return b * b - 4 * a * c


def is_reduced(f: Form) -> bool:
    a, b, c = f
    return abs(b) <= a <= c and not (b < 0 and (abs(b) == a or a == c))


def reduce_form(f: Form) -> Form:
    """Gauss reduction of a positive definite form (Cohen, Alg. 5.4.2)."""
    a, b, c = f
    while True:
        # normalise: -a < b <= a
        if not (-a < b <= a):
            r = (a - b) // (2 * a)          # b' = b + 2 a r lies in (-a, a]
            b_new = b + 2 * a * r
            c = c + r * (b + a * r)  # c' = (b'^2 - D)/(4a) = c + r b + a r^2
            b = b_new
        if a > c:
            a, b, c = c, -b, a
            continue
        if a == c and b < 0:
            b = -b
        return (a, b, c)


def _xgcd(x: int, y: int) -> Tuple[int, int, int]:
    """Return (u, v, g) with u x + v y = g = gcd(x, y) >= 0."""
    u0, v0, u1, v1 = 1, 0, 0, 1
    while y:
        q, r = divmod(x, y)
        x, y = y, r
        u0, u1 = u1, u0 - q * u1
        v0, v1 = v1, v0 - q * v1
    if x < 0:
        x, u0, v0 = -x, -u0, -v0
    return u0, v0, x


def compose(f1: Form, f2: Form) -> Form:
    """Gaussian composition followed by reduction (Cohen, Alg. 5.4.7)."""
    D = disc(f1)
    a1, b1, c1 = f1
    a2, b2, c2 = f2
    if a1 > a2:
        a1, b1, c1, a2, b2, c2 = a2, b2, c2, a1, b1, c1
    s = (b1 + b2) // 2
    n = b2 - s
    if a2 % a1 == 0:
        y1, d = 0, a1
    else:
        u, v, d = _xgcd(a2, a1)
        y1 = u
    if s % d == 0:
        y2, x2, d1 = -1, 0, d
    else:
        x2, y2, d1 = _xgcd(s, d)
        y2 = -y2
    v1 = a1 // d1
    v2 = a2 // d1
    r = (y1 * y2 * n - x2 * c2) % v1
    b3 = b2 + 2 * v2 * r
    a3 = v1 * v2
    num = b3 * b3 - D
    assert num % (4 * a3) == 0
    return reduce_form((a3, b3, num // (4 * a3)))


def inverse(f: Form) -> Form:
    a, b, c = f
    return reduce_form((a, -b, c))


def identity(D: int) -> Form:
    if D % 4 == 0:
        return (1, 0, -D // 4)
    return (1, 1, (1 - D) // 4)


def power(f: Form, e: int) -> Form:
    D = disc(f)
    if e < 0:
        f, e = inverse(f), -e
    r = identity(D)
    base = reduce_form(f)
    while e:
        if e & 1:
            r = compose(r, base)
        e >>= 1
        if e:
            base = compose(base, base)
    return r


def small_primes(bound: int) -> List[int]:
    """Primes <= bound by the sieve of Eratosthenes (public small primes)."""
    if bound < 2:
        return []
    s = bytearray([1]) * (bound + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(bound) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(bound + 1) if s[i]]


def sqrt_mod_prime(a: int, ell: int) -> int:
    """A square root of a modulo the small odd prime ell (Tonelli-Shanks); a must be a square mod ell."""
    a %= ell
    if a == 0:
        return 0
    if pow(a, (ell - 1) // 2, ell) != 1:
        raise ValueError("not a square")
    if ell % 4 == 3:
        return pow(a, (ell + 1) // 4, ell)
    q, s = ell - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (ell - 1) // 2, ell) != ell - 1:
        z += 1
    m, c, t, r = s, pow(z, q, ell), pow(a, q, ell), pow(a, (q + 1) // 2, ell)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % ell
            i += 1
        b = pow(c, 1 << (m - i - 1), ell)
        m, c, t, r = i, b * b % ell, t * b * b % ell, r * b % ell
    return r


def prime_form(D: int, ell: int) -> Form:
    """The reduced form of the prime ideal above the odd prime ell with (D/ell) = 1, b the least root >= 0
    with b == D (mod 2)."""
    b = sqrt_mod_prime(D, ell)
    if (b - D) % 2:
        b = ell - b if b else b
    if (b - D) % 2:
        b += ell
    # b^2 == D mod 4 ell
    assert (b * b - D) % (4 * ell) == 0
    return reduce_form((ell, b, (b * b - D) // (4 * ell)))


def split_primes(D: int, count: int, start_bound: int = 1000) -> List[int]:
    """The first `count` odd primes ell not dividing D with (D/ell) = 1."""
    out: List[int] = []
    bound = start_bound
    while True:
        out = []
        for ell in small_primes(bound):
            if ell == 2 or D % ell == 0:
                continue
            if pow(D % ell, (ell - 1) // 2, ell) == 1:
                out.append(ell)
                if len(out) == count:
                    return out
        bound *= 4


def is_ambiguous(f: Form) -> bool:
    """A reduced form is ambiguous iff b = 0, a = b, or a = c (the class has order <= 2)."""
    a, b, c = f
    return b == 0 or a == b or a == c


def ambiguous_split(N: int, f: Form) -> Optional[int]:
    """For a reduced ambiguous form, return a nontrivial divisor of N read off from (a, b, c) by gcd, if any.
    Works for any form (gcd checks only)."""
    a, b, c = f
    cands = [a, b, c]
    if a == c:
        cands += [a - b // 2, a + b // 2]   # b^2 - 4a^2 = D  =>  (2a - b)(2a + b) = -D
    for x in cands:
        g = gcd(x, N)
        if 1 < g < N:
            return g
    return None


def public_ambiguous(N: int, k: int = 1) -> List[Form]:
    """Ambiguous reduced forms of D = -4kN whose first coefficient divides 4k (computable from k alone)."""
    D = -4 * k * N
    out = []
    for a in range(1, 4 * k + 1):
        if (4 * k) % a:
            continue
        for b in range(0, a + 1):
            if b not in (0, a):
                continue
            num = b * b - D
            if num % (4 * a):
                continue
            f = (a, b, num // (4 * a))
            if is_reduced(f) and gcd(gcd(f[0], f[1]), f[2]) == 1:  # primitive forms only (codex_2 finding 11)
                out.append(f)
    return out
