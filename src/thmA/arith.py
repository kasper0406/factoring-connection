"""Integer and small-field helpers for the toy Theorem-A pipeline (agent thmA).  OPERATIONAL.

Nothing here factors N or uses factors of N.  Primality: deterministic Miller-Rabin with the first 12 prime
bases (exact for n < 3.3e24; a probable-prime test above that; the theorem uses AKS - disclosed).  The
Pollard-Brent routine is used ONLY on the split halves |u0|, |v0| <= sqrt(2N) as a toy stand-in for the
Lenstra-Pomerance shape test, never on N.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Dict, List, Tuple

_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n + 1) if sieve[i]]


def iroot(n: int, e: int) -> int:
    """floor(n^(1/e)) exactly."""
    if n < 0:
        raise ValueError
    if n < 2:
        return n
    x = 1 << ((n.bit_length() + e - 1) // e)
    while True:
        y = ((e - 1) * x + n // x ** (e - 1)) // e
        if y >= x:
            break
        x = y
    while x ** e > n:
        x -= 1
    while (x + 1) ** e <= n:
        x += 1
    return x


def perfect_power(n: int) -> Tuple[int, int] | None:
    """(m, e) with n = m^e, e >= 2 maximal-first search over prime e; None if n is not a perfect power."""
    for e in range(2, n.bit_length() + 1):
        m = iroot(n, e)
        if m > 1 and m ** e == n:
            return m, e
    return None


def pollard_brent(n: int, rng: random.Random) -> int:
    """A nontrivial factor of composite n (Las Vegas).  Toy shape-test stand-in only."""
    if n % 2 == 0:
        return 2
    while True:
        y, c, m = rng.randrange(1, n), rng.randrange(1, n), 64
        g, r, q = 1, 1, 1
        x = ys = y
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = gcd(q, n)
                k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = gcd(abs(x - ys), n)
        if g != n:
            return g


def factor_small(n: int, rng: random.Random, small_primes: List[int]) -> Dict[int, int]:
    """Complete factorisation of a SPLIT HALF n (|u0| or |v0|, at most ~sqrt(2N)).  Shape-test stand-in."""
    out: Dict[int, int] = {}
    for p in small_primes[:60]:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            out[m] = out.get(m, 0) + 1
            continue
        pp = perfect_power(m)
        if pp:
            stack.extend([pp[0]] * pp[1])
            continue
        f = pollard_brent(m, rng)
        stack.extend([f, m // f])
    return out


def gauss_reduce(u: Tuple[int, int], v: Tuple[int, int]) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Lagrange-Gauss reduction of a 2D lattice basis (exact). Returns (b1, b2) with |b1| <= |b2|."""
    def n2(x):
        return x[0] * x[0] + x[1] * x[1]
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        nu = n2(u)
        if nu == 0:
            raise ValueError("degenerate basis")
        dot = u[0] * v[0] + u[1] * v[1]
        # nearest integer to dot/nu, exactly
        m = (2 * dot + nu) // (2 * nu)
        v = (v[0] - m * u[0], v[1] - m * u[1])
        if n2(v) >= nu:
            return u, v
        u, v = v, u


def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def sqrt_mod_p(a: int, p: int, rng: random.Random) -> int | None:
    """Tonelli-Shanks in F_p (Las Vegas), p odd prime; None if a is a non-residue."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    while True:
        z = rng.randrange(2, p)
        if pow(z, (p - 1) // 2, p) == p - 1:
            break
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r
