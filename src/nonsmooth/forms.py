"""Positive definite binary quadratic forms (a, b, c), discriminant b^2 - 4ac < 0.

Reduction and Dirichlet composition are polynomial time (OPERATIONAL-safe: they use only the
form coefficients and gcds). Class-number enumeration is O(|D|) and only for small diagnostic
inputs.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt


def disc(f):
    a, b, c = f
    return b * b - 4 * a * c


def reduce_form(f):
    """Gauss reduction of a positive definite form: |b| <= a <= c, b >= 0 if |b| = a or a = c."""
    a, b, c = f
    while True:
        if c < a:
            a, b, c = c, -b, a
            continue
        if b > a or b <= -a:
            # bring b into (-a, a]
            k = (a - b) // (2 * a)
            # b' = b + 2ka in (-a, a]
            b2 = b + 2 * k * a
            c = (b2 * b2 - (b * b - 4 * a * c)) // (4 * a)
            b = b2
            continue
        if a == c and b < 0:
            b = -b
        return (a, b, c)


def _xgcd(a, b):
    """Return (g, x, y) with a x + b y = g = gcd(a, b) >= 0."""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        qq = a // b
        a, b = b, a - qq * b
        x0, x1 = x1, x0 - qq * x1
        y0, y1 = y1, y0 - qq * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0


def compose(f1, f2):
    """Dirichlet/Shanks composition of two primitive forms of the same discriminant."""
    a1, b1, c1 = f1
    a2, b2, c2 = f2
    D = disc(f1)
    assert disc(f2) == D
    s = (b1 + b2) // 2
    n = b2 - s
    # solve: g = gcd(a1, a2, s) with u a1 + v a2 + w s = g
    g1, u1, v1 = _xgcd(a1, a2)
    g, x, w = _xgcd(g1, s)
    u, v = x * u1, x * v1
    a1g, a2g = a1 // g, a2 // g
    A = a1g * a2g
    # B = b2 + 2 a2/g * (v (s - b2) - w c2)  (mod 2A)
    B = b2 + 2 * a2g * (v * (s - b2) - w * c2)
    B %= 2 * A
    C = (B * B - D) // (4 * A)
    return reduce_form((A, B, C))


def identity(D):
    """Principal form of discriminant D < 0."""
    if D % 4 == 0:
        return (1, 0, -D // 4)
    return (1, 1, (1 - D) // 4)


def power(f, e):
    D = disc(f)
    r = identity(D)
    base = reduce_form(f)
    while e:
        if e & 1:
            r = compose(r, base)
        base = compose(base, base)
        e >>= 1
    return r


def reduced_forms(D):
    """All reduced primitive forms of discriminant D < 0 (O(|D|) enumeration; diagnostic)."""
    assert D < 0 and D % 4 in (0, 1)
    out = []
    amax = isqrt(-D // 3) + 1
    for a in range(1, amax + 1):
        for b in range(-a + 1, a + 1):
            if (b - D) % 2:
                continue
            num = b * b - D
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if (a == c) and b < 0:
                continue
            if gcd(gcd(a, b), c) != 1:
                continue
            out.append((a, b, c))
    return out


def class_number(D):
    return len(reduced_forms(D))


def hurwitz_H(n):
    """Hurwitz class number H(n) for n >= 0 (n = 0 or 3 mod 4), as a Fraction.

    H(n) = number of SL2(Z)-classes of positive definite forms (not necessarily primitive) of
    discriminant -n, with forms equivalent to a(x^2+y^2) weighted 1/2 and a(x^2+xy+y^2)
    weighted 1/3; H(0) = -1/12. Enumeration O(n); diagnostic.
    """
    if n == 0:
        return Fraction(-1, 12)
    if n % 4 not in (0, 3):
        return Fraction(0)
    D = -n
    total = Fraction(0)
    amax = isqrt(n // 3) + 1
    for a in range(1, amax + 1):
        for b in range(-a + 1, a + 1):
            if (b - D) % 2:
                continue
            num = b * b - D
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a or (a == c and b < 0):
                continue
            if a == b == c:
                total += Fraction(1, 3)
            elif b == 0 and a == c:
                total += Fraction(1, 2)
            else:
                total += 1
    return total
