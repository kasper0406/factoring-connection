"""Algebraic intersection of two matrix loops on the compact curve X_0(n)."""

from math import gcd


def _word(matrix, n):
    """Write a Gamma_0(n) matrix as powers of T and copies of S, in PSL_2."""
    a, b, c, d = matrix
    if a*d - b*c != 1 or c % n:
        raise ValueError("both matrices must have determinant one and lower left divisible by n")
    while c:
        # Nearest-integer division ensures the next lower-left entry is at
        # most half as large. A one-sided division can take a long time.
        power = (2*a + c) // (2*c)
        if power:
            yield "T", power
        yield "S", 1
        a, b, c, d = c, d, power*c - a, power*d - b
    if b:
        yield "T", b*d  # the remaining upper triangular matrix, up to sign


def _crossing_runs(matrix, n):
    """Yield the prefix row, inclusive turn interval, and signed crossing of [0,infinity]."""
    c, d = 0, 1
    for generator, power in _word(matrix, n):
        if generator == "S":
            yield c, d, 0, 0, 1
            c, d = d, -c % n
        else:
            # A full orbit bounds the filled cusp disk and pairs to zero.
            width = n // gcd(c*c, n)
            remainder = abs(power) % width
            if remainder:
                if power > 0:
                    yield c, d, 1, remainder, -1
                else:
                    yield c, d, 1 - remainder, 0, 1
            d = (d + power*c) % n


def _primal_edges(matrix, n):
    # T fixes infinity. An S step in the cusp path infinity -> g*infinity
    # contributes the negative of the prefix edge [0,infinity].
    c, d = 0, 1
    for generator, power in _word(matrix, n):
        if generator == "S":
            yield c, d
            c, d = d, -c % n
        else:
            d = (d + power*c) % n


def unit_chord_reduction(n, g):
    """Return (proper divisor, None) or (None, unit chord indices) for the loop g.

    The unit loop indexed by k runs from 0 to -1/k, with its endpoints
    identified on X_0(n). The sum of the returned loops represents [g].
    """
    if n < 1:
        raise ValueError("n must be positive")
    chords = []
    for c, d in _primal_edges(g, n):
        for value in (c, d):
            divisor = gcd(value, n)
            if 1 < divisor < n:
                return divisor, None
        if n > 1 and c % n and d % n:
            chords.append(d * pow(c, -1, n) % n)
    return None, tuple(chords)


def _edge_matches(edge, c, d, low, high, n):
    """Count j in [low,high] with edge projectively equal to (c,d+j*c) modulo n."""
    x, y = edge
    coefficient = x*c
    rhs = y*c - x*d
    common = gcd(coefficient, n)
    if rhs % common:
        return 0
    modulus = n // common
    if modulus == 1:
        return high - low + 1
    first = (rhs // common) * pow(coefficient // common, -1, modulus) % modulus
    return (high - first) // modulus - (low - 1 - first) // modulus


def level_intersection(n, g, h):
    """Return [g]·[h] on X_0(n) for row-major SL_2(Z) tuples from Gamma_0(n).

    Orient the upper half-plane by dx wedge dy. The two loops run from z
    to g*z and from z to h*z; the choice of z leaves their classes unchanged.
    """
    if n < 1:
        raise ValueError("n must be positive")
    primal = list(_primal_edges(g, n))
    total = 0
    for c, d, low, high, sign in _crossing_runs(h, n):
        for x, y in primal:
            forward = _edge_matches((x, y), c, d, low, high, n)
            reverse = _edge_matches((y, -x), c, d, low, high, n)
            total += sign * (reverse - forward)
    return total
