"""Enumerate x*x+y*y=N in the sector 0 <= x <= y using rational tangents."""

from math import isqrt


def tangent_direction(numerator, denominator, bound):
    """Return the last continued-fraction convergent with denominator <= bound."""
    old_p, p = 0, 1
    old_q, q = 1, 0
    while denominator:
        quotient, remainder = divmod(numerator, denominator)
        next_q = quotient * q + old_q
        if next_q > bound:
            break
        old_p, p = p, quotient * p + old_p
        old_q, q = q, next_q
        numerator, denominator = denominator, remainder
    return p, q


def circle_points(n):
    if n < 1:
        raise ValueError("expected a positive integer")
    # This power of two is within a factor of two of the sixth root.
    width = 1 << ((n.bit_length() - 1) // 6)
    last = isqrt(n // 2)
    for start in range(0, last + 1, width):
        radicand = n - start * start
        p, q = tangent_direction(start, isqrt(radicand), width)
        center = p * start + isqrt(q * q * radicand)
        norm = p * p + q * q
        for level in range(center - 7, center + 8):
            delta = norm * n - level * level
            if delta < 0:
                continue
            root = isqrt(delta)
            if root * root != delta:
                continue
            for signed in {root, -root}:
                x, rx = divmod(p * level + q * signed, norm)
                y, ry = divmod(q * level - p * signed, norm)
                if not (rx or ry) and start <= x < start + width and 0 <= x <= y:
                    yield x, y


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    for point in circle_points(parser.parse_args().N):
        print(*point)
