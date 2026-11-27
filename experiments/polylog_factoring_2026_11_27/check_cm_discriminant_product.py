"""Decode a supplied square class of (4p-1)(4q-1) for balanced factors."""

from math import isqrt, prod


def decode(n, square_class):
    if n < 2 or square_class < 1:
        raise ValueError("Expected n >= 2 and a positive square-class representative")
    top = isqrt(16 * n // square_class)
    for z in (top, top - 1):
        if z < 1:
            continue
        total, remainder = divmod(16 * n + 1 - square_class * z * z, 4)
        if remainder or total < 1:
            continue
        delta = total * total - 4 * n
        if delta < 0:
            continue
        gap = isqrt(delta)
        if gap * gap != delta or (total - gap) % 2:
            continue
        p, q = (total - gap) // 2, (total + gap) // 2
        if 1 < p <= q <= 2 * p and p * q == n:
            return p, q
    return None


if __name__ == "__main__":
    from sympy import factorint, primerange

    def kernel(value):
        return prod(int(p) for p, e in factorint(value).items() if e % 2)

    primes = list(primerange(3, 500))
    checked = 0
    for i, p in enumerate(primes):
        for q in primes[i + 1:]:
            if q > 2 * p:
                break
            h = kernel((4 * p - 1) * (4 * q - 1))
            assert decode(p*q, h) == (p, q)
            checked += 1
    assert decode(703, 1) == (19, 37)
    assert decode(703, 2) is None
    print(f"Verified {checked} balanced prime pairs")
