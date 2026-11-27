"""Split an odd modulus using parity counts on ordinary binary prefixes."""

from math import gcd, isqrt

from sympy import primerange

from bounded_quadratic_roots import quadratic_roots


def progression_parity(n, intervals, length):
    """Count roots of one using the largest easily extracted prime power."""
    stride = 1
    for prime in primerange(3, n.bit_length() + 1):
        rest = n
        while rest % prime == 0:
            rest //= prime
        stride = max(stride, n // rest)
    modulus = n // stride
    if stride == 1 or length > stride * isqrt(modulus):
        return None
    inverse = pow(stride, -1, modulus) if modulus > 1 else 0
    block = max(1, isqrt(modulus) // 8) + 1
    parity = 0
    for left, right in intervals:
        for sign in (-1, 1):
            first = left + (sign - left) % stride
            count = max(0, (right - first + stride - 1) // stride)
            if modulus == 1:
                parity ^= count % 2
                continue
            for offset in range(0, count, block):
                base = (first * inverse + offset) % modulus
                span = min(block, count - offset) - 1
                roots = quadratic_roots(modulus, 2*base, base**2-inverse**2, span)
                parity ^= len(roots) % 2
    return parity


def reflected_parity(n, y, start, width, half_parity):
    """For a unit square on odd n, return cyclic-window parity or None.

    half_parity is one for a prime power, zero for multiple prime factors.
    Use reflection and short-root lattices, optionally on prime-power progressions.
    """
    if n <= 1 or n % 2 == 0 or not 0 <= width <= n or half_parity not in (0, 1):
        raise ValueError("Expected an odd modulus, valid width and parity bit")
    stop = (n+1)//2
    end = start + width
    boundaries = {1, stop}
    for point in (start, end, 1-end, 1-start):
        point %= n
        if 1 < point < stop:
            boundaries.add(point)
    groups = [[], []]
    points = sorted(boundaries)
    for left, right in zip(points, points[1:]):
        inside = (left-start) % n < width
        reflected = (-left-start) % n < width
        groups[inside ^ reflected].append((left, right))
    lengths = [sum(right-left for left, right in group) for group in groups]
    side = int(lengths[1] <= lengths[0])
    cutoff = isqrt(n)
    if lengths[side] > cutoff:
        if y % n != 1:
            return None
        bit = progression_parity(n, groups[side], lengths[side])
        return None if bit is None else bit ^ (0 if side else half_parity)
    parity = 0 if side else half_parity
    block = max(1, cutoff//8) + 1
    for left, right in groups[side]:
        for origin in range(left, right, block):
            span = min(block, right-origin) - 1
            parity ^= len(quadratic_roots(n, 2*origin, origin**2-y, span)) % 2
    return parity


def split(n, prefix_parity):
    """Assume n is odd and has at least two distinct prime divisors.

    prefix_parity(start, size) counts roots of 1 in an aligned dyadic block mod 2.
    """
    stop = (n + 1) // 2
    start = 0
    size = 1 << (stop - 1).bit_length()
    cutoff = isqrt(n)
    while min(size, stop-start) > cutoff:
        size //= 2
        middle = start + size
        if middle >= stop:
            continue
        parity = reflected_parity(n, 1, start, size, half_parity=0)
        if parity is None:
            parity = prefix_parity(start, size)
        parity ^= int(start <= 1 < middle)
        if not parity:
            start = middle
    end = min(start + size, stop)
    block = max(1, cutoff // 8) + 1
    for origin in range(start, end, block):
        span = min(block, end-origin) - 1
        for offset in quadratic_roots(n, 2*origin, origin**2-1, span):
            root = origin + offset
            divisor = gcd(root - 1, n)
            if 1 < divisor < n and root**2 % n == 1:
                return divisor
    raise ValueError("No verified mixed root in terminal interval")


if __name__ == "__main__":
    from sympy import factorint

    cases = 0
    for n in range(3, 2000, 2):
        if len(factorint(n)) < 2:
            continue
        def oracle(start, size):
            assert start % size == 0 and start + size <= (n+1)//2
            return sum(x*x % n == 1 for x in range(start, start+size)) % 2
        assert n % split(n, oracle) == 0
        cases += 1
    print(f"Verified splits for {cases} odd moduli.")
