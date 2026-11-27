"""Split an odd modulus using parity counts on ordinary binary prefixes."""

from math import gcd, isqrt

from bounded_quadratic_roots import quadratic_roots


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
        parity = prefix_parity(start, size) ^ int(start <= 1 < middle)
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
