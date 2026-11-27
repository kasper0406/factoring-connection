"""Split an odd modulus using parity counts on ordinary binary prefixes."""

from math import gcd


def split(n, prefix_parity):
    """Assume n is odd and has at least two distinct prime divisors.

    prefix_parity(start, size) counts roots of 1 in an aligned dyadic block mod 2.
    """
    stop = (n + 1) // 2
    start = 0
    size = 1 << (stop - 1).bit_length()
    while size > 1:
        size //= 2
        middle = start + size
        if middle >= stop:
            continue
        parity = prefix_parity(start, size) ^ int(start <= 1 < middle)
        if not parity:
            start = middle
    divisor = gcd(start - 1, n)
    if not (1 < divisor < n and start**2 % n == 1):
        raise ValueError("Decoded candidate failed verification")
    return divisor


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
