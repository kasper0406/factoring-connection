"""Compare balanced wide-window queries against the existing small-multiplier baseline."""

from math import gcd, isqrt
from random import Random

from sympy import nextprime
from sympy.ntheory.modular import crt


def window_answer(n, left, right):
    """Assume n is a product of two distinct odd primes larger than 32."""
    for multiplier in range(1, 33):
        target = 4 * multiplier * n
        first = isqrt(target - 1) + 1
        for center in range(first, first + 8):
            gap = center * center - target
            delta = isqrt(gap)
            if delta * delta != gap:
                continue
            p = gcd(center + delta, n)
            if 1 < p < n:
                q = n // p
                root = (1 - 2 * p * pow(p, -1, q)) % n
                return sum(left <= x < right for x in (1, n - 1, root, n - root)) % 2
    return None


def benchmark():
    rng = Random(20261127)
    for bits in (16, 32, 48):
        for planted in (True, False):
            solved = odd = 0
            for _ in range(100):
                p = int(nextprime(rng.randrange(2**(bits-1), 2**bits)))
                center = 3 * p // 2 if planted else rng.randrange(135*p//100, 165*p//100)
                q = int(nextprime(center))
                n = p * q
                left, right = n // 4, n // 2
                answer = window_answer(n, left, right)
                if answer is None:
                    continue
                roots = [int(crt((p, q), (a, b))[0]) for a in (-1, 1) for b in (-1, 1)]
                assert answer == sum(left <= r < right for r in roots) % 2
                solved += 1
                odd += answer
            print(bits, 'planted' if planted else 'control', solved, odd)


if __name__ == '__main__':
    benchmark()
