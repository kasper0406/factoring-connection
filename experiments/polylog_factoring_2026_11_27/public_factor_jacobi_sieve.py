"""Common public Jacobi tests for candidate balanced factors, with a scan benchmark."""

from math import gcd, isqrt
from random import Random, SystemRandom

if __package__:
    from .quadratic_density import _jacobi as jacobi
else:
    from quadratic_density import _jacobi as jacobi


def public_moduli(n, failure_bits=20, *, rng=None):
    """Return common moduli; uniform independent draws give failure <= 2^-failure_bits."""
    if n < 15 or n % 2 == 0 or failure_bits < 0:
        raise ValueError("expected an odd balanced semiprime and nonnegative failure bits")
    upper = isqrt(2*n)
    height = 100*(isqrt(2*n-1)+1)
    generator = SystemRandom() if rng is None else rng
    rounds = 2*(upper.bit_length()+failure_bits)
    return [1+4*n*generator.randrange(height) for _ in range(rounds)]


def interval_scan(n, moduli):
    """Enumerate the interval. The counts include the initial odd and nonsquare stages."""
    low, high = isqrt((n-1)//2)+1, isqrt(2*n)
    candidates = list(range(low | 1, high+1, 2))
    counts = [len(candidates)]
    candidates = [x for x in candidates if isqrt(x)**2 != x]
    counts.append(len(candidates))
    for d in moduli:
        candidates = [x for x in candidates if jacobi(d, x) == 1]
        counts.append(len(candidates))
    return candidates, counts


def self_check():
    periods = 0
    for n in (15, 35, 143, 323):
        for x in range(3, 80, 2):
            if gcd(x, n) != 1 or isqrt(x)**2 == x:
                continue
            values = [jacobi(1+4*t*n, x) for t in range(x)]
            assert values == [jacobi(x, 1+4*t*n) for t in range(x)]
            assert values.count(1) == values.count(-1) <= x//2
            periods += 1
    n = 2989601243
    moduli = public_moduli(n, rng=Random(20261127))
    survivors, counts = interval_scan(n, moduli)
    assert survivors == [42697, 70019] and survivors[0]*survivors[1] == n
    first_exact = counts[1:].index(2)
    print("Complete odd nonsquare character periods:", periods)
    print("N, odd candidate count, nonsquare candidate count:", n, *counts[:2])
    print("Draws and survivor counts:", [(i, counts[i+1]) for i in (1, 4, 8, 12, first_exact, len(moduli))])
    print("Survivors:", survivors)


if __name__ == "__main__":
    self_check()
