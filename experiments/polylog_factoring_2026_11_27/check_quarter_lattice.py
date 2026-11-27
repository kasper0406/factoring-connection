"""Small fixed-dimension diagnostic for the final translated root search."""
from math import gcd, isqrt
from random import Random

from sympy import Matrix, Poly, nextprime, symbols

T = symbols("T")


def candidates(n, center, radius):
    """Try a six-dimensional lattice using public interval data alone."""
    f = Poly(T**2 + (2 * center % n) * T + ((center**2 - 1) % n), T)
    basis = []
    for power in range(3):
        for shift in range(2):
            g = Poly(T**shift, T) * f**power * n**(2 - power)
            basis.append([g.nth(i) * radius**i for i in range(6)])
    for row in Matrix(basis).lll().tolist():
        g = Poly(sum(entry // radius**i * T**i for i, entry in enumerate(row)), T)
        for offset in g.ground_roots():
            if offset.is_Integer and abs(offset) <= radius:
                a = center + int(offset)
                if a * a % n == 1 and 1 < gcd(a - 1, n) < n:
                    yield a


def check():
    rng = Random(421)
    hits = [0, 0, 0]
    trials = 0
    while trials < 24:
        p, q = (int(nextprime(rng.randrange(1024, 2047))) for _ in range(2))
        if p == q:
            continue
        n = p * q
        root = 1 + p * ((-2 * pow(p, -1, q)) % q)
        root = min(root, n - root)
        for i, divisor in enumerate((32, 8, 2)):
            radius = isqrt(n) // divisor
            center = root + rng.randint(-radius, radius)
            hits[i] += root in set(candidates(n, center, radius))
        trials += 1
    print("24 planted intervals per radius; recoveries:", hits)


if __name__ == "__main__":
    check()
