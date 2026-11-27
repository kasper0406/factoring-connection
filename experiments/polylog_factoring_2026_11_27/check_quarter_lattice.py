"""Small fixed-dimension diagnostic for the final translated root search."""
from math import gcd, isqrt
from random import Random

from sympy import Matrix, Poly, nextprime, symbols

T = symbols("T")


def candidates(n, center, radius, *, certified_only=False, max_level=0):
    """Try a six-dimensional lattice using public interval data alone."""
    f = Poly(T**2 + (2 * center % n) * T + ((center**2 - 1) % n), T)
    basis = []
    for power in range(3):
        for shift in range(2):
            g = Poly(T**shift, T) * f**power * n**(2 - power)
            basis.append([g.nth(i) * radius**i for i in range(6)])
    for row in Matrix(basis).lll().tolist():
        level = sum(abs(entry) for entry in row) // n**2 if certified_only else 0
        if level > max_level or not any(row[1:]):
            continue
        g = Poly(sum(entry // radius**i * T**i for i, entry in enumerate(row)), T)
        for multiple in range(-level, level + 1):
            for offset in (g - multiple * n**2).ground_roots():
                if offset.is_Integer and abs(offset) <= radius:
                    a = center + int(offset)
                    if a * a % n == 1 and 1 < gcd(a - 1, n) < n:
                        yield a
        if certified_only:
            return
    if certified_only:
        raise ValueError("No completeness certificate")


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


def check_certificates(max_level=0):
    rng = Random(724)
    counts = [[0, 0] for _ in range(3)]
    for trial in range(48):
        p = int(nextprime(rng.randrange(1024, 1500)))
        q = int(nextprime(rng.randrange(1500, 2047)))
        n = p * q
        u = 1 + p * ((-2 * pow(p, -1, q)) % q)
        for i, divisor in enumerate((32, 8, 2)):
            radius = isqrt(n) // divisor
            center = rng.randrange(radius + 2, n - radius - 2)
            if trial % 2:
                center = u + rng.randint(-radius, radius)
            try:
                roots = set(candidates(n, center, radius, certified_only=True, max_level=max_level))
            except ValueError as error:
                if str(error) != "No completeness certificate":
                    raise
                continue
            expected = {a for a in range(center-radius, center+radius+1)
                        if a*a % n == 1 and 1 < gcd(a-1, n) < n}
            assert roots == expected
            counts[i][bool(roots)] += 1
    print("Maximum level", max_level, "certified empty and nonempty counts:", counts)


if __name__ == "__main__":
    check()
    check_certificates()
    check_certificates(max_level=8)
