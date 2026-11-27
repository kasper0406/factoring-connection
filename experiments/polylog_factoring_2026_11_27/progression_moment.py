"""Evaluate the fixed-cap absolute moment using roots modulo the discriminant."""

from collections import Counter
from math import gcd, isqrt

from sympy.ntheory.residue_ntheory import sqrt_mod


def factor_public(value, work):
    factors = []
    prime = 2
    while prime * prime <= value:
        work["trial_divisors"] += 1
        exponent = 0
        while value % prime == 0:
            exponent += 1
            value //= prime
        if exponent:
            factors.append((prime, exponent))
        prime = 3 if prime == 2 else prime + 2
    if value > 1:
        factors.append((value, 1))
    return factors


def root_classes(value, prime, exponent):
    """Return a common step and disjoint unit-lifted root classes."""
    valuation = 0
    while valuation < exponent and value % prime == 0:
        value //= prime
        valuation += 1
    if valuation == exponent:
        return prime ** ((exponent + 1) // 2), [0]
    if valuation % 2:
        return 1, []
    power = exponent - valuation
    if prime == 2:
        modulus, roots = 2, [1]
        for _ in range(1, power):
            roots = [
                candidate
                for root in roots
                for candidate in (root, root + modulus)
                if (candidate * candidate - value) % (2 * modulus) == 0
            ]
            modulus *= 2
    else:
        roots = sqrt_mod(value, prime, all_roots=True)
        modulus = prime
        for _ in range(1, power):
            roots = [
                root + modulus * (
                    ((value - root * root) // modulus) * pow(2 * root, -1, prime)
                    % prime
                )
                for root in roots
            ]
            modulus *= prime
    scale = prime ** (valuation // 2)
    return modulus * scale, [root * scale for root in roots]


def classes(value, factors, work):
    modulus, residues = 1, [0]
    for prime, exponent in factors:
        work["prime_power_queries"] += 1
        step, roots = root_classes(value, prime, exponent)
        if not roots:
            return modulus, []
        inverse = pow(modulus, -1, step)
        residues = [
            residue + modulus * ((root - residue) * inverse % step)
            for residue in residues
            for root in roots
        ]
        modulus *= step
    return modulus, residues


def progression_points(n, delta, cap, work):
    scale, d = (2, delta // 4) if delta % 4 == 0 else (1, delta)
    factors = factor_public(d, work)
    for j in range(1, cap + 1):
        work["levels"] += 1
        height = 4 * n * j // (scale * scale)
        maximum = isqrt(height - d) if height >= d else 0
        step, residues = classes(height, factors, work)
        work["progressions"] += len(residues)
        for residue in residues:
            for x in range(residue or step, maximum + 1, step):
                work["square_tests"] += 1
                square = (height - x * x) // d
                y = isqrt(square)
                if y * y == square:
                    yield scale * x, y


def fixed_cap(n, work, points=progression_points):
    c = isqrt(4 * n)
    if 4 * n - c * c <= 4:
        c -= 1
    divisor = gcd(c, n)
    if divisor > 1:
        return divisor, None
    delta = 4 * n - c * c
    cap = isqrt(delta // 3)
    trace = sum(w for w, _ in points(n, delta, cap, work)) % n
    size = isqrt(cap)
    z = (trace - c * size * (size + 1) // 2) * pow(c, -1, n) % n
    divisor = gcd(z, n)
    if divisor == 1:
        square = z * z % n
        t = isqrt(square)
        assert t * t == square
        divisor = gcd(z - t, n)
    assert 1 < divisor < n and n % divisor == 0
    return divisor, trace


if __name__ == "__main__":
    from argparse import ArgumentParser

    cli = ArgumentParser()
    cli.add_argument("n", type=int)
    args = cli.parse_args()
    work = Counter()
    divisor, trace = fixed_cap(args.n, work)
    print({"divisor": divisor, "trace": trace, **work})
