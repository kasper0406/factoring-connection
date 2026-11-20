"""Blind CM derivative experiment. Only the scoring loop knows the factors."""

from collections import Counter
from itertools import combinations
from math import gcd, isqrt


def derivative(coefficients, has_y, modulus):
    # Normal form is y^has_y times a polynomial in x, with y² = x³ - x.
    result = [0] * (len(coefficients) + (2 if has_y else 0))
    for k, coefficient in enumerate(coefficients):
        if has_y:
            result[k + 2] += (2 * k + 3) * coefficient
            result[k] -= (2 * k + 1) * coefficient
        elif k:
            result[k - 1] += 2 * k * coefficient
    result = [c % modulus for c in result]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def probes(n):
    coefficients = [0, 1]
    results = []
    for order in range(1, 256):
        coefficients = derivative(coefficients, (order - 1) % 2, n * n)
        if order in (31, 127, 255):
            modulus = n if order < 255 else n * n
            results.append(gcd(modulus, *coefficients))
    return tuple(results)


primes = [p for p in range(67, 128) if all(p % d for d in range(2, isqrt(p) + 1))]
counts = Counter()
for p, q in combinations(primes, 2):
    inert = [r for r in (p, q) if r % 4 == 3]
    product = (p if p in inert else 1) * (q if q in inert else 1)
    values = probes(p * q)
    assert values == (1, product, product**2), (p, q, values)
    counts[len(inert)] += 1
print(f"{sum(counts.values())} pairs, all predictions matched")
print(f"0 inert: {counts[0]}, 1 inert: {counts[1]}, 2 inert: {counts[2]}")
print(f"67 * 73: {probes(67 * 73)}")
