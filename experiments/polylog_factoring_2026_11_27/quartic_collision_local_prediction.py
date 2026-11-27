"""Exact best single exponent given only the semiprime residue modulo 4*ell."""

from collections import Counter
from fractions import Fraction
from math import gcd


def local_prediction(n, ell):
    """Select an optimal exponent for the conditional uniform residue model."""
    modulus = 4 * ell
    if gcd(n, modulus) != 1 or n % ell == 1:
        raise ValueError("Expected an eligible odd product and auxiliary prime")
    for slope in range(min(8, ell)):
        if slope % ell in (2, 1, ell - 1):
            continue
        targets = []
        for sign in (1, -1):
            residue = (slope + sign) * pow((slope - 2) % ell, -1, ell) % ell
            targets.append(residue + ell * ((sign - residue) * pow(ell, -1, 4) % 4))
        partners = {n * pow(target, -1, modulus) % modulus for target in targets}
        if partners.isdisjoint(targets):
            return slope * pow(5, -1, ell) % ell
    raise ArithmeticError("No disjoint target pair")


def optimal_local_scores(ell):
    """Count ordered unit-class pairs, conditionally on each eligible product."""
    modulus = 4 * ell
    units = [x for x in range(1, modulus) if gcd(x, modulus) == 1]
    exponents = {}
    for x in units:
        if x % ell != 1:
            character = 1 if x % 4 == 1 else -1
            exponents[x] = (2*x + character) * pow(5*(x-1), -1, ell) % ell

    best = {}
    for product in units:
        if product % ell == 1:
            continue
        scores = Counter()
        for left in units:
            right = product * pow(left, -1, modulus) % modulus
            hits = {exponents[x] for x in (left, right) if x in exponents}
            scores.update(hits)
        best[product] = max(scores.values())
    return Counter(best.values()), Fraction(sum(best.values()), len(best) * len(units))


if __name__ == "__main__":
    for prime in (7, 11, 17, 31, 67, 127, 257):
        histogram, success = optimal_local_scores(prime)
        print(prime, sorted(histogram.items()), success)
