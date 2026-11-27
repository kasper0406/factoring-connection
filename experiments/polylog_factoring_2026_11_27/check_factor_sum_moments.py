"""Exact binomial rejection moments from two CRT residue lists per subset."""

from bisect import bisect_left
from itertools import combinations
from math import comb

from sympy import Matrix


def rejection_classes(primes, masks):
    modulus, residues = 1, [0]
    for p in primes:
        inverse = pow(modulus, -1, p)
        residues = [a + modulus * ((b - a) * inverse % p)
                    for a in residues for b in masks[p]]
        modulus *= p
    return modulus, residues


def intersection_size(primes, masks, lo, hi):
    left, residues = rejection_classes(primes[::2], masks)
    right, targets = rejection_classes(primes[1::2], masks)
    inverse = pow(left, -1, right)
    positions = sorted(b * inverse % right for b in targets)

    def prefix(end):
        periods, rem = divmod(end, right)
        return periods * len(positions) + bisect_left(positions, rem)

    count = 0
    for a in residues:
        # t = a + left*j. Shift j so its allowed residues are positions.
        shift = a * inverse % right
        start = (lo - 1 - a) // left + 1 + shift
        stop = (hi - a) // left + 1 + shift
        count += prefix(stop) - prefix(start)
    return count


def main():
    n, lo, hi = 2989601243, 54678, 56335
    primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
    masks = {p: {t for t in range(p)
                 if pow((t*t - n) % p, (p-1)//2, p) == p-1} for p in primes}
    moments = [sum(intersection_size(s, masks, lo, hi)
                   for s in combinations(primes, k)) for k in range(11)]
    # Independent direct control on the small benchmark only.
    counts = [sum(t % p in masks[p] for p in primes) for t in range(lo, hi+1)]
    assert moments == [sum(comb(r, k) for r in counts) for k in range(11)]
    print("Binomial moments", moments)
    for degree in range(1, 6):
        # Convert binomial to ordinary moments by forward differences.
        def power_moment(k):
            differences = [r**k for r in range(k+1)]
            value = 0
            for coefficient in moments[:k+1]:
                value += coefficient * differences[0]
                differences = [b-a for a, b in zip(differences, differences[1:])]
            return value

        gram = Matrix(degree+1, degree+1, lambda i, j: power_moment(int(i+j)))
        bound = 1 / gram.inv()[0, 0]
        print(degree, float(bound))


if __name__ == "__main__":
    main()
