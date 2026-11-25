"""Compare a norm form, an eta product, and independent local elliptic counts."""
from collections import Counter
from itertools import combinations
from math import isqrt

from cypari2 import Pari


def norm_representations(n):
    limit = isqrt(4 * n // 23)
    count = 0
    for t in range(-limit, limit + 1):
        disc = 4 * n - 23 * t * t
        root = isqrt(disc)
        if root * root == disc and (root - t) % 2 == 0:
            count += 1 if root == 0 else 2
    return count


def pentagonal(limit):
    bound = isqrt(2 * limit) + 1
    return [k * (3 * k - 1) // 2 for k in range(-bound, bound + 1)
            if k * (3 * k - 1) // 2 <= limit]


def eta_parity_bits(limit):
    first_factor = sum(1 << i for i in pentagonal(limit))
    product = 0
    for j in pentagonal((limit - 1) // 23):
        product ^= first_factor << (1 + 23 * j)
    return product


def main():
    pari = Pari()
    E = pari.ellinit([0, 0, 0, -1, 1])
    local_orders = {r: r + 1 - int(pari.ellap(E, r)) for r in range(3, 312)
                    if r != 23 and pari.isprime(r)}
    bits = eta_parity_bits(311 ** 2)
    totals = Counter()
    for p, q in combinations(local_orders, 2):
        n = p * q
        correct = (local_orders[p] * local_orders[q]) % 2
        count = norm_representations(n)
        assert count in (0, 4, 8)
        assert correct == count // 4 % 2 == (bits >> n) & 1
        totals[count] += 1
    print(sum(totals.values()), sorted(totals.items()))


if __name__ == "__main__":
    main()
