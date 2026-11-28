"""Exhaustive checks of the native binary polynomial and reversal arguments."""

from collections import defaultdict
from itertools import combinations
from math import gcd

from sympy import Poly, Symbol, primerange


X = Symbol("X")


def reverse_word(word):
    return int(f"{word:b}"[::-1], 2)


def reciprocal(poly):
    return Poly.from_list(poly.all_coeffs()[::-1], X, domain="ZZ")


def check_semiprimes(bits):
    primes = list(primerange(1 << (bits - 1), 1 << bits))
    pairs = reducible = reversal_hits = twins_only = 0
    for p, q in combinations(primes, 2):
        pairs += 1
        n = p * q
        poly = Poly.from_list([int(bit) for bit in f"{n:b}"], X, domain="ZZ")
        assert reciprocal(poly).eval(2) == reverse_word(n)
        content, factors = poly.factor_list()
        assert content == 1
        assert all(factor.eval(2) > abs(factor.eval(1)) >= 1 for factor, _ in factors)
        if len(factors) == 1 and factors[0][1] == 1:
            continue
        reducible += 1
        assert len(factors) == 2 and all(power == 1 for _, power in factors)
        left, right = [factor for factor, _ in factors]
        assert {left.eval(2), right.eval(2)} == {p, q}
        assert left.TC() == right.TC() == 1
        twin = left * reciprocal(right)
        assert set(twin.all_coeffs()) <= {0, 1}
        assert twin * reciprocal(twin) == poly * reciprocal(poly)
        trivial = twin == poly or twin == reciprocal(poly)
        assert trivial == (left == reciprocal(left) or right == reciprocal(right))
        has_common_prime = gcd(n, reverse_word(n)) > 1
        assert has_common_prime or not trivial
        reversal_hits += has_common_prime
        twins_only += not has_common_prime
    return pairs, reducible, reversal_hits, twins_only


def check_reversal_fibres(bits, degree):
    """Check each conditional word family and each eligible prime separately."""
    length = degree + 1
    middle = bits - 1 - ((length - (bits - 1)) % 2)
    offset = (length - middle) // 2
    mask = ((1 << middle) - 1) << offset
    primes = tuple(primerange(1 << (bits - 1), 1 << bits))
    groups = defaultdict(list)
    for word in range((1 << degree) + 1, 1 << length, 2):
        sums = tuple(((word >> (offset + j)) & 1) + ((word >> (degree - offset - j)) & 1)
                     for j in range(middle // 2))
        centre = (word >> (length // 2)) & 1 if middle % 2 else 0
        groups[word & ~mask, sums, centre].append(word)
    prime_fibres = 0
    for (_, sums, _), words in groups.items():
        assert len(words) == 1 << sums.count(1)
        total = words[0] + reverse_word(words[0])
        assert all(word + reverse_word(word) == total for word in words)
        for prime in primes:
            if total % prime == 0:
                prime_fibres += 1
                assert len({word % prime for word in words}) == len(words)
    return len(groups), prime_fibres


if __name__ == "__main__":
    counts = [check_semiprimes(bits) for bits in range(3, 9)]
    print("Pairs, reducible, shared reversal prime, twin without a shared prime:",
          tuple(map(sum, zip(*counts))))
    fibres = [check_reversal_fibres(bits, degree)
              for bits in range(3, 9) for degree in (2 * bits - 2, 2 * bits - 1)]
    print("Conditional word families, candidate-prime checks:", tuple(map(sum, zip(*fibres))))
