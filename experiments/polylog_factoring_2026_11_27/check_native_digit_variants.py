"""Finite checks for ordinary-binary and negabinary integer-polynomial extraction."""

from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import combinations
from math import gcd

from sympy import Poly, Symbol, primerange


X = Symbol('X')


def encode(value, base):
    """Canonical binary digits for a nonnegative integer, or negabinary for any integer."""
    assert base in (2, -2) and (base < 0 or value >= 0)
    word = 0
    i = 0
    while value:
        bit = value % 2
        word |= bit << i
        value = (value-bit)//base
        i += 1
    return word


def evaluate(word, base):
    value = 0
    for char in bin(word)[2:]:
        value = value*base + int(char)
    return value


def reverse_word(word):
    return int(bin(word)[2:][::-1], 2)


def reciprocal(poly):
    return Poly.from_list(poly.all_coeffs()[::-1], X, domain='ZZ')


@lru_cache(maxsize=2048)
def polynomial_factors(word):
    poly = Poly.from_list([int(bit) for bit in bin(word)[2:]], X, domain='ZZ')
    _, factors = poly.factor_list()
    return poly, tuple(f for f, _ in factors)


def check_encoding():
    total = 0
    for base in (2, -2):
        for z in range(-4096 if base < 0 else 0, 4097):
            word = encode(z, base)
            assert evaluate(word, base) == z
            if word:
                degree = word.bit_length()-1
                if base == -2:
                    assert (1 << degree) <= 3*abs(z) <= 4 << degree
                    assert (z < 0) == bool(degree % 2)
            total += 1
    return total


def check_orientations():
    groups_checked = common_words = 0
    for b in (4, 5, 6):
        radius = 1 << (b-1)
        primes = tuple(map(int, primerange(radius, 2*radius)))
        for degree in (2*b-2, 2*b-1, 2*b):
            # The centered interval [s, degree-s] has at most b-1 positions.
            s = (degree-b+3)//2
            assert degree-2*s+1 <= b-1
            pairs = [(i, degree-i) for i in range(s, (degree+1)//2)]
            free = sum((1 << i) | (1 << j) for i, j in pairs)
            groups = defaultdict(list)
            for word in range((1 << degree)+1, 2 << degree, 2):
                sums = tuple(((word >> i) & 1) + ((word >> j) & 1) for i, j in pairs)
                groups[word & ~free, sums].append(word)
            weights = Fraction(0)
            for (_, sums), words in groups.items():
                j = sums.count(1)
                assert len(words) == 2**j
                weights += Fraction(len(words), 2**j)
                for base in (2, -2):
                    values = [evaluate(word, base) for word in words]
                    mirrors = [evaluate(reverse_word(word), base) for word in words]
                    totals = {a+c for a, c in zip(values, mirrors)}
                    assert len(totals) == 1 and 0 not in totals
                    h = totals.pop()
                    assert (max(values)-min(values)) % (2**s) == 0
                    assert (max(values)-min(values))//(2**s) < radius
                    assert len(set(values)) == len(words)
                    passing = set()
                    candidate_count = 0
                    for r in primes:
                        marked = {idx for idx, val in enumerate(values) if val % r == 0}
                        assert len(marked) <= 1
                        if h % r == 0:
                            candidate_count += 1
                            assert all(mirrors[idx] % r == 0 for idx in marked)
                            passing.update(marked)
                    common_words += len(passing)
                    assert len(passing) <= candidate_count
                    groups_checked += 1
            assert weights / 2**(degree-1) == Fraction(3, 4)**len(pairs)
    return groups_checked, common_words


def check_factorizations():
    counts = {base: Counter() for base in (2, -2)}
    witnesses = {}
    for b in (4, 5, 6):
        radius = 1 << (b-1)
        for p, q in combinations(map(int, primerange(radius, 2*radius)), 2):
            n = p*q
            for base in (2, -2):
                seen = set()
                for multiplier in (1, -1, 2, -2, 3, -3, 5, -5, 6, -6):
                    z = multiplier*n
                    if base == 2:
                        z = abs(z)
                    original = encode(z, base)
                    word = original
                    while word % 2 == 0:
                        word >>= 1
                        z //= base
                    assert encode(z, base) == word and z % n == 0
                    if word in seen:
                        continue
                    seen.add(word)
                    counts[base]['evaluations'] += 1
                    poly, factors = polynomial_factors(word)
                    mirror = int(reciprocal(poly).eval(base))
                    assert poly.eval(base) == z
                    assert (int(mirror) > 0) == (z > 0)
                    direct = False
                    for f in factors:
                        divisor = gcd(n, int(f.eval(base)))
                        if divisor == 1 or divisor == n:
                            continue
                        direct = True
                        g = poly.exquo(f)
                        if f == reciprocal(f):
                            assert mirror % divisor == 0
                        elif g == reciprocal(g):
                            assert mirror % (n//divisor) == 0
                        else:
                            twin = f*reciprocal(g)
                            assert set(twin.all_coeffs()) <= {0, 1}
                            assert twin != poly and twin != reciprocal(poly)
                            assert twin*reciprocal(twin) == poly*reciprocal(poly)
                            counts[base]['nonreciprocal factor splits'] += 1
                    if direct:
                        counts[base]['successful evaluations'] += 1
                        if gcd(n, mirror) > 1:
                            counts[base]['common prime'] += 1
                        else:
                            counts[base]['twin only'] += 1
                            witnesses.setdefault(base, (n, z//n, word))
    return counts, witnesses


def check_twin_fixtures():
    # The second example was selected separately to exercise this less frequent
    # branch in negative base. It is excluded from the exhaustive table.
    for base, p, q, multiplier in ((2, 19, 29, 5), (-2, 37, 61, -17)):
        word = encode(multiplier*p*q, base)
        poly, factors = polynomial_factors(word)
        assert gcd(p*q, int(reciprocal(poly).eval(base))) == 1
        for target in (p, q):
            f = next(f for f in factors if gcd(int(f.eval(base)), p*q) == target)
            g = poly.exquo(f)
            twin = f*reciprocal(g)
            assert set(twin.all_coeffs()) <= {0, 1}
            assert twin != poly and twin != reciprocal(poly)
            assert twin*reciprocal(twin) == poly*reciprocal(poly)


if __name__ == '__main__':
    print(f'{check_encoding()} canonical encodings checked.')
    groups, hits = check_orientations()
    print(f'{groups} base-specific orientation groups checked, with {hits} common-prime words.')
    results, witnesses = check_factorizations()
    check_twin_fixtures()
    for base, counts in results.items():
        print(f'Base {base}: {dict(counts)}, first twin-only fixture {witnesses.get(base)}')
