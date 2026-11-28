"""Exact controls for binary convolution, eligible divisors, and noisy masks."""

from fractions import Fraction
from itertools import combinations, product
from math import comb, prod

from sympy import primerange

if __package__:
    from .carry_parity import (
        _degree_counts, _divisors, binary_product, factor_from_carry_parity,
        polynomial_factors,
    )
else:
    from carry_parity import (
        _degree_counts, _divisors, binary_product, factor_from_carry_parity,
        polynomial_factors,
    )


def carries(p, q):
    """Compute incoming carries directly from schoolbook multiplication."""
    width = (p * q).bit_length()
    incoming = [0] * (width + 1)
    for j in range(width):
        raw = sum(((p >> i) & 1) * ((q >> (j - i)) & 1)
                  for i in range(j + 1))
        outgoing, digit = divmod(raw + incoming[j], 2)
        assert digit == (p * q >> j) & 1
        incoming[j + 1] = outgoing
    assert incoming[-1] == 0
    return incoming[:-1]


def carry_mask(values):
    return sum((value % 2) << j for j, value in enumerate(values))


def divides(polynomial, divisor):
    while polynomial.bit_length() >= divisor.bit_length():
        polynomial ^= divisor << (polynomial.bit_length() - divisor.bit_length())
    return polynomial == 0


def check_divisors():
    trials = 0
    for degree in range(1, 9):
        total = 0
        population = 1 << (degree - 1)
        for bits in range((1 << degree) + 1, 2 << degree, 2):
            factors = polynomial_factors(bits)
            counts = _degree_counts(factors, degree)
            direct = {g for g in range(1, bits + 1, 2) if divides(bits, g)}
            decoded = set()
            for d in range(degree + 1):
                if counts[0][d]:
                    values = list(_divisors(factors, counts, d))
                    assert counts[0][d] == len(values)
                    decoded.update(values)
            assert direct == decoded and len(direct) == prod(e + 1 for _, e in factors)
            total += len(direct)
            trials += 1
        assert Fraction(total, population) == Fraction(degree + 3, 2)
    return trials


def check_odd_carry_bound():
    """First m free columns stochastically dominate Binomial(m, 1/4)."""
    m = 5
    histogram = [0] * (m + 1)
    odd_words = range(1, 1 << (m + 1), 2)
    for p, q in product(odd_words, repeat=2):
        incoming = carries(p, q)
        mask = carry_mask(incoming)
        assert mask == p * q ^ binary_product(p, q)
        observed = (mask >> 2) & ((1 << m) - 1)
        histogram[observed.bit_count()] += 1
        if not mask:
            assert not any(incoming)
    assert sum(histogram) == 4**m
    for k in range(m + 1):
        assert sum(histogram[:k + 1]) <= sum(comb(m, t) * 3**(m - t)
                                             for t in range(k + 1))
    return histogram


def check_decoder():
    pairs = list(combinations(primerange(64, 128), 2))
    true_covered = [0, 0, 0]
    public_splits = [0, 0, 0]
    for p, q in pairs:
        n = p * q
        mask = carry_mask(carries(p, q))
        assert factor_from_carry_parity(n, (mask,)) == p
        assert factor_from_carry_parity(n, (mask ^ 1,), errors=1) == p
        for errors in range(3):
            result = factor_from_carry_parity(n, errors=errors)
            assert result in (None, p)
            public_splits[errors] += result is not None
            if mask.bit_count() <= errors:
                assert result == p
                true_covered[errors] += 1
    assert true_covered == [0, 1, 9] and public_splits == [1, 16, 60]
    p, q = 881, 967
    n = p * q
    incoming = carries(p, q)
    assert [j for j, x in enumerate(incoming) if x % 2] == [6, 19]
    assert sum(x > 0 and x % 2 == 0 for x in incoming) == 12
    stats = {}
    assert factor_from_carry_parity(n, errors=2, stats=stats) == p
    cap_stats = {}
    mask = carry_mask(incoming)
    assert factor_from_carry_parity(n, (mask, mask), cap=0, stats=cap_stats) is None
    assert cap_stats == {"masks": 1, "skipped": 1, "divisors": 0}
    return len(pairs), true_covered, public_splits, stats


if __name__ == "__main__":
    print("Odd polynomials checked by independent long division:", check_divisors())
    print("Odd carry counts across the first five free columns:", check_odd_carry_bound())
    pairs, true_covered, public_splits, stats = check_decoder()
    print(f"True masks recovered all {pairs} pairs in [64,128). "
          f"Zero-template error radii 0,1,2 cover {true_covered} true masks "
          f"and split {public_splits} products.")
    print("Public two-error run on 881 * 967:", stats)
