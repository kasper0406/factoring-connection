"""Check the two local cyclic-fold formulas for fixed prime gaps."""

from fractions import Fraction
from math import gcd, isqrt, lcm


def cyclic_multiply(left, right):
    width = len(left)
    answer = [0] * width
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[(i + j) % width] += a * b
    return answer


def cyclic_power(value, exponent, modulus=None):
    result = [1] + [0] * (len(value) - 1)
    while exponent:
        if exponent & 1:
            result = cyclic_multiply(result, value)
            if modulus is not None:
                result = [x % modulus for x in result]
        exponent //= 2
        if exponent:
            value = cyclic_multiply(value, value)
            if modulus is not None:
                value = [x % modulus for x in value]
    return result


def rational_at_q(width, h):
    if width % 2:
        inverse = [Fraction((-1)**i, 2) for i in range(width)]
        return cyclic_power(inverse, 2 * h)[0]
    projector = [Fraction((-1)**i, width) for i in range(width)]
    inverse = []
    previous = 0
    for i in range(width):
        previous = int(i == 0) - projector[i] - previous
        inverse.append(previous)
    assert inverse[-1] == 0
    power = cyclic_power(inverse, 2 * h)
    at_minus_one = sum((-1)**i * x for i, x in enumerate(power))
    return power[0] - at_minus_one * projector[0]


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def check(widths, limit):
    progression = 2 * lcm(*widths)
    primes = [p for p in range(progression + 1, limit, progression) if is_prime(p)]
    nonunit_cases = 0
    cases = 0
    for p, q in zip(primes, primes[1:]):
        h = (q - p) // 2
        if q >= 2 * p or h * h >= p + q - 1:
            continue
        n = p * q
        b = isqrt(n)
        assert b == p + h - 1
        all_units = True
        for width in widths:
            assert b % width == 0
            base = [0] * width
            base[0] += 1
            base[1 % width] += 1
            actual = cyclic_power(base, 2 * b, n)[0]
            at_p = cyclic_power(base, 2 * h)[0]
            at_q = rational_at_q(width, h)
            assert at_p > 0 and at_q > 0
            assert actual % p == at_p % p
            assert actual % q == at_q.numerator * pow(at_q.denominator, -1, q) % q
            all_units &= gcd(actual, n) == 1
        cases += 1
        nonunit_cases += not all_units
    print(f"{widths}: {cases} pairs, {cases - nonunit_cases} with all folds units")
    return cases


if __name__ == "__main__":
    assert sum(check(widths, 6000) for widths in [(1, 2, 3), (2, 4, 8), (3, 5, 6)]) >= 100
