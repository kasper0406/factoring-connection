"""Try gcds of successive digits of (1+x)^n - x^n - 1."""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm


def defect(n, x, modulus):
    return (pow(x + 1, n, modulus) - pow(x, n, modulus) - 1) % modulus


def split_additive_power(n, x, precision):
    """Return (factor, digit), or None. Digit zero denotes a denominator factor."""
    if n < 3 or n % 2 == 0 or precision < 1:
        raise ValueError("expected odd n >= 3 and positive precision")
    x = Fraction(x)
    divisor = gcd(x.denominator, n)
    if divisor == n:
        raise ValueError("denominator is zero modulo n")
    if divisor != 1:
        return divisor, 0
    modulus = n**precision
    base = x.numerator * pow(x.denominator, -1, modulus) % modulus
    value = defect(n, base, modulus)
    for digit in range(1, precision + 1):
        divisor = gcd(value, n)
        if 1 < divisor < n:
            return divisor, digit
        if divisor == 1 or value == 0:
            break
        value //= n
    return None


def valuation(value, prime, cap):
    for e in range(cap):
        if value % prime:
            return e
        value //= prime
    return cap


def prime_control(n, x, p, precision):
    modulus = p**precision
    base = x.numerator * pow(x.denominator, -1, modulus) % modulus
    value = defect(n, base, modulus)
    return valuation(value, p, precision)


def predicted_split(n, x, p, q, precision):
    vp = prime_control(n, x, p, precision)
    vq = prime_control(n, x, q, precision)
    if vp == vq:
        return None
    return (p, vq + 1) if vp > vq else (q, vp + 1)


def endpoint_or_cubic_split(n, x):
    for value in (x, x + 1, x * x + x + 1):
        while value and value % n == 0:
            value //= n
        if 1 < gcd(value, n) < n:
            return True
    return False


def check_exponent_differences():
    from sympy import n_order, primerange

    primes = list(primerange(2, 54))
    counts = Counter()
    for r in primes:
        for x in range(1, r - 1):
            period = lcm(n_order(x, r), n_order(x + 1, r))
            hits = {
                e for e in range(period)
                if (pow(x + 1, e, r) - pow(x, e, r)) % r == 1
            }
            differences = [(a - b) % period for a in hits for b in hits if a != b]
            assert 0 not in differences
            assert len(differences) == len(set(differences))
            counts["exponent_sets"] += 1
            counts["exponents_tested"] += period
            counts["ordered_differences"] += len(differences)
            for s in primes:
                if s != r and s != 2:
                    assert (s % period in hits) == (defect(r * s, x, r) == 0)
                    counts["prime_partner_checks"] += 1
    return counts


def check():
    from sympy import primerange

    counts = check_exponent_differences()
    precision = 5
    for p, q in combinations(primerange(3, 38), 2):
        if q >= 2 * p:
            continue
        n = p * q
        counts["balanced_moduli"] += 1
        for r, s in ((p, q), (q, p)):
            inverses = [pow(j, -1, r) for j in range(1, r)]
            for x in range(r):
                value = defect(n, x, r * r)
                if value % r:
                    continue
                b = pow(x, s, r)
                logarithm = -sum(inv * pow(-b, j, r)
                                 for j, inv in enumerate(inverses, 1))
                assert value // r == logarithm % r
                counts["finite_log"] += 1
        modulus = n**precision
        for x in range(n):
            result = split_additive_power(n, x, precision)
            assert result == predicted_split(n, x, p, q, precision)
            counts[f"digit_{result[1]}" if result else "no_split"] += 1
            if result and endpoint_or_cubic_split(n, x):
                counts[f"digit_{result[1]}_with_endpoint_or_cubic_split"] += 1
            value = defect(n, x, modulus)
            translated = defect(n, x + n, modulus)
            slope = (pow(x + 1, n - 1, n) - pow(x, n - 1, n)) % n
            assert (translated - value) % (n**3) == n * n * slope
            after_translation = split_additive_power(n, x + n, precision)
            if result != after_translation:
                assert all(hit is None or hit[1] >= 3 for hit in (result, after_translation))
                counts["changed_on_translation"] += 1
                if endpoint_or_cubic_split(n, x) or endpoint_or_cubic_split(n, x + n):
                    counts["translation_with_endpoint_or_cubic_split"] += 1
            for r in (p, q):
                cubic = x * x + x + 1
                if p == 3 or cubic % r:
                    continue
                multiplicity = 2 if n % 3 == 1 else 1
                expected = 1 + multiplicity * valuation(cubic, r, precision)
                assert valuation(value, r, precision) == min(precision, expected)
                counts["cubic_valuation"] += 1
        for x in (Fraction(-1, 2), Fraction(2, 3), Fraction(5, 7), Fraction(n + 2, 2)):
            result = split_additive_power(n, x, precision)
            if gcd(n, x.denominator) > 1:
                assert result == (gcd(n, x.denominator), 0)
            else:
                assert result == predicted_split(n, x, p, q, precision)
                a, b = x.numerator, x.denominator
                homogeneous = pow(a + b, n, modulus) - pow(a, n, modulus) - pow(b, n, modulus)
                base = a * pow(b, -1, modulus)
                affine = defect(n, base, modulus)
                assert homogeneous % modulus == affine * pow(b, n, modulus) % modulus
            counts["rational_controls"] += 1
    # Separately selected witnesses, excluded from the exhaustive cohort counts.
    fixtures = ((79, 157, 11, 79), (2137, 3881, 113370, 3881),
                (2003, 5333, Fraction(439, 723), 2003))
    for p, q, x, expected in fixtures:
        n = p * q
        residue = x.numerator * pow(x.denominator, -1, n) % n
        for divisor in (p, q):
            u = residue % divisor
            orbit_tests = (u, u + 1, 2 * u + 1, u - 1, u + 2, u * u + u + 1)
            assert all(test % divisor for test in orbit_tests)
        assert split_additive_power(n, x, 1) is None
        assert split_additive_power(n, x, 2) == (expected, 2)
    return dict(counts)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("n", nargs="?", type=int)
    parser.add_argument("x", nargs="?", type=Fraction)
    parser.add_argument("precision", nargs="?", type=int, default=3)
    args = parser.parse_args()
    if args.n is None:
        print(check())
    elif args.x is None:
        parser.error("provide both n and x")
    else:
        print(split_additive_power(args.n, args.x, args.precision))
