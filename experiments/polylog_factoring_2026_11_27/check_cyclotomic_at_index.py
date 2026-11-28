"""Check cyclotomic values and first derivatives modulo their own index."""

from collections import Counter
from itertools import combinations
from math import gcd

from sympy import Poly, Symbol, cyclotomic_poly, factorint, isprime, n_order, perfect_power, primerange


def value_or_factor(n, a):
    """Return ('value', residue) or ('factor', proper divisor), for n > 1."""
    if n < 2:
        raise ValueError("the index must exceed one")
    divisor = gcd(pow(a, n, n) - 1, n)
    if divisor == 1:
        return "value", 1
    if divisor == n:
        divisor = gcd(a - 1, n)
        if divisor == n:
            power = perfect_power(n, factor=False)
            base = int(power[0]) if power else n
            return "value", base % n if isprime(base) else 1
    assert 1 < divisor < n
    return "factor", divisor


def evaluate_factored_semiprime(p, q, a):
    """Distinct primes p < q. Only q can have the exceptional order."""
    return q if a % q != 1 and pow(a, p, q) == 1 else 1


def generic_derivative(p, q):
    n = p * q
    denominator = (pow(2, n, n) - 1) % n
    coefficient = pow(2, n-1, n) * pow(denominator, -1, n) % n
    return (1 - (p+q) * coefficient) % n


def geometric_sum(base, length, modulus):
    power, total = 1, 0
    for bit in bin(length)[2:]:
        total = (1 + power) * total % modulus
        power = power * power % modulus
        if bit == "1":
            total = (total + power) % modulus
            power = power * base % modulus
    return total


def factored_derivative(p, q):
    n = p*q
    modulus = n*n
    values = []
    for base in (2, 2+n):
        numerator = geometric_sum(pow(base, p, modulus), q, modulus)
        denominator = geometric_sum(base, q, modulus)
        values.append(numerator * pow(denominator, -1, modulus) % modulus)
    difference = (values[1] - values[0]) % modulus
    assert difference % n == 0
    return difference // n


def run_checks():
    x = Symbol("x")
    counts = Counter()
    for n in range(2, 241):
        factors = factorint(n)
        polynomial = Poly(cyclotomic_poly(n, x), x)
        for a in range(-4, 13):
            actual = int(polynomial.eval(a)) % n
            expected = 1
            for p, e in factors.items():
                power = p**e
                exceptional = a % p != 0 and n_order(a % p, p) == n // power
                local = p if exceptional else 1
                assert actual % power == local % power
                if exceptional:
                    expected = p % n
                counts["local"] += 1
            assert actual == expected
            kind, result = value_or_factor(n, a)
            if kind == "value":
                assert result == actual
            else:
                assert 1 < result < n and n % result == 0
            if len(factors) == 2 and all(e == 1 for e in factors.values()):
                p, q = sorted(factors)
                assert evaluate_factored_semiprime(p, q, a) == actual
            counts[kind] += 1
    for p, q in combinations(primerange(3, 50), 2):
        n = p*q
        actual = Poly(cyclotomic_poly(n, x), x).diff().eval(2) % n
        assert factored_derivative(p, q) == actual
        if gcd(pow(2, n, n) - 1, n) > 1:
            assert value_or_factor(n, 2)[0] == "factor"
            counts["base_two_split"] += 1
            continue
        assert actual == generic_derivative(p, q)
        counts["derivative"] += 1
    return dict(counts)


if __name__ == "__main__":
    print(run_checks())
