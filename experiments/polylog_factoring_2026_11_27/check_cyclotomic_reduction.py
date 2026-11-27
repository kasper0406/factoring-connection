"""End-to-end control for factoring from one fixed-base cyclotomic residue."""

from collections import Counter
from functools import cache
from math import gcd
from random import Random

from sympy import divisors, factorint, isprime, perfect_power


def factor_with_query(n, residue, rng):
    """Factor a positive integer using a correct Φ_index(4) mod modulus oracle."""
    if n < 1:
        raise ValueError("expected a positive integer")
    factors = Counter()
    for p in (2, 3):
        while n % p == 0:
            factors[p] += 1
            n //= p
    if n == 1:
        return factors

    scale = 3
    while scale <= n:
        scale *= 3
    modulus = scale**2
    value = residue(scale*n, modulus) % modulus
    y = pow(4, scale//3, modulus)
    delta = (y**3 - y) % modulus
    if value % scale != 1:
        raise ValueError("residue is inconsistent with cyclotomic linearization")
    phi = (2*((value-1)//scale)*pow(delta//scale, -1, scale)) % scale
    if not 0 < phi < n or phi % 2:
        raise ValueError("invalid decoded exponent")
    odd_part = phi
    while odd_part % 2 == 0:
        odd_part //= 2
    squarings = (phi//odd_part).bit_length()-1

    pending = [(n, 1)]
    while pending:
        part, multiplicity = pending.pop()
        if isprime(part):
            factors[part] += multiplicity
            continue
        power = perfect_power(part, factor=False)
        if power:
            base, exponent = map(int, power)
            pending.append((base, exponent*multiplicity))
            continue
        while True:
            base = rng.randrange(1, part)
            divisor = gcd(base, part)
            if divisor == 1:
                term = pow(base, odd_part, part)
                for _ in range(squarings+1):
                    divisor = gcd(term-1, part)
                    if divisor != 1:
                        break
                    term = term**2 % part
            if 1 < divisor < part:
                break
        pending.extend(((divisor, multiplicity), (part//divisor, multiplicity)))
    return factors


def reference_residue(index, modulus, base=4):
    """Slow exact evaluation via base^k−1 = product of Φ_d(base) over d dividing k."""
    @cache
    def value(k):
        result = base**k-1
        for d in divisors(k)[:-1]:
            result //= value(d)
        return result
    return value(index) % modulus


if __name__ == "__main__":
    rng = Random(994)
    inputs = [*range(1, 301), 385, 539, 625, 770, 1155]
    queries = 0
    for n in inputs:
        calls = []

        def query(index, modulus):
            calls.append(index)
            return reference_residue(index, modulus)

        assert factor_with_query(n, query, rng) == factorint(n)
        cofactor = n
        for p in (2, 3):
            while cofactor % p == 0:
                cofactor //= p
        assert len(calls) == int(cofactor > 1)
        queries += len(calls)
    print(f"{len(inputs)} factorizations checked, {queries} residue queries")
