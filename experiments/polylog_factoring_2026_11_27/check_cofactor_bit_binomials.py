"""Recover all prime powers from the supports of a fixed dyadic binomial menu."""

from itertools import combinations
from math import comb, factorial, gcd

from sympy import factorint, integer_nthroot, primerange


def direct_ratio(top, modulus):
    half = top//2
    coefficient = comb(2*half, half)
    if top % 2:
        coefficient *= top
    return coefficient % modulus


def factor_from_dyadic_coefficients(n, evaluate=direct_ratio):
    """Split an odd n when evaluate(top,n) supplies floor(top)! / floor(top/2)!² mod n."""
    if n <= 1 or n % 2 == 0:
        raise ValueError("expected an odd integer above one")
    parts = [n]
    for j in range(1, n.bit_length()):
        top = n >> j
        if top < 2:
            break
        residue = evaluate(top, n)
        support = gcd(pow(residue, n.bit_length(), n), n)
        if support == 1 or support == n:
            continue
        refined = []
        for part in parts:
            left = gcd(part, support)
            if 1 < left < part:
                refined.extend((left, part//left))
            else:
                refined.append(part)
        parts = refined
    factors = {}
    for part in parts:
        for exp in range(part.bit_length()-1, 1, -1):
            base, exact = integer_nthroot(part, exp)
            if exact:
                factors[base] = exp
                break
        else:
            factors[part] = 1
    return factors


def valuation(value, p):
    power = 0
    while value % p == 0:
        value //= p
        power += 1
    return power


def evaluate_integer(top):
    half = top//2
    return (top if top % 2 else 1)*comb(2*half, half)


def check_balanced():
    cases = 0
    for p, q in combinations(primerange(3, 80), 2):
        if q >= 2*p:
            continue
        n = p*q
        for c in (2, 3, 4, 5, 7, 13, 16, 17):
            ratio = evaluate_integer(n//c)
            assert valuation(ratio, p) == (q//c) % 2
            assert valuation(ratio, q) == (p//c) % 2
            cases += 1
    return cases


def check_general():
    checks = support_pairs = 0
    for n in [*range(3, 1200, 2), 3*5*7*11, 3**5*5*7, 7**3*11, 3*4001, 1009*1013]:
        factors = factorint(n)
        assert factor_from_dyadic_coefficients(n) == factors
        signatures = {}
        for p in factors:
            quotient, bits = n//p, 0
            while quotient:
                bits |= quotient
                quotient //= p
            signatures[p] = bits >> 1
        for p, q in combinations(sorted(factors), 2):
            j = ((n//p) ^ (n//q)).bit_length()-1
            assert j > 0 and (1 << j)*q*q > n
            first = signatures[p] >> (j-1)
            second = signatures[q] >> (j-1)
            assert first % 2 == 1 and second % 2 == 0 and first > second
            support_pairs += 1
        if n < 900:
            product = 1
            for j in range(1, n.bit_length()):
                product *= evaluate_integer(n >> j)**(1 << (j-1))
            assert product == factorial(n//2)
        checks += 1
    return checks, support_pairs


if __name__ == "__main__":
    print("Exact balanced cutoff checks", check_balanced())
    print("All-input reconstructions and ordered prime pairs", check_general())
