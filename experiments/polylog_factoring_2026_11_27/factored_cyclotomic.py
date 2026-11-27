"""Evaluate a cyclotomic residue near 1 from a supplied prime factorization."""

from math import comb, prod


def factored_residue(factors, base, precision):
    """Compute Φ_n(base) mod 3^precision for base congruent to 1 mod 3.

    factors maps each distinct prime dividing n to its positive multiplicity.
    """
    if precision < 1 or base % 3 != 1:
        raise ValueError("Expected positive precision and a base congruent to 1 mod 3")
    modulus = 3**precision
    if not factors:
        return (base - 1) % modulus
    # Absorb repeated prime powers into the evaluation point.
    exponent = prod(p**(e-1) for p, e in factors.items())
    point = pow(base, exponent, modulus)
    primes = sorted(factors, key=lambda p: p == 3)
    first = primes[0]
    coefficients = [comb(first, k+1) % modulus for k in range(precision)]
    for prime in primes[1:]:
        # Substitute (1+T)^prime-1 by modular Horner composition.
        shift = [0] + [comb(prime, k) % modulus for k in range(1, precision)]
        numerator = [0] * precision
        for coefficient in reversed(coefficients):
            numerator = [coefficient] + [
                sum(numerator[k-j]*shift[j] for j in range(1, k+1)) % modulus
                for k in range(1, precision)
            ]
        inverse = pow(coefficients[0], -1, modulus)
        quotient = []
        for k in range(precision):
            known = sum(coefficients[j]*quotient[k-j] for j in range(1, k+1))
            quotient.append((numerator[k]-known)*inverse % modulus)
        coefficients = quotient
    result = 0
    for coefficient in reversed(coefficients):
        result = (result*(point-1) + coefficient) % modulus
    return result


if __name__ == "__main__":
    from sympy import factorint
    from check_cyclotomic_reduction import reference_residue

    cases = 0
    for n in range(1, 501):
        factors = factorint(n)
        for precision in (1, 4, 9):
            modulus = 3**precision
            for base in (4, 16, pow(4, -1, modulus)):
                # The exact divisor recursion requires base > 1.
                representative = base if base > 1 else base + modulus
                assert factored_residue(factors, base, precision) == reference_residue(
                    n, modulus, representative
                )
                cases += 1
    print(f"{cases} factored residue evaluations verified")
