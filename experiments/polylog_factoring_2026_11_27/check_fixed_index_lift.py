"""Check the fixed-index ternary decoder against exact divisor recursion."""

from sympy import primerange

from check_cyclotomic_reduction import reference_residue


def unit_quotient(exponent, modulus):
    """Return (4^(2u) + 4^u + 1)/3 mod modulus with exact division."""
    power = pow(4, exponent, 3 * modulus)
    return ((power * power + power + 1) // 3) % modulus


def factor_digits(n, residue, digits):
    """Yield q mod 3^k for k=1,...,digits on the mixed semiprime promise."""
    if digits < 1 or n % 3 != 2 or residue % 9 != 1:
        raise ValueError("Expected positive precision and mixed-semiprime residues")
    factor = 2
    period = 3
    yield factor
    for _ in range(1, digits):
        next_period = 3 * period
        modulus = 3 * next_period
        target = unit_quotient(n, modulus) * unit_quotient(1, modulus) % modulus
        matches = []
        for candidate in range(factor, next_period, period):
            other = n * pow(candidate, -1, next_period) % next_period
            product = unit_quotient(candidate, modulus) * unit_quotient(other, modulus)
            if residue * product % modulus == target:
                matches.append(candidate)
        if len(matches) != 1:
            raise ArithmeticError("Expected exactly one ternary lift")
        factor, = matches
        period = next_period
        yield factor


if __name__ == "__main__":
    pairs = levels = 0
    primes = list(primerange(5, 80))
    for p in primes:
        for q in primes:
            if p % 3 != 1 or q % 3 != 2:
                continue
            n = p * q
            digits = 1
            while 3**digits <= n:
                digits += 1
            residue = reference_residue(3 * n, 3 ** (digits + 1))
            for k, result in enumerate(factor_digits(n, residue, digits), 1):
                assert result == q % (3**k)
                levels += 1
            assert result == q
            pairs += 1
    print(f"{pairs} mixed semiprimes, {levels} ternary levels verified")
