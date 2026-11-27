"""Decode p+q from two reference values at index 3pq."""

from itertools import combinations
from math import factorial
from sympy import primerange

from check_cyclotomic_reduction import reference_residue


def norm_quotient(n, total, power, precision):
    """Compute the norm of G(power*T) in T²-total*T+n=0."""
    modulus = 3**precision

    def multiply(left, right):
        a, b = left
        c, d = right
        return ((a*c - n*b*d) % modulus, (a*d + b*c + total*b*d) % modulus)

    # Expand (4^(power*T)-1)/3 via the falling factorial.
    falling = (1, 0)
    constant = linear = 0
    for k in range(1, 2*precision + 2):
        falling = multiply(falling, (1-k, power))
        denominator = factorial(k)
        coefficient = 3**(k-1)
        while denominator % 3 == 0:
            denominator //= 3
            coefficient //= 3
        coefficient *= pow(denominator, -1, modulus)
        constant += coefficient * falling[0]
        linear += coefficient * falling[1]
    square = multiply((constant, linear), (constant, linear))
    a = (1 + 3*(constant + square[0])) % modulus
    b = (3*(linear + square[1])) % modulus
    return (a*a + total*a*b + n*b*b) % modulus


def decode(n, residues, digits):
    modulus = 3**(digits + 1)
    targets = []
    for power, residue in zip((1, 2), residues, strict=True):
        def g(u):
            y = pow(4, power*u, 3*modulus)
            return (y*y + y + 1)//3
        targets.append(g(n)*g(1)*pow(residue, -1, modulus) % modulus)
    total, period = 0, 1
    for precision in range(2, digits + 2):
        matches = []
        for candidate in range(total, 3*period, period):
            difference = (norm_quotient(n, candidate, 1, precision)
                          - norm_quotient(n, candidate, 2, precision))
            if (difference - targets[0] + targets[1]) % (3**precision) == 0:
                matches.append(candidate)
        if len(matches) != 1:
            raise ArithmeticError("Expected a unique sum digit")
        total, = matches
        period *= 3
    return total


if __name__ == "__main__":
    checked = 0
    for p, q in combinations(primerange(5, 80), 2):
        n = p*q
        digits = 1
        while 3**digits <= n:
            digits += 1
        values = [reference_residue(3*n, 3**(digits+1), base) for base in (4, 16)]
        assert decode(n, values, digits) == p+q
        checked += 1
    print(f"{checked} distinct-prime sums verified")
