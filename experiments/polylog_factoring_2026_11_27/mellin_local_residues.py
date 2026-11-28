"""Divisor-sum residues forced by the product of two unknown unit residues."""

from math import gcd, isqrt


def local_h_residue(n):
    """Return (H(n) mod m, m) using n alone on the promise n=pq, p,q>3 distinct.

    H(pq) = (p²+p+1)(q²+q+1). The modulus is maximal for a rule based on
    multiplication of local unit residues. The semiprime promise is supplied
    by the caller. Inputs divisible by 2 or 3 already have a known divisor.
    """
    if n <= 1 or gcd(n, 6) != 1:
        raise ValueError("expected a semiprime with distinct prime factors greater than three")
    if n % 4 == 1:
        modulus = 8
    elif n % 8 == 3:
        modulus = 32
    elif n % 16 == 15:
        modulus = 64
    else:  # n = 7 mod 16
        modulus = 128
    if n % 3 == 2:
        modulus *= 9
    x = n % modulus
    return 3 * (x*x + x + 1) % modulus, modulus


def _local_factor(x, j, modulus):
    # After removing x^(j-1) from H(x^j), the remaining polynomial is
    # x^(j+1) + (x+1)(1+x+...+x^(j-1)). Division by x-1 need not exist.
    geometric_sum = 0
    for _ in range(j):
        geometric_sum = (x * geometric_sum + 1) % modulus
    return (pow(x, j + 1, modulus) + (x + 1) * geometric_sum) % modulus


def forced_power_residue(n, j):
    """Return (residue, largest coprime modulus) forced for H((pq)^j).

    The promise is n=pq with distinct primes. At each modulus the only
    information used about the factors is that they are units of product n.
    Work is polynomial in the numerical exponent j and in the bit length of n.
    """
    if not isinstance(n, int) or n < 1 or not isinstance(j, int) or j < 1:
        raise ValueError("expected positive integers n and j")
    degree = 2 * j + 2
    modulus = 1
    for prime in range(2, degree + 2):
        if n % prime == 0 or any(prime % k == 0 for k in range(2, isqrt(prime) + 1)):
            continue

        # Unit residues recur in blocks of phi(prime^i) consecutive units.
        phi = prime - 1
        valuation = 0
        while phi <= degree:
            valuation += degree // phi
            phi *= prime
        cap = prime**valuation
        reference = (2*j + 1) * _local_factor(n % cap, j, cap)
        remaining = cap
        x = count = 0
        while count <= degree and remaining > 1:
            x += 1
            if x % prime == 0:
                continue
            count += 1
            y = n * pow(x, -1, remaining) % remaining
            value = _local_factor(x, j, remaining) * _local_factor(y, j, remaining)
            remaining = gcd(remaining, value - reference)
        modulus *= remaining

    residue = pow(n, j - 1, modulus) * (2*j + 1) * _local_factor(n % modulus, j, modulus)
    return residue % modulus, modulus
