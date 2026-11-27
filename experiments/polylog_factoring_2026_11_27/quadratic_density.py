"""Certified bounds for Pr[gcd(x²-d, n)=1] for uniformly random x modulo n."""

from fractions import Fraction
from math import gcd


def _jacobi(a, n):
    """Jacobi symbol for odd positive n, by quadratic reciprocity."""
    a %= n
    sign = 1
    while a:
        twos = (a & -a).bit_length() - 1
        a >>= twos
        if twos % 2 and n % 8 in (3, 5):
            sign = -sign
        if a % 4 == n % 4 == 3:
            sign = -sign
        a, n = n % a, a
    return sign if n == 1 else 0


def _local_density(prime, d):
    if prime == 2 or d % prime == 0:
        roots = 1
    else:
        roots = 2 if pow(d, (prime - 1) // 2, prime) == 1 else 0
    return Fraction(prime - roots, prime)


def density_interval(n, d, cutoff=None):
    """Return rational lower and upper bounds, trial-dividing up to cutoff."""
    if n < 1:
        raise ValueError("n must be positive")
    if cutoff is None:
        cutoff = 2 * n.bit_length() ** 2
    if cutoff < 2:
        raise ValueError("cutoff must be at least 2")

    remaining = n
    known = Fraction(1)
    p = 2
    while p <= cutoff and p * p <= remaining:
        if remaining % p == 0:
            known *= _local_density(p, d)
            while remaining % p == 0:
                remaining //= p
        p = 3 if p == 2 else p + 2

    if remaining == 1:
        return known, known
    if p * p > remaining:
        density = known * _local_density(remaining, d)
        return density, density

    # Every residual prime is at least p. Bound their number by log_p(remaining).
    count_bound = 0
    size = remaining
    while size >= p:
        size //= p
        count_bound += 1
    upper = known
    if _jacobi(d, remaining) == -1:
        # At least one prime with odd multiplicity has no roots.
        if count_bound == 2:
            witness = pow(d, (remaining - 1) // 2, remaining)
            for sign in (-1, 1):
                factor = gcd(witness - sign, remaining)
                if 1 < factor < remaining:
                    # At most two primes counting multiplicity. Negative
                    # Jacobi excludes a square, so both pieces are primes.
                    density = known * _local_density(factor, d)
                    density *= _local_density(remaining // factor, d)
                    return density, density
            if witness != remaining - 1:
                # Euler's criterion proves compositeness and opposite signs.
                upper *= 1 - Fraction(2 * p, remaining)
        count_bound -= 1
    lower = known * Fraction(p - 2, p) ** count_bound
    # The true probability is a multiple of 1/n.
    return Fraction(-((-lower * n) // 1), n), upper


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("d", type=int)
    parser.add_argument("--cutoff", type=int)
    args = parser.parse_args()
    print(*density_interval(args.n, args.d, args.cutoff))
