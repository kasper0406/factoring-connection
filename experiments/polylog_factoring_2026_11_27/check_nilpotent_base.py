"""Cyclotomic evaluation at multiples of three using short divisor scans."""

from itertools import combinations

from sympy import mobius, primerange

from check_cyclotomic_reduction import reference_residue


def nilpotent_residue(n, base, precision, mu):
    """Evaluate Φ_n(base) mod 3^precision with a Möbius-value callback."""
    if n < 1 or precision < 1 or base % 3:
        raise ValueError("Expected positive index and precision and a multiple of 3")
    modulus = 3**precision
    if n == 1:
        return (base - 1) % modulus
    result = 1
    for divisor in range(1, min(n+1, precision)):
        if n % divisor == 0:
            term = 1 - pow(base, divisor, modulus)
            result = result * pow(term, int(mu(n//divisor)), modulus) % modulus
    return result


def semiprime_residue(n, base, precision):
    """Use the promise n=pq for distinct primes, including when n is even."""
    def quotient_mu(quotient):
        if quotient == n or quotient == 1:
            return 1
        return -1

    return nilpotent_residue(n, base, precision, quotient_mu)


if __name__ == "__main__":
    general = mixed = 0
    for n in range(2, 501):
        assert reference_residue(n, 9, 3) == (1-3*mobius(n)) % 9
        for base in (3, 6, 9):
            for precision in (1, 4, 12):
                assert nilpotent_residue(n, base, precision, mobius) == reference_residue(
                    n, 3**precision, base
                )
                general += 1
    for p, q in combinations(primerange(2, 80), 2):
        for precision in (p, p+1):
            value = semiprime_residue(p*q, 3, precision)
            assert value == reference_residue(p*q, 3**precision, 3)
            assert (value == (-2) % (3**precision)) == (precision == p)
            mixed += 1
    print(f"{general} general residues and {mixed} semiprime boundary cases verified")
