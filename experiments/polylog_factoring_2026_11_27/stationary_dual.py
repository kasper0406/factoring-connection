"""Certified nonzero Fourier coefficients of log magnitude on the full grid.

For odd M, extend the stationary trace to all t modulo M by
F(t) = sum(e_M(2*t*u) for u^2 = 1 mod M). For unit t and M=N² this
is the normalized square-modulus Kloosterman sum. This tool evaluates the
Fourier coefficients of log|F| at nonzero frequencies from M and k alone.
The mean is outside this evaluator.

SymPy supplies exact integer roots and the modular arithmetic for a
deterministic seven-base Miller–Rabin test below 2**64. Above that range,
possible prime bases are proved with python-flint.
The polynomial bit-time result uses deterministic polynomial-time primality.
"""

from argparse import ArgumentParser
from dataclasses import dataclass
from decimal import ROUND_CEILING, Decimal, localcontext
from fractions import Fraction
from math import gcd, isqrt, prod

from sympy import integer_nthroot, isprime
from sympy.ntheory.primetest import mr


def proven_prime(value):
    if value < 2**64:
        if value == 2:
            return True
        return value > 2 and value % 2 == 1 and mr(
            value, [2, 325, 9375, 28178, 450775, 9780504, 1795265022],
        )
    if not isprime(value):
        return False
    try:
        from flint import fmpz
    except ImportError as error:
        raise RuntimeError('Install python-flint to prove primality of larger candidates') from error
    return bool(fmpz(value).is_prime())


def prime_power_base(value):
    """Return p if value=p^e for prime p and e>=1, or None."""
    for exponent in range(value.bit_length() - 1, 1, -1):
        base, exact = integer_nthroot(value, exponent)
        if exact:
            value = base  # This is the maximal exact exponent.
            break
    return value if proven_prime(value) else None


@dataclass(frozen=True)
class Coefficient:
    conductor: int
    local_modulus: int | None
    alias: int | None
    sign: int
    estimate: Fraction
    error: Fraction


def coefficient(modulus, frequency, terms=40):
    """Return an exact rational estimate, one-sided error bound and exact sign."""
    if modulus <= 1 or modulus % 2 == 0 or terms < 0:
        raise ValueError('Expected an odd modulus greater than one and nonnegative term count')
    frequency %= modulus
    if frequency == 0:
        raise ValueError('Use a nonzero frequency modulo the modulus')
    conductor = modulus // gcd(frequency, modulus)
    prime = prime_power_base(conductor)
    if prime is None:
        return Coefficient(conductor, None, None, 0, Fraction(0), Fraction(0))

    q, complement = 1, modulus
    while complement % prime == 0:
        q *= prime
        complement //= prime
    # The first local harmonic is 4/complement mod q. Divide by that
    # unit to translate the requested local frequency into an alias index.
    r = (frequency // complement) * complement * pow(4, -1, q) % q
    s = min(r, q - r)
    sign = -1 if s % 2 == 0 else 1

    near = Fraction(1, s)
    far = Fraction(1, q - s)
    error = (near - far) / (1 << (terms + 1))
    estimate = Fraction(0)
    for index in range(terms):
        estimate += (near - far) / (1 << (index + 2))
        step = (index + 1) * q
        near *= Fraction(step, s + step)
        far *= Fraction(step, q - s + step)
    return Coefficient(conductor, q, s, sign, sign * estimate, error)


def check():
    import numpy as np

    comparisons = 0
    for factors in (((3, 1),), ((3, 2),), ((3, 1), (5, 1)), ((3, 2), (5, 1)),
                    ((3, 1), (5, 1), (7, 1)), ((5, 2), (7, 1))):
        modulus = prod(prime**(2 * exponent) for prime, exponent in factors)
        n = isqrt(modulus)
        grid = np.arange(modulus)
        lower_grid = np.arange(n)
        signal = np.zeros(modulus)
        lift_average = np.zeros(n)
        support = np.zeros(modulus, dtype=bool)
        mean = 0.0
        for prime, exponent in factors:
            q = prime**(2 * exponent)
            complement = modulus // q
            a = pow(complement, -1, q)
            angles = 2 * np.pi * ((2 * a * grid) % q) / q
            signal += np.log(np.abs(2 * np.cos(angles)))
            p_power = isqrt(q)
            lower_angles = 2 * np.pi * ((2 * a * lower_grid) % p_power) / p_power
            lift_average += np.log(np.abs(2 * np.cos(lower_angles))) / p_power
            support[complement::complement] = True
            mean += np.log(2) / q
        assert np.max(abs(signal.reshape(n, n).mean(axis=0) - lift_average)) < 2e-12
        spectrum = np.fft.fft(signal) / modulus
        assert abs(spectrum[0] - mean) < 2e-12
        assert np.array_equal(abs(spectrum[1:]) > 2e-12, support[1:])
        for k in range(1, modulus):
            precise = coefficient(modulus, k, 42)
            assert bool(precise.sign) == support[k]
            assert abs(float(precise.estimate) - spectrum[k]) <= float(precise.error) + 2e-12
            if precise.sign:
                assert precise.sign * spectrum[k].real > 0
                coarse = coefficient(modulus, -k, 3)
                delta = precise.sign * (precise.estimate - coarse.estimate)
                assert 0 < delta <= coarse.error
            comparisons += 1
    print(f'{comparisons} nonzero-frequency coefficients checked against full-grid DFTs')


def display_bounds(result, digits):
    with localcontext() as context:
        context.prec = digits
        estimate = Decimal(result.estimate.numerator) / Decimal(result.estimate.denominator)
        error = result.error + abs(Fraction(estimate) - result.estimate)
        context.rounding = ROUND_CEILING
        displayed_error = Decimal(error.numerator) / Decimal(error.denominator)
        return str(estimate), str(displayed_error)


def exact_string(value):
    # Decimal converts arbitrary-size integers without the int-to-string limit.
    numerator = str(Decimal(value.numerator))
    return numerator if value.denominator == 1 else f'{numerator}/{Decimal(value.denominator)}'


def main():
    parser = ArgumentParser(description=__doc__)
    moduli = parser.add_mutually_exclusive_group()
    moduli.add_argument('--n', type=int, help='Set M=N²')
    moduli.add_argument('--modulus', type=int, help='Set M directly')
    parser.add_argument('frequency', nargs='*', type=int, metavar='K')
    parser.add_argument('--terms', type=int, default=40, help='Number of exact rational terms')
    parser.add_argument('--digits', type=int, default=18, help='Displayed decimal precision')
    parser.add_argument('--exact', action='store_true', help='Display exact rational values')
    args = parser.parse_args()
    if args.n is None and args.modulus is None and not args.frequency:
        check()
        return
    if (args.n is not None and args.n <= 1) or args.digits < 1:
        parser.error('Expected N greater than one and positive display precision')
    modulus = args.n**2 if args.n is not None else args.modulus
    if modulus is None or not args.frequency:
        parser.error('Supply --n or --modulus and at least one nonzero frequency')
    try:
        for k in args.frequency:
            result = coefficient(modulus, k, args.terms)
            if args.exact:
                estimate, error = exact_string(result.estimate), exact_string(result.error)
            else:
                estimate, error = display_bounds(result, args.digits)
            print(f'k={k % modulus}, conductor={result.conductor}, '
                  f'local modulus={result.local_modulus}, sign={result.sign:+d}, '
                  f'estimate={estimate}, absolute error <= {error}')
    except (RuntimeError, ValueError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
