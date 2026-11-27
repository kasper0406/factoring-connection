"""Native zero-residue bits for the complete idempotent count, using lattices.

For prime-power components q of N, the count at A is product(q+A). A
component with residue r=-A mod ell either appears as a unitary divisor
d <= sqrt(N) in residue r, or has a unitary complement d <= sqrt(N) in
residue N/r. Check the appropriate side for a prime power at each hit.
The input uses a supplied prime ell and no factorization of N.

In a band [2**j, 2**(j+1)), parameterize one candidate progression as a
monic linear polynomial f(x) modulo N, with unknown divisor d >= 2**j.
Use shifts N**(k-i)*f**i for 0<=i<k and x**i*f**k for 0<=i<D-k.
After scaling x by block width X, the determinant is
N**(k*(k+1)/2)*X**(D*(D-1)/2). Its LLL L1 bound is below (2**j)**k
for D=n=N.bit_length(), k=j, and X<=2**(j*j/n)/16. The executable
chooses a smaller certified lattice when possible and checks the reduced
row explicitly. There are soft-O(1+N**(1/4)/ell) total blocks since
j-j*j/n <= n/4. This is polynomial bit time for ell >= N**(1/4)
with polynomial-time primality and LLL. It gives one bit per A, and
does not speed up fixed small auxiliary primes. The CLI uses synthetic data.
"""

import argparse
import json
import sys
from collections import Counter
from math import comb, gcd, isqrt, prod
from pathlib import Path
from time import perf_counter

from flint import fmpz_mat
from sympy import Matrix, Poly, Rational, Symbol, nextprime

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.stationary_dual import prime_power_base


def certified_parameters(modulus, band, width, work):
    """Find a dimension and multiplicity with a provable scaled L1 cutoff."""
    for dimension in range(2, modulus.bit_length() + 1):
        common = (dimension**(2 * dimension) << (dimension * (dimension - 1)))
        common *= width**(2 * dimension * (dimension - 1))
        for multiplicity in range(1, dimension):
            work['parameter_checks'] += 1
            # Raise the L1 bound to 4*dimension to compare exact integers.
            left = common * modulus**(2 * multiplicity * (multiplicity + 1))
            right = 1 << (4 * dimension * multiplicity * band)
            if left < right:
                return multiplicity, dimension
    raise ArithmeticError('The selected block exceeds the certified lattice width')


def bounded_roots(modulus, band, offset, width, table, work):
    """Cover every bounded root modulo any divisor at least 2**band."""
    if width <= 1:
        return range(width + 1)
    key = band, width
    if key not in table:
        table[key] = certified_parameters(modulus, band, width, work)
    multiplicity, dimension = table[key]
    scale = [width**index for index in range(dimension)]
    rows = []
    for power in range(multiplicity):
        row = [modulus**(multiplicity - power) * comb(power, index)
               * offset**(power - index) * scale[index] for index in range(power + 1)]
        rows.append(row + [0] * (dimension - len(row)))
    for power in range(dimension - multiplicity):
        row = [0] * power + [comb(multiplicity, index) * offset**(multiplicity - index)
                             * scale[power + index] for index in range(multiplicity + 1)]
        rows.append(row + [0] * (dimension - len(row)))

    work['lattice_calls'] += 1
    work['largest_rank'] = max(work['largest_rank'], dimension)
    reduced = fmpz_mat(rows).lll()
    vector = [int(reduced[0, index]) for index in range(dimension)]
    cutoff = 1 << (band * multiplicity)
    if sum(abs(value) for value in vector) >= cutoff:
        work['exact_fallbacks'] += 1
        vector = list(Matrix(rows).lll(delta=Rational(3, 4))[0, :])
    if (not any(vector) or sum(abs(value) for value in vector) >= cutoff
            or any(value % weight for value, weight in zip(vector, scale))):
        raise ArithmeticError('Reduced row did not certify all short roots')
    coefficients = [value // weight for value, weight in zip(vector, scale)]
    polynomial = Poly.from_list(coefficients[::-1], Symbol('x'))
    return [int(root) for root in polynomial.ground_roots() if root.is_Integer and 0 <= root <= width]


def zero_residue_bit(modulus, prime, selector_count):
    """Return the complete idempotent count's zero bit modulo a supplied prime."""
    if modulus < 1 or prime < 2:
        raise ValueError('Expected a positive input modulus and a prime answer modulus')
    work = Counter()
    target = -selector_count % prime
    residual = modulus
    while residual % prime == 0:
        residual //= prime
        work['removed_prime_exponent'] += 1
    if target == 0:
        return residual != modulus, work
    if residual == 1:
        return False, work
    if residual % prime == target:
        work['prime_power_tests'] += 1
        if prime_power_base(residual) is not None:
            return True, work

    limit = isqrt(residual)
    inverse_step = pow(prime, -1, residual)
    companion = residual * pow(target, -1, prime) % prime
    table = {}
    bit_length = residual.bit_length()
    for residue, complement in ((target, False), (companion, True)):
        for band in range(1, limit.bit_length()):
            low = 1 << band
            high = min(limit, 2 * low - 1)
            first = low + (residue - low) % prime
            if first > high:
                work['empty_bands'] += 1
                continue
            last_index = (high - first) // prime
            block_limit = 1 << max(0, band * band // bit_length - 4)
            work['nonempty_bands'] += 1
            for start in range(0, last_index + 1, block_limit + 1):
                width = min(block_limit, last_index - start)
                offset = (first * inverse_step + start) % residual
                roots = bounded_roots(residual, band, offset, width, table, work)
                work['blocks'] += 1
                work['represented_terms'] += width + 1
                for root in roots:
                    work['root_tests'] += 1
                    candidate = first + prime * (start + root)
                    if residual % candidate:
                        continue
                    cofactor = residual // candidate
                    if gcd(candidate, cofactor) != 1:
                        continue
                    work['unitary_candidates'] += 1
                    work['prime_power_tests'] += 1
                    if prime_power_base(cofactor if complement else candidate) is not None:
                        work['complement_hit' if complement else 'component_hit'] += 1
                        return True, work
    return False, work


def check():
    fixtures = ((), (7,), (27,), (3, 5, 7), (3, 7, 125), (9, 25, 7, 11), (8, 9, 11),
                (101, 109), (101, 1009), (101**2, 109**3, 127),
                (int(nextprime(2**22 + 81)), int(nextprime(2**29 + 818))))
    checked = 0
    for components in fixtures:
        modulus = prod(components)
        for prime in (2, 3, 7, 11, 31, 137, 1009, 32771):
            if modulus > 10**9 and prime < 1000:
                continue
            for selector_count in range(min(7, prime)):
                observed, _ = zero_residue_bit(modulus, prime, selector_count)
                assert observed == (prod(q + selector_count for q in components) % prime == 0)
                checked += 1

    # 125 is too large to appear in the direct bands. Its complement is 21.
    positive, work = zero_residue_bit(3 * 7 * 125, 31, -125)
    assert positive and work['complement_hit'] == 1
    # Complement 15 here has the target residue, but is not a component.
    negative, work = zero_residue_bit(3 * 5 * 7, 11, -15)
    assert not negative and work['unitary_candidates'] > 0
    print(f'{checked} synthetic bits matched, including repeated powers and prime auxiliary divisors')


def benchmark():
    components = (int(nextprime(2**17 + 341)), int(nextprime(2**48 + 251)))
    modulus = prod(components)
    for prime in [int(nextprime(2**exponent)) for exponent in (21, 18, 15)] + [1009]:
        positive = -components[-1] % prime
        negative = next(a for a in range(1, prime) if prod(q+a for q in components) % prime)
        for selector_count in (positive, negative):
            started = perf_counter()
            zero, work = zero_residue_bit(modulus, prime, selector_count)
            assert zero == (prod(q + selector_count for q in components) % prime == 0)
            print(json.dumps({'modulus_bits': modulus.bit_length(), 'prime': prime,
                              'selector_count': selector_count, 'zero': zero,
                              'seconds': round(perf_counter() - started, 6),
                              'work': dict(work)}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    check()
    if not args.check_only:
        benchmark()


if __name__ == '__main__':
    main()
