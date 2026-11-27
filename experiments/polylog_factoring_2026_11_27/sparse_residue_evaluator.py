"""Evaluate one zero-residue bit of the idempotent count on semiprimes.

For N=pq with distinct primes p<q and public prime ell, ell divides
(p+A)(q+A) exactly when p is -A or N/(-A) modulo ell, provided A and N
are units modulo ell. Test the corresponding arithmetic progressions up
to floor(sqrt(N)). Multiplying by ell**(-1) modulo N turns each into an
ordinary interval, so the existing Pollard–Strassen product/remainder
trees apply. A nonunit in these intervals can only come from p, because
the unscaled entries are smaller than q. Check p's residue if a multiple
of p generated the hit. No nonunit certifies a negative bit directly.

The charged work is soft-O(1+N**(1/4)/sqrt(ell)) times polynomial bit
factors. This saves work for sparse large-prime residue tests. It does
not improve the exponent for small ell or evaluate the full residue.
The CLI uses only synthetic semiprimes.
"""

import argparse
import json
import sys
from math import gcd, isqrt
from pathlib import Path

from sympy import nextprime

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.metrics import Metrics
from src.trackC.baselines import ps_interval


def zero_residue_bit(modulus, prime, selector_count, backend='flint'):
    """Return the bit and counters, given promised distinct-prime semiprime N."""
    if modulus < 6 or prime < 2:
        raise ValueError('Expected a semiprime modulus and a prime answer modulus')
    work = Metrics()

    def from_divisor(divisor):
        return work.to_record(
            zero=(divisor + selector_count) * (modulus // divisor + selector_count) % prime == 0,
            factor_recovered=True,
        )

    if modulus % 2 == 0:
        return from_divisor(2)
    work.inc('gcd')
    common = gcd(modulus, prime)
    if common != 1:
        return from_divisor(common)

    root = -selector_count % prime
    if root == 0 or prime == 2:
        return work.to_record(zero=root != 0, factor_recovered=False)
    other_root = modulus * pow(root, -1, prime) % prime
    inverse_step = pow(prime, -1, modulus)
    work.inc('inv', 2)
    limit = isqrt(modulus)
    for start in sorted({root, other_root}):
        if start > limit:
            work.inc('empty_progressions')
            continue
        length = (limit - start) // prime + 1
        normalized_start = start * inverse_step % modulus
        work.inc('progressions')
        work.inc('represented_terms', length)
        work.hi('largest_progression', length)
        result = ps_interval(modulus, normalized_start, normalized_start + length,
                             m=work, backend=backend)
        if result['nonunit']:
            divisor = result['divisor']
            if not 1 < divisor <= limit or modulus % divisor:
                raise ValueError('Progression probe did not return a valid smaller prime')
            return from_divisor(divisor)
    return work.to_record(zero=False, factor_recovered=False)


def check(backend):
    cases = [(2, 11), (3, 53), (101, 107), (103, 109), (101, 1009), (503, 937),
             (int(nextprime(2**17 + 1234)), int(nextprime(2**17 + 4321)))]
    checked = negative_without_factor = 0
    for p, q in cases:
        for prime in (2, 3, 5, 7, 11, 31, 101, 1009):
            for selector_count in range(min(prime, 10)):
                result = zero_residue_bit(p*q, prime, selector_count, backend)
                expected = (p + selector_count) * (q + selector_count) % prime == 0
                assert result['zero'] == expected
                checked += 1
                negative_without_factor += not expected and not result['factor_recovered']

    # Here 3p is in a tested progression even though p and q miss -A.
    result = zero_residue_bit(101 * 1009, 7, 5, backend)
    assert result['factor_recovered'] and not result['zero']
    print(f'{checked} zero bits matched, including {negative_without_factor} '
          'negative results without factor recovery')


def benchmark(backend):
    p, q = (int(nextprime(2**34 + offset)) for offset in (12345, 43210))
    modulus = p * q
    primes = [int(nextprime(2**exponent)) for exponent in (35, 30, 25, 20, 15, 10)] + [1009]
    for prime in primes:
        positive = -p % prime
        negative = next(a for a in range(1, prime) if (p+a)*(q+a) % prime != 0)
        for selector_count in (positive, negative):
            result = zero_residue_bit(modulus, prime, selector_count, backend)
            assert result['zero'] == ((p+selector_count)*(q+selector_count) % prime == 0)
            print(json.dumps({'modulus_bits': modulus.bit_length(), 'prime': prime,
                              'selector_count': selector_count, **result}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend', default='flint', choices=('flint', 'py'))
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    check(args.backend)
    if not args.check_only:
        benchmark(args.backend)


if __name__ == '__main__':
    main()
