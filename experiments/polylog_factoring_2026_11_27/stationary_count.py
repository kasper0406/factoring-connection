"""Count idempotents using a cube-root bounded prime search and an unsplit residual."""

import argparse
import json
import sys
from math import gcd, isqrt, prod
from pathlib import Path

from sympy import integer_nthroot

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.stationary_dual import proven_prime
from src.common.metrics import Metrics
from src.specialfam.smoothcopp import smooth_part_ps


def perfect_power_base(n, work):
    for exponent in range(n.bit_length() - 1, 1, -1):
        work.inc('power_tests')
        root, exact = integer_nthroot(n, exponent)
        if exact:
            return root, exponent
    return n, 1


def count_stationary(n):
    """Return exact unweighted counts and the work of the bounded prime search."""
    if n < 1:
        raise ValueError('use a positive integer modulus')
    work = Metrics()
    base, exponent = perfect_power_base(n, work)
    cutoff = 0
    stripped = {}
    residual = base
    if base == 1:
        residual_count, residual_case = 0, 'empty'
    else:
        work.inc('primality_tests')
        if proven_prime(base):
            residual_count, residual_case = 1, 'prime'
        else:
            cutoff, _ = integer_nthroot(base, 3)
            work.inc('bounded_searches')
            stripped = smooth_part_ps(base, cutoff, work)
            residual //= prod(prime ** multiplicity for prime, multiplicity in stripped.items())
            # Each residual prime exceeds the original cube-root cutoff.
            # Three prime factors counted with multiplicity would exceed base.
            if residual == 1:
                residual_count, residual_case = 0, 'empty'
            elif isqrt(residual) ** 2 == residual:
                residual_count, residual_case = 1, 'prime-square'
            else:
                if residual != base:
                    work.inc('primality_tests')
                if residual != base and proven_prime(residual):
                    residual_count, residual_case = 1, 'prime'
                else:
                    residual_count, residual_case = 2, 'unsplit-two-prime'

    distinct = len(stripped) + residual_count
    return work.to_record(
        n=n, base=base, power=exponent, cutoff=cutoff, stripped_distinct=len(stripped),
        residual=residual, residual_case=residual_case, distinct_prime_count=distinct,
        idempotents=2 ** distinct,
        cubic_critical_roots=2 ** distinct if gcd(n, 6) == 1 else None,
        quartic_critical_roots=3 ** distinct if n % 2 else None,
    )


def check():
    for n in range(1, 513):
        result = count_stationary(n)
        assert result['idempotents'] == sum((value * value - value) % n == 0 for value in range(n))
        if result['cubic_critical_roots'] is not None:
            assert result['cubic_critical_roots'] == sum((3 * value * value - 3) % n == 0
                                                         for value in range(n))
        if result['quartic_critical_roots'] is not None:
            assert result['quartic_critical_roots'] == sum((4 * value ** 3 - 4 * value) % n == 0
                                                           for value in range(n))
    print('512 integer moduli checked by direct enumeration')

    for n, distinct, case, residual in (
        (10007 * 10009, 2, 'unsplit-two-prime', 10007 * 10009),
        (3 ** 2 * 5 ** 2 * 10007 * 10009, 4, 'unsplit-two-prime', 10007 * 10009),
        (3 * 10007 ** 2, 2, 'prime-square', 10007 ** 2),
        ((17 * 19) ** 6, 2, 'unsplit-two-prime', 17 * 19),
        (3 ** 5 * 7 ** 3, 2, 'empty', 1),
        (32416190071 * 32416187567, 2, 'unsplit-two-prime', 32416190071 * 32416187567),
    ):
        result = count_stationary(n)
        assert result['distinct_prime_count'] == distinct
        assert result['residual_case'] == case
        assert result['residual'] == residual
        print(json.dumps(result, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('modulus', nargs='?', type=int, help='evaluate from the modulus alone')
    args = parser.parse_args()
    if args.modulus is None:
        check()
        return
    try:
        print(json.dumps(count_stationary(args.modulus), sort_keys=True))
    except (RuntimeError, ValueError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
