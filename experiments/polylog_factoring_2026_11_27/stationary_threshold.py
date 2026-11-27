"""Test fixed magnitude and doubling comparisons on the complete trace."""

import argparse
import json
from math import cos, fsum, gcd, log, prod, tau
from random import Random

from stationary_log import crt_parts, spectrum_at


def logarithm_tables(modulus, parts):
    tables = []
    for q, idempotent in parts:
        slope = 2 * (idempotent // (modulus // q))
        tables.append([log(abs(2 * cos(tau * ((slope * r) % q) / q))) for r in range(q)])
    return tables


def comparison_tables(logs, mode):
    if mode == 'threshold':
        return logs
    return [[values[(2 * r) % len(values)] - value for r, value in enumerate(values)]
            for values in logs]


def worst_covariance_bit(total, mean, gap):
    if total > gap:
        return 1
    if total < -gap:
        return 0
    # On undecided inputs choose the bit that minimizes (total − mean) * bit.
    return int(total < mean)


def divisor_from_frequency(n, modulus, frequency):
    return gcd(n, modulus // gcd(modulus, frequency))


def dense_fixture(components, mode, gap):
    n = prod(components)
    modulus = n * n
    parts = crt_parts([q * q for q in components])
    tables = comparison_tables(logarithm_tables(modulus, parts), mode)
    mean = fsum(log(2) / q for q, _ in parts) if mode == 'threshold' else 0.0
    totals = [fsum(values[t % len(values)] for values in tables) for t in range(modulus)]
    assert abs(fsum(totals) / modulus - mean) < 1e-12
    bits = [worst_covariance_bit(total, mean, gap) for total in totals]
    covariance = fsum((total - mean) * bit for total, bit in zip(totals, bits)) / modulus
    assert covariance > 0.1
    # At every proper nonunit a gcd probe could stop before the comparison.
    # Keep those bits to test the spectral route without the early success.
    for t in range(0, modulus, n):
        bits[t] = 0
    frequencies = range(1, modulus // 2 + 1)
    ranked = sorted(zip(spectrum_at(bits, frequencies), frequencies), reverse=True)
    for rank, (peak, frequency) in enumerate(ranked, 1):
        divisor = divisor_from_frequency(n, modulus, frequency)
        if 1 < divisor < n:
            assert peak > 0.03
            return {'n': n, 'mode': mode, 'gap': gap, 'covariance': covariance,
                    'ambiguous': sum(abs(total) <= gap for total in totals),
                    'frequency_rank': rank, 'peak': peak, 'divisor': divisor}
    raise AssertionError('No factor-revealing frequency was found')


def sampled_covariance(components, mode, gap, samples, rng):
    n = prod(components)
    parts = crt_parts([q * q for q in components])
    tables = comparison_tables(logarithm_tables(n * n, parts), mode)
    means = [log(2) / q if mode == 'threshold' else 0.0 for q, _ in parts]
    mean = fsum(means)
    totals = [0.0] * len(tables)
    ambiguous = 0
    for _ in range(samples):
        values = [rng.choice(table) for table in tables]
        total = fsum(values)
        ambiguous += abs(total) <= gap
        if worst_covariance_bit(total, mean, gap):
            for i, value in enumerate(values):
                totals[i] += value - means[i]
    covariances = [value / samples for value in totals]
    assert fsum(covariances) > 0.1
    return {'prime_components': len(parts), 'mode': mode, 'gap': gap, 'samples': samples,
            'ambiguous_fraction': ambiguous / samples, 'total_covariance': fsum(covariances),
            'largest_local_covariance': max(covariances)}


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--gap', type=float, default=1 / 64, help='Symmetric gap in logarithmic units')
    cli.add_argument('--samples', type=int, default=32768)
    args = cli.parse_args()
    if not 0 <= args.gap <= 1 / 4 or args.samples < 1024:
        cli.error('use a gap in [0,1/4] and at least 1,024 samples')
    rng = Random(261127)
    for mode in ('threshold', 'doubling'):
        for powers in ((3, 5), (3, 7), (9, 5), (5, 7), (3, 5, 7)):
            print(json.dumps(dense_fixture(powers, mode, args.gap), sort_keys=True))
        for powers in ((9, 25, 7, 11, 13), (9, 25, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)):
            print(json.dumps(sampled_covariance(powers, mode, args.gap, args.samples, rng), sort_keys=True))
