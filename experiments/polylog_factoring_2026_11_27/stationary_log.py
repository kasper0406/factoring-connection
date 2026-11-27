"""Check factor-revealing frequencies of the transformed full stationary trace.

The fixtures supply prime powers to synthesize a trace and verify its CRT
identity. Dense spectral recovery on the small grids uses only trace values.
"""

import json
from cmath import exp
from itertools import product
from math import cos, gcd, log1p, pi, prod
from time import perf_counter


def crt_parts(prime_powers):
    modulus = prod(prime_powers)
    parts = []
    for q in prime_powers:
        complement = modulus // q
        inverse = pow(complement, -1, q)
        parts.append((q, complement * inverse))
    return parts


def stationary_roots(modulus, parts):
    return [
        sum(sign * e for sign, (_, e) in zip(signs, parts)) % modulus
        for signs in product((-1, 1), repeat=len(parts))
    ]


def root_trace(modulus, roots, t):
    return sum(cos(2 * pi * ((2 * r * t) % modulus) / modulus) for r in roots)


def product_trace(modulus, parts):
    factors = []
    for q, e in parts:
        inverse = e // (modulus // q)
        factors.append([2 * cos(2 * pi * ((2 * inverse * r) % q) / q) for r in range(q)])
    return [prod(table[t % len(table)] for table in factors) for t in range(modulus)]


def kloosterman(n, parameters):
    """Literal complete K(1,t²;N²), evaluated only on small fixtures."""
    square = n * n
    units = [(x, pow(x, -1, square)) for x in range(1, square) if gcd(x, n) == 1]
    phases = [exp(2j * pi * a / square) for a in range(square)]
    return [sum(phases[(x + t * t * inverse) % square] for x, inverse in units)
            for t in parameters]


def transformed_trace(raw):
    return [log1p(abs(value)) for value in raw]


def spectrum_at(values, frequencies):
    modulus = len(values)
    phases = [exp(-2j * pi * a / modulus) for a in range(modulus)]
    return [
        abs(sum(value * phases[(k * t) % modulus] for t, value in enumerate(values))) / modulus
        for k in frequencies
    ]


def benchmark(prime_powers, exponent):
    n = prod(prime_powers)
    modulus = n**exponent
    started = perf_counter()
    parts = crt_parts([q**exponent for q in prime_powers])
    roots = stationary_roots(modulus, parts)
    assert len(set(roots)) == 2**len(parts)
    assert all(r * r % modulus == 1 for r in roots)
    raw = product_trace(modulus, parts)
    crt_error = max(abs(value - root_trace(modulus, roots, t)) for t, value in enumerate(raw))
    assert crt_error < 1e-11
    fixture_seconds = perf_counter() - started

    direct_parameters = []
    kloosterman_error = None
    if exponent == 2 and n <= 45:
        direct_parameters = sorted({0, 1, 2, 3, 5, n - 1, n, n + 1, modulus - 1})
        direct = kloosterman(n, direct_parameters)
        expected = [raw[t] if gcd(t, n) == 1 else 0 for t in direct_parameters]
        kloosterman_error = max(abs(value / n - target) for value, target in zip(direct, expected))
        assert kloosterman_error < 1e-11

    clean = transformed_trace(raw)
    # Quantizing the signed raw trace keeps an absolute error of at most 1/128.
    rounded = [round(value * 64) / 64 for value in raw]
    coarse = transformed_trace(rounded)
    raw_error = max(abs(a - b) for a, b in zip(raw, rounded))
    log_error = max(abs(a - b) for a, b in zip(clean, coarse))
    assert raw_error <= 1 / 128 and log_error <= raw_error + 1e-14

    # Multiples of N have gcd N and are supplied as zero in the reduction.
    # At other nonunits a gcd already splits N. Leaving these trace values in
    # the diagnostic checks spectral recovery without using that early exit.
    for t in range(0, modulus, n):
        coarse[t] = 0.0

    started = perf_counter()
    local_frequencies = [min(4 * e % modulus, -4 * e % modulus) for _, e in parts]
    clean_peaks = spectrum_at(clean, local_frequencies)
    coarse_peaks = spectrum_at(coarse, local_frequencies)
    assert min(coarse_peaks) > 0.02
    local_divisors = [gcd(k, n) for k in local_frequencies]
    assert all(divisor == n // q for q, divisor in zip(prime_powers, local_divisors))
    local_seconds = perf_counter() - started

    dense = None
    if modulus <= 2_025:
        started = perf_counter()
        frequencies = range(1, modulus // 2 + 1)
        peak, frequency = max(zip(spectrum_at(coarse, frequencies), frequencies))
        divisor = gcd(frequency, n)
        assert 1 < divisor < n
        dense = {
            "frequency": frequency,
            "peak": round(peak, 8),
            "divisor": divisor,
            "products": modulus * len(frequencies),
            "seconds": round(perf_counter() - started, 6),
        }

    return {
        "n": n,
        "modulus": modulus,
        "stationary_roots": len(roots),
        "crt_root_sum_error": crt_error,
        "literal_kloosterman_queries": len(direct_parameters),
        "literal_kloosterman_error": kloosterman_error,
        "raw_quantization_error": raw_error,
        "transformed_quantization_error": log_error,
        "min_clean_local_coefficient": round(min(clean_peaks), 8),
        "min_coarse_local_coefficient": round(min(coarse_peaks), 8),
        "local_frequencies": local_frequencies,
        "local_divisors": local_divisors,
        "known_frequency_products": 2 * modulus * len(parts),
        "fixture_seconds": round(fixture_seconds, 6),
        "known_frequency_seconds": round(local_seconds, 6),
        "dense": dense,
    }


if __name__ == "__main__":
    fixtures = [((3, 5), 1), ((9, 5), 2), ((3, 5), 2), ((9, 25), 2),
                ((9, 5, 7), 1), ((9, 25, 7), 1), ((3, 5, 7, 11), 1)]
    for components, exponent in fixtures:
        print(json.dumps(benchmark(components, exponent), sort_keys=True))
