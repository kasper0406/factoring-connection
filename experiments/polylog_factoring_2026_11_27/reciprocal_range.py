"""Approximate one reciprocal quotient range with bounded quadratic roots."""

from argparse import ArgumentParser
from collections import Counter
from math import cos, fsum, gcd, isqrt, pi, sin, tau
from time import perf_counter

from lattice_window import block_square_roots, multiplicity


def phase(value, modulus):
    centered = (value + modulus // 2) % modulus - modulus // 2
    angle = tau * (centered / modulus)
    return complex(cos(angle), sin(angle))


def quotient_interval(modulus, quotient):
    """Positive centered representatives where round(modulus/x) = quotient."""
    low = 2 * modulus // (2 * quotient + 1) + 1
    high = min(modulus // 2, 2 * modulus // (2 * quotient - 1)) + 1
    return low, high


def digit_intersection(n, digit, low, high):
    first = (low - digit + n - 1) // n
    stop = (high - digit + n - 1) // n
    return first, stop - first


def lift_coefficient(n, c, digit, derivative, first, count):
    if derivative == 0:
        weight = count / n
    else:
        centered = (derivative + n // 2) % n - n // 2
        ratio = centered / n
        sine_ratio = sin(pi * ratio) / ratio if ratio else pi
        # Evaluate the reciprocal of N(1−e_N(d)) without a small subtraction.
        reciprocal = (0.5j * (1 / centered) / sine_ratio) * phase(-centered, 2 * n)
        weight = phase(derivative * first, n) * (1 - phase(derivative * count, n)) * reciprocal
    modulus = n * n
    inverse = pow(digit, -1, modulus)
    return (2 * c * inverse % n, phase(digit + c * c * inverse, modulus) * weight)


def scan_spectrum(n, c, start, width, low, high):
    """Sum the upper lifts analytically, with one visit per lower digit."""
    work = Counter()
    terms = {}
    for index in range(start, start + width):
        work['lower_digits'] += 1
        digit = index % n
        if gcd(digit, n) != 1:
            continue
        first, count = digit_intersection(n, digit, low, high)
        if count:
            derivative = (1 - (c * pow(digit, -1, n))**2) % n
            frequency, coefficient = lift_coefficient(n, c, digit, derivative, first, count)
            terms[frequency] = (derivative, coefficient)
            work['represented_lifts'] += count
    return terms, work


def truncated_spectrum(n, c, start, width, low, high, cutoff):
    """Compute exactly the coefficients at the retained derivative classes."""
    work = Counter()
    terms = {}
    limit = min(max(1, isqrt(n) // 8), width - 1)
    m = multiplicity(n, limit) if limit > 1 else None
    for d in range(-min(cutoff, n // 2), min(cutoff, n // 2) + 1):
        if gcd(1 - d, n) != 1:
            continue  # A nonunit times a square cannot equal the unit c².
        square = c * c * pow(1 - d, -1, n) % n
        for origin in range(start, start + width, limit + 1):
            stop = min(limit, start + width - 1 - origin)
            work['certified_blocks'] += 1
            if stop > 1:
                indices = block_square_roots(n, origin, stop, m, square)
                work['lattice_calls'] += 1
            else:
                indices = [i for i in range(origin, origin + stop + 1)
                           if (i * i - square) % n == 0]
            for index in indices:
                digit = index % n
                first, count = digit_intersection(n, digit, low, high)
                if count:
                    frequency, coefficient = lift_coefficient(n, c, digit, d, first, count)
                    terms[frequency] = (d % n, coefficient)
                    work['represented_lifts'] += count
    return terms, work


def value_at(n, digit, spectrum):
    terms = [coefficient * phase(k * digit, n) for k, (_, coefficient) in spectrum.items()]
    return complex(fsum(term.real for term in terms), fsum(term.imag for term in terms))


def direct_values(n, c, digit, start, width, low, high):
    """Enumerate literal terms and apply both inverse-reciprocity steps."""
    modulus = n * n
    a = (c + n * digit)**2 % modulus
    direct, reciprocal, twice, coarse, correction = 0j, 0j, 0j, 0j, 0.0
    for x in range(low, high):
        if (x - start) % n >= width or gcd(x, n) != 1:
            continue
        direct += phase(x + a * pow(x, -1, modulus), modulus)
        sign = 1 if x > 0 else -1
        v = abs(x)
        smooth = phase(sign * a, modulus * v)
        base = phase(x, modulus)
        first = phase(-sign * a * pow(modulus, -1, v), v)
        reciprocal += base * first * smooth
        coarse += base * first
        correction += tau * a / (modulus * v)
        q = (2 * modulus + v) // (2 * v)
        r = modulus - q * v
        if r == 0 or gcd(q, modulus) != 1:
            twice += base * first * smooth
            continue
        r_sign = 1 if r > 0 else -1
        denominator = abs(r)
        second = phase(sign * r_sign * q * a * pow(modulus, -1, denominator), denominator)
        second *= phase(-sign * r_sign * a, v * denominator)
        twice += base * second * smooth
    return direct / n, reciprocal / n, twice / n, coarse / n, correction / n


def truncation_rms(n, full, band, cutoff):
    expected = {k for k, (d, _) in full.items() if min(d, n - d) <= cutoff}
    assert band.keys() == expected
    assert all(abs(band[k][1] - full[k][1]) < 1e-12 for k in band)
    return fsum(abs(value)**2 for k, (_, value) in full.items() if k not in band)**0.5


def check():
    checks = 0
    for n in (15, 45, 77, 105, 315):
        for q in (2, 3, 7):
            left, right = quotient_interval(n * n, q)
            for low, high in ((left, right), (1 - right, 1 - left)):
                start, width = n - 4, min(n, 29)
                full, _ = scan_spectrum(n, 2, start, width, low, high)
                band, _ = truncated_spectrum(n, 2, start, width, low, high, 4)
                tail_rms = truncation_rms(n, full, band, 4)
                observed = fsum(abs(value_at(n, j, full) - value_at(n, j, band))**2
                                for j in range(n)) / n
                assert abs(observed - tail_rms**2) < 1e-12
                for j in (0, 1, n // 2):
                    direct, first, second, coarse, correction = direct_values(
                        n, 2, j, start, width, low, high,
                    )
                    assert abs(direct - value_at(n, j, full)) < 1e-9
                    assert abs(direct - first) < 1e-9
                    assert abs(direct - second) < 1e-9
                    assert abs(direct - coarse) <= correction + 1e-12
                    checks += 1
    print(f'{checks} direct quotient-range checks passed')


def benchmark(n, c, start, width, q, cutoff):
    low, high = quotient_interval(n * n, q)
    started = perf_counter()
    full, scan = scan_spectrum(n, c, start, width, low, high)
    scan_seconds = perf_counter() - started
    started = perf_counter()
    band, work = truncated_spectrum(n, c, start, width, low, high, cutoff)
    band_seconds = perf_counter() - started
    tail_rms = truncation_rms(n, full, band, cutoff)
    print(f'N={n}, q={q}, width={width}, cutoff={cutoff}, tail RMS={tail_rms:.6g}')
    print(f'Full scan {scan_seconds:.3f}s {dict(scan)}')
    print(f'Truncated evaluator {band_seconds:.3f}s {dict(work)}')
    print(f'Retained {len(band)} of {len(full)} lower phase coefficients')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('n', nargs='?', type=int)
    parser.add_argument('--scale', type=int, default=1)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--width', type=int)
    parser.add_argument('--quotient', type=int, default=2)
    parser.add_argument('--cutoff', type=int, default=8)
    parser.add_argument('--digit', type=int, help='Evaluate one upper digit using only the truncated spectrum')
    args = parser.parse_args()
    if args.n is None:
        check()
        for p, q in ((101, 109), (503, 937), (101**2, 109)):
            n = p * q
            mixed = 1 + p * (-2 * pow(p, -1, q) % q)
            width = n // 128
            benchmark(n, 1, mixed - width // 2, width, 2, 8)
    else:
        width = args.width if args.width is not None else max(1, args.n // 128)
        if (args.n < 5 or args.n % 2 == 0 or gcd(args.scale, args.n) != 1
                or not 0 < width <= args.n or args.quotient < 2 or args.cutoff < 0):
            parser.error('Expected odd n >= 5, unit scale, width in [1,n], quotient >= 2 and cutoff >= 0')
        if args.digit is None:
            benchmark(args.n, args.scale % args.n, args.start % args.n, width,
                      args.quotient, args.cutoff)
        else:
            low, high = quotient_interval(args.n * args.n, args.quotient)
            spectrum, work = truncated_spectrum(
                args.n, args.scale % args.n, args.start % args.n, width, low, high, args.cutoff,
            )
            print(f'Approximate normalized value {value_at(args.n, args.digit, spectrum)}')
            print(f'Truncated evaluator {dict(work)}')
