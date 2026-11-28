"""Check continuation against independent interval sums over the supplied divisors."""

from fractions import Fraction
from math import ceil, floor, isqrt
from random import Random

from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import to_rational

if __package__:
    from .check_four_square_incidence import prime_powers
    from .check_four_square_prefix import divisor_weights
    from .check_norm_one_exponent import split_norm_one_exponent
    from .mellin_continuation import (
        critical_magnitude_bounds, exponent_from_right_line, interpolate_mellin_product,
    )
    from .mellin_right_line import right_line_magnitude, right_line_squared_magnitude
else:
    from check_four_square_incidence import prime_powers
    from check_four_square_prefix import divisor_weights
    from check_norm_one_exponent import split_norm_one_exponent
    from mellin_continuation import (
        critical_magnitude_bounds, exponent_from_right_line, interpolate_mellin_product,
    )
    from mellin_right_line import right_line_magnitude, right_line_squared_magnitude


def endpoints(interval):
    return tuple(Fraction(*to_rational(bound)) for bound in interval._mpi_)


class DivisorSum:
    """Use cosines and sines for every divisor, independently of local Horner factors."""

    def __init__(self, n):
        self.terms = divisor_weights(prime_powers(n))
        self.ctx = MPIntervalContext()
        self.logs = {}

    def prepare(self, bits):
        if self.ctx.prec != bits or not self.logs:
            self.ctx.prec = bits
            self.logs = {d: self.ctx.log(self.ctx.mpf(d)) for d, _ in self.terms}

    def rational(self, q):
        q = Fraction(q)
        return self.ctx.mpf(q.numerator) / q.denominator

    def value(self, sigma, t):
        real = imag = self.ctx.mpf(0)
        for d, weight in self.terms:
            angle = t * self.logs[d]
            amplitude = self.ctx.mpf(weight) / d**sigma
            real += amplitude * self.ctx.cos(angle)
            imag -= amplitude * self.ctx.sin(angle)
        return real, imag

    def noisy_reader(self, frequency):
        # Keep alternating errors near the opposite edges of the promised
        # interval. Precision includes guard bits for all the small fixtures.
        calls = []

        def read(x, delta):
            self.prepare(delta.denominator.bit_length() + 80)
            real, imag = self.value(3, self.rational(frequency + x))
            low, high = endpoints(real**2 + imag**2)
            assert high - low < delta / 4
            direction = 1 if len(calls) % 2 else -1
            calls.append(x)
            return (low + high) / 2 + direction * 3 * delta / 4

        return read, calls


def check_native_critical():
    fixtures = ((1, 2), (27, 3), (81, 5), (105, 5), (11 * 17, 13), (1009 * 1231, 1009))
    eps = Fraction(1, 8)
    for n, a in fixtures:
        control = DivisorSum(n)
        control.prepare(256)
        real, imag = control.value(2, control.ctx.pi / control.ctx.log(control.ctx.mpf(a)))
        true_low, true_high = endpoints(control.ctx.sqrt(real**2 + imag**2))
        info = {}
        low, high = critical_magnitude_bounds(n, a, eps, stats=info)
        assert 0 < low <= true_low <= true_high <= high <= (1 + eps) * low
        if n == 1:
            assert true_low == true_high == 1
    return len(fixtures), info


def check_complex_interpolation():
    eps = Fraction(1, 9)
    unit = interpolate_mellin_product(1, eps, lambda x, delta: 1)
    assert (unit.real, unit.imag) == (1, 0)

    fixtures = ((9, 2), (9, -2), (81, Fraction(1, 2)), (105, 1), (187, Fraction(3, 2)))
    quadrants = set()
    total_calls = 0
    for n, t in fixtures:
        control = DivisorSum(n)
        reader, calls = control.noisy_reader(t)
        disk = interpolate_mellin_product(n, eps, reader)
        assert len(calls) == disk.samples
        total_calls += disk.samples
        control.prepare(256)
        x, y = control.value(2, control.rational(t))
        u, v = control.value(4, control.rational(-t))
        real, imag = endpoints(x * u - y * v), endpoints(x * v + y * u)
        real_error = max(abs(z - disk.real) for z in real)
        imag_error = max(abs(z - disk.imag) for z in imag)
        assert real_error**2 + imag_error**2 <= disk.error**2
        assert disk.error <= eps / (50 * n)
        assert real[0] * real[1] > 0 and imag[0] * imag[1] > 0
        quadrants.add((real[0] > 0, imag[0] > 0))
    assert len(quadrants) >= 3
    return len(fixtures), total_calls


def sigma_square(n):
    square = n * n
    return sum(d + (square // d if d * d != square else 0)
               for d in range(1, isqrt(square) + 1) if square % d == 0)


def exact_exponent(n, terms):
    value = n**5 * sigma_square(n) * sum(Fraction(weight, d**6) for d, weight in terms)
    assert value.denominator == 1
    return value.numerator


def check_exponents():
    inputs = [*range(1, 42, 2), 45, 77]
    noisy_calls = native_calls = 0
    splits = []
    for n in inputs:
        control = DivisorSum(n)
        if n in (45, 77):
            def native(x, delta):
                nonlocal native_calls
                native_calls += 1
                return right_line_squared_magnitude(n, None, delta, offset=x)

            value = exponent_from_right_line(n, native)
            divisor = split_norm_one_exponent(n, 3, value, Random(20261128 + n))
            assert 1 < divisor < n and n % divisor == 0
            splits.append((n, divisor))
        else:
            reader, calls = control.noisy_reader(0)
            value = exponent_from_right_line(n, reader)
            noisy_calls += len(calls)
        assert value == exact_exponent(n, control.terms)
    return len(inputs), noisy_calls, native_calls, splits


def check_zero_frequency():
    inputs = [n for n in range(1, 100, 2) if all(e == 1 for _, e in prime_powers(n))]
    inputs.extend((3 * 5 * 7, 3 * 5 * 7 * 11, 17 * 23))
    pairs = 0
    for n in inputs:
        low, high = right_line_magnitude(n, None, Fraction(1, 8 * n**3))
        exact = sigma_square(n)
        assert ceil(n**2 * low) == exact == floor(n**2 * high)
        factors = prime_powers(n)
        if len(factors) == 2:
            s = sum(p for p, _ in factors)
            assert exact == n**2 + n * s + s**2 - n + s + 1
            pairs += 1
    # On a repeated prime the same public measurement gives a different integer.
    low, high = right_line_magnitude(9, None, Fraction(1, 8 * 9**3))
    assert ceil(81 * low) == 129 == floor(81 * high) and sigma_square(9) == 121
    return len(inputs), pairs


if __name__ == "__main__":
    critical, info = check_native_critical()
    complex_cases, product_calls = check_complex_interpolation()
    exponents, noisy_calls, native_calls, splits = check_exponents()
    squarefree, pairs = check_zero_frequency()
    print(f"{critical} native critical enclosures; N=1242079 used {info['samples']} samples, "
          f"{info['max_precision_bits']} working bits, {info['total_trial_candidates']} trial candidates")
    print(f"{complex_cases} noisy complex enclosures from {product_calls} samples plus exact N=1")
    print(f"{exponents} exact integer exponents from {noisy_calls} injected and {native_calls} native samples; "
          f"native exponent splits {splits}")
    print(f"{squarefree} exact squarefree zero-frequency counts, including {pairs} semiprime identities")
