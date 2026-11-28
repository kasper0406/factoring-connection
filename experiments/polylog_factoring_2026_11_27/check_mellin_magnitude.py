"""Factor-aware interval arithmetic checks for the real Mellin magnitude."""

from dataclasses import dataclass
from decimal import Context, Decimal
from fractions import Fraction
from functools import cache
from math import isqrt, prod
from random import Random

if __package__:
    from . import mellin_magnitude as decoder
    from .mellin_right_line import right_line_magnitude
else:
    import mellin_magnitude as decoder
    from mellin_right_line import right_line_magnitude


SCALE = 1 << 256


def ceildiv(n, d):
    return -((-n) // d)


@dataclass(frozen=True)
class Interval:
    """Endpoints in units of 2^-256. All arithmetic rounds outwards."""

    lo: int
    hi: int

    @classmethod
    def rationals(cls, a, b):
        a, b = Fraction(a), Fraction(b)
        return cls(a.numerator * SCALE // a.denominator, ceildiv(b.numerator * SCALE, b.denominator))

    def __add__(self, other):
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other):
        return Interval(self.lo - other.hi, self.hi - other.lo)

    def __mul__(self, other):
        products = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(products) // SCALE, ceildiv(max(products), SCALE))

    def divide(self, positive):
        assert positive.lo > 0
        reciprocal = Interval(SCALE**2 // positive.hi, ceildiv(SCALE**2, positive.lo))
        return self * reciprocal

    def over_integer(self, n):
        assert n > 0
        return Interval(self.lo // n, ceildiv(self.hi, n))

    def square(self):
        smallest = 0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2)
        return Interval(smallest // SCALE, ceildiv(max(self.lo**2, self.hi**2), SCALE))

    def sqrt(self):
        assert self.lo >= 0
        high = isqrt(self.hi * SCALE)
        return Interval(isqrt(self.lo * SCALE), high + (high**2 != self.hi * SCALE))


@cache
def log_integer(n):
    # Decimal.ln is correctly rounded (nearest, ties to even). The adjacent
    # representables enclose the exact logarithm even at a decimal boundary.
    context = Context(prec=90)
    value = context.ln(Decimal(n))
    return Interval.rationals(context.next_minus(value), context.next_plus(value))


def arctan_inverse(n):
    term, result = Fraction(1, n), Fraction()
    j = 0
    while True:
        result += term / (2 * j + 1)
        j += 1
        term /= -(n * n)
        next_term = term / (2 * j + 1)
        if abs(next_term) * SCALE < 1:
            return Interval.rationals(min(result, result + next_term), max(result, result + next_term))


@cache
def pi():
    # Machin's identity. Both arctangents have alternating, decreasing terms.
    return Interval.rationals(16, 16) * arctan_inverse(5) - Interval.rationals(4, 4) * arctan_inverse(239)


def sine(angle):
    turns = round(Fraction(angle.lo + angle.hi, pi().lo + pi().hi))
    angle = angle - integer(turns) * pi()
    assert max(abs(angle.lo), abs(angle.hi)) <= 2 * SCALE
    square = angle.square()
    term = result = angle
    j = 0
    while True:
        j += 1
        magnitude = (term * square).over_integer((2 * j) * (2 * j + 1))
        term = Interval(-magnitude.hi, -magnitude.lo)
        bound = max(abs(term.lo), abs(term.hi))
        if bound <= 2:
            # Throughout [-2,2] the sine series error is at most the first
            # omitted term in absolute value, also when the input crosses 0.
            result = Interval(result.lo - bound, result.hi + bound)
            return Interval(-result.hi, -result.lo) if turns % 2 else result
        result = result + term


def integer(value):
    return Interval(value * SCALE, value * SCALE)


def local_magnitude(r, log_a):
    """|1+1/r+exp(i*pi*log(r)/log(a))| via the small phase difference."""
    angle = (pi() * (log_integer(r) - log_a)).divide(log_a).over_integer(2)
    inverse = integer(1).over_integer(r)
    return (inverse.square() + integer(4) * (integer(1) + inverse) * sine(angle).square()).sqrt()


def target_bounds(p, q, a):
    n = p * q
    log_a = log_integer(a)
    product = local_magnitude(p, log_a) * local_magnitude(q, log_a)
    scale = n * n.bit_length()**2
    lo, hi = Fraction(product.lo * scale, SCALE), Fraction(product.hi * scale, SCALE)
    assert lo > 0 and hi - lo < lo / 10**18
    return lo, hi


def lower_endpoint(n):
    a = isqrt(n // 2)
    return a + (2 * a * a < n)


def check_comparison(p, q, a, lo, hi):
    n = p * q
    b = n.bit_length()
    wp = b + 10 * abs(a - p)
    wq = b + 10 * abs(a - q)
    assert Fraction(wp * wq, 4) <= lo <= hi <= Fraction(9 * wp * wq, 4), (p, q, a)
    if a <= p:
        ceiling_root = isqrt(n) + (isqrt(n)**2 != n)
        proxy = b + 10 * (ceiling_root - a)
        assert proxy <= wq <= Fraction(5 * proxy, 2)


def check_center(p, q):
    n = p * q
    # At t=2*pi/log(N), the factor-blind central approximation is zero.
    log_center = log_integer(n).over_integer(2)
    magnitude = local_magnitude(p, log_center) * local_magnitude(q, log_center)
    limit = Fraction(16, n.bit_length()**2) + Fraction(7, isqrt(n))
    assert 0 <= magnitude.lo <= magnitude.hi and Fraction(magnitude.hi, SCALE) <= limit


def right_line_control(factors, a):
    exact = Fraction(1)
    value = integer(1)
    for r, exponent in factors:
        if r == a:
            # Sum the alternating geometric series term by term, independently
            # of the simplified prime-power formula in the public routine.
            exact *= 1 + (1 + Fraction(1, r)) * sum(Fraction(-1, r)**j for j in range(1, exponent + 1))
            continue
        assert exponent == 1
        amplitude = (integer(1) + integer(1).over_integer(r)).over_integer(r)
        phase = (pi() * (log_integer(r) - log_integer(a))).divide(log_integer(a)).over_integer(2)
        size = ((integer(1) - amplitude).square() + integer(4) * amplitude * sine(phase).square()).sqrt()
        value = value * size
    return exact * Fraction(value.lo, SCALE), exact * Fraction(value.hi, SCALE)


def check_right_line():
    rough = [(671101003, 1), (939530897, 1)]
    large = prod(p for p, _ in rough)
    fixtures = [
        ([], 2, 1),
        ([(3, 1), (5, 1), (13, 1)], 17, 1),
        (rough, 900000000, large),
        ([(3, 2), *rough], 3, large),
        ([(5, 1), *rough], 131, large),
        ([(21990232567889, 1), (30786325584527, 1)], 21990232567889, 21990232567889 * 30786325584527),
    ]
    fixtures.extend(([(7, e)], 7, 1) for e in (1, 2, 3, 7))
    checked = retained = 0
    for factors, a, expected_residual in fixtures:
        n = prod(p**e for p, e in factors)
        ref_lo, ref_hi = right_line_control(factors, a)
        for tolerance in (Fraction(1), Fraction(1, 4), Fraction(1, 100)):
            stats = {}
            lower, upper = right_line_magnitude(n, a, tolerance, stats=stats)
            assert 0 < lower <= ref_lo <= ref_hi <= upper <= (1 + tolerance) * lower
            assert stats["residual"] == expected_residual
            assert prod(p**e for p, e in stats["known_prime_powers"]) * expected_residual == n
            if not factors or (len(factors) == 1 and factors[0][0] == a):
                assert lower == upper == ref_lo == ref_hi
            retained += expected_residual > 1
            checked += 1
    return checked, retained


def check_grid(p, q):
    n = p * q
    first, last = lower_endpoint(n), isqrt(n)
    # Cover both sides of the lower factor, including the exact cancellation
    # of the oscillating pair up to the remaining 1/p at a=p.
    if last - first < 300:
        positions = range(first, last + 1)
    else:
        positions = {first, p, last}
        positions.update(first + j * (last - first) // 32 for j in range(33))
        positions.update(a for a in (p - 1, p + 1, p - n.bit_length(), p + n.bit_length())
                         if first <= a <= last)
    for a in positions:
        check_comparison(p, q, a, *target_bounds(p, q, a))
    exact_local = local_magnitude(p, log_integer(p))
    assert Fraction(exact_local.lo, SCALE) <= Fraction(1, p) <= Fraction(exact_local.hi, SCALE)
    check_center(p, q)
    return len(positions)


def check_decoder(p, q, mode, rng, distortion=2):
    n = p * q
    queries = 0

    def oracle(modulus, a):
        nonlocal queries
        assert modulus == n and lower_endpoint(n) <= a <= p
        lo, hi = target_bounds(p, q, a)
        check_comparison(p, q, a, lo, hi)
        # Every value between hi/K and K*lo satisfies the promised error
        # for every real value inside [lo,hi]. Try both boundaries and noise.
        left, right = hi / distortion, distortion * lo
        assert left <= right
        if mode == "low":
            response = left
        elif mode == "high":
            response = right
        else:
            response = left + (right - left) * Fraction(rng.randrange(1025), 1024)
        assert hi / distortion <= response <= distortion * lo
        queries += 1
        return response

    assert decoder.factor_from_mellin_magnitude(n, oracle, distortion) == p
    assert queries > 0
    return queries


def main():
    rng = Random(20261128)
    fixtures = [
        (257, 509), (271, 277), (283, 401),
        (671101003, 939530897), (1000000007, 1000000009),
        (21990232567889, 30786325584527),
    ]
    grid = calls = recovered = 0
    for p, q in fixtures:
        grid += check_grid(p, q)
        if (p * q).bit_length() >= 60:
            for mode in ("low", "high", "random"):
                calls += check_decoder(p, q, mode, rng)
                recovered += 1
    for distortion in (Fraction(3, 2), 3):
        p, q = fixtures[-1]
        calls += check_decoder(p, q, "random", rng, distortion)
        recovered += 1
    right_checks, retained = check_right_line()
    print(f"{grid} certified comparison points, {len(fixtures)} central bounds, "
          f"{recovered} verified splits from {calls} bounded-error magnitude queries, "
          f"{right_checks} public right-line enclosures ({retained} retaining rough cofactors)")


if __name__ == "__main__":
    main()
