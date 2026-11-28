"""Continue squared Mellin magnitudes to reflected products at heights 1 and 3."""

from dataclasses import dataclass
from fractions import Fraction
from math import factorial, isqrt
from numbers import Rational

if __package__:
    from .mellin_right_line import right_line_magnitude, right_line_squared_magnitude
else:
    from mellin_right_line import right_line_magnitude, right_line_squared_magnitude


@dataclass(frozen=True)
class ContinuationEstimate:
    """A complex rational center and a certified absolute error in the plane."""

    real: Fraction
    imag: Fraction
    error: Fraction
    samples: int
    requested_accuracy: Fraction

    def magnitude_bounds(self):
        """Rational bounds for the absolute value of the true product."""
        scale = (self.error.denominator + self.error.numerator - 1) // self.error.numerator
        squared = (self.real**2 + self.imag**2) * scale**2
        floor_norm = isqrt(squared.numerator // squared.denominator)
        lower = max(Fraction(), Fraction(floor_norm, scale) - self.error)
        upper = Fraction(floor_norm + 1, scale) + self.error
        return lower, upper


def _weights(m, height=1):
    """Lagrange weights at i*height for samples at 0,1/m,...,(m-1)/m."""
    target = m * height
    real, imag = 1, 0
    for k in range(m):
        real, imag = -k * real - target * imag, target * real - k * imag

    # Divide the common Gaussian integer ∏(target*i-k) by its j-th factor.
    # ∏_{k!=j}(j-k) = (-1)^(m-1-j) j! (m-1-j)!.
    denominator = (-1)**(m - 1) * factorial(m - 1)
    for j in range(m):
        norm = j * j + target * target
        numerator_real = (target * imag - j * real) // norm
        numerator_imag = (-target * real - j * imag) // norm
        yield Fraction(numerator_real, denominator), Fraction(numerator_imag, denominator)
        if j + 1 < m:
            denominator = -(denominator * (j + 1) // (m - 1 - j))


def _check_modulus(n):
    if not isinstance(n, int) or n < 1 or n % 2 == 0:
        raise ValueError("expected a positive odd modulus")


def _interpolate(n, height, budget, sample_squared):
    # For F(z)=D_n(3+i(t+z)) D_n(3-i(t+z)), expansion into exponentials
    # gives |F^(m)| <= n^(2+height) b^m on 0 <= Im(z) <= height, b=bits(n).
    # Frequencies are log(e/d), d,e|n. Their magnitudes are at most log n,
    # and the absolute coefficients at height zero sum to D_n(3)^2 <= n².
    # The divided-difference simplex integral gives the remainder bound
    # sup |F^(m)|/m! times ∏ |height*i-j/m|, with each distance <= height+1.
    growth = n**(2 + height)
    derivative_scale = (height + 1) * n.bit_length()
    m = 1
    remainder = Fraction(growth * derivative_scale)
    while remainder > budget:
        m += 1
        remainder *= Fraction(derivative_scale, m)

    # The sum of the absolute weights is at most
    # (2 sqrt(1+height²) m)^(m-1)/(m-1)!; bases 10 and 24 suffice.
    amplification = (10 if height == 1 else 24)**m
    delta = budget / amplification
    real = imag = Fraction()
    for j, (weight_real, weight_imag) in enumerate(_weights(m, height)):
        value = sample_squared(Fraction(j, m), delta)
        if not isinstance(value, Rational):
            raise TypeError("squared-magnitude samples must be rational")
        real += weight_real * value
        imag += weight_imag * value
    return ContinuationEstimate(real, imag, remainder + amplification * delta, m, delta)


def interpolate_mellin_product(n, epsilon, sample_squared):
    """Return a certified disk around D_n(2+it) D_n(4-it) for a fixed real t.

    At rational x, sample_squared(x, delta) must return a rational within
    absolute delta of |D_n(3+i(t+x))|². No phase observations are needed.
    The returned error is at most epsilon/(50*n).
    """
    _check_modulus(n)
    if not isinstance(epsilon, Rational) or epsilon <= 0:
        raise ValueError("epsilon must be a positive rational")
    return _interpolate(n, 1, Fraction(epsilon) / (100 * n), sample_squared)


def exponent_from_right_line(n, sample_squared):
    """Recover the exact integer n^5 D_n(0) D_n(6) from t=0 squared readings.

    sample_squared(x, delta) must approximate |D_n(3+ix)|² by a rational
    within absolute error delta. The extrapolation uses O(log(2n)) samples
    of O(log(2n)) bits of absolute accuracy. The public right-line evaluator's
    cost depends on 1/delta; a source polynomial in the precision bit count
    remains unknown.
    """
    _check_modulus(n)
    scale = n**5
    product = _interpolate(n, 3, Fraction(1, 8 * scale), sample_squared)
    # For d|n, n^5 ψ(d²)/d^6=(n/d)^5 ψ(d) is integral, as is D_n(0).
    # Scaling the certified disk at 3i leaves radius at most 1/4.
    return round(scale * product.real)


def critical_magnitude_bounds(n, a, epsilon, *, stats=None):
    """Enclose |D_n(2+i*pi/log(a))| within relative width 1+epsilon.

    Uses public right-line evaluations. Their running time depends on the
    inverse sample tolerance, so the end-to-end cost can be exponential.
    """
    if not isinstance(a, int) or a < 2:
        raise ValueError("expected an integer a at least two")
    if not isinstance(epsilon, Rational) or not 0 < epsilon <= 1:
        raise ValueError("epsilon must be rational and in (0,1]")
    eps = Fraction(epsilon)
    trials = 0
    precision = 0

    def reading(offset, delta):
        nonlocal trials, precision
        info = {} if stats is not None else None
        value = right_line_squared_magnitude(n, a, delta, offset=offset, stats=info)
        if info is not None:
            trials += info["total_trial_candidates"]
            precision = max(precision, info["precision_bits"])
        return value

    product = interpolate_mellin_product(n, eps, reading)
    product_low, product_high = product.magnitude_bounds()
    info = {} if stats is not None else None
    fourth_low, fourth_high = right_line_magnitude(n, a, eps / 8, real_part=4, stats=info)
    lower, upper = product_low / fourth_high, product_high / fourth_low
    # On the critical line |D_n| >= 1/n. At real part four the omitted
    # distances from one sum to at most ∑_{k>=3} 1/[k(k-1)] = 1/2.
    # Hence the product is at least 1/(2n). The product interval has
    # width at most 3*product.error <= 3*eps/(50n).
    if not 0 < lower <= upper <= (1 + eps) * lower:
        raise ArithmeticError("the continued enclosure missed the relative width target")
    if stats is not None:
        stats.update(samples=product.samples, requested_accuracy=product.requested_accuracy,
                     product_error=product.error, total_trial_candidates=trials + info["trial_candidates"],
                     max_precision_bits=max(precision, info["precision_bits"]))
    return lower, upper
