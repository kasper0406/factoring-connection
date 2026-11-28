"""Certified readings of the four-square transform on Re(s)=3 and Re(s)=4."""

from fractions import Fraction
from math import ceil, isqrt

from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import to_rational

if __package__:
    from .four_square_prefix import strip_small_prime_powers
else:
    from four_square_prefix import strip_small_prime_powers


def _known_magnitude(factors, a, offset, radial_power, bits):
    # A private interval context avoids changing the precision of other users.
    ctx = MPIntervalContext()
    ctx.prec = bits
    magnitude = ctx.mpf(1)
    angle_scale = ctx.mpf(offset.numerator) / offset.denominator
    if a is not None:
        angle_scale += ctx.pi / ctx.log(ctx.mpf(a))
    for p, exponent in factors:
        r = ctx.mpf(p)
        angle = angle_scale * ctx.log(r)
        radial = r ** radial_power
        real_step = ctx.cos(angle) / radial
        imag_step = -ctx.sin(angle) / radial
        # Horner evaluation of z + z² + ... + z^e, z=p^(2-s).
        real = imag = ctx.mpf(0)
        for _ in range(exponent):
            shifted = 1 + real
            real, imag = (real_step * shifted - imag_step * imag,
                          imag_step * shifted + real_step * imag)
        weight = 1 + 1 / r
        real, imag = 1 + weight * real, weight * imag
        magnitude *= ctx.sqrt(real ** 2 + imag ** 2)
    # Read the outward-rounded binary endpoints as exact rationals, with no
    # intermediate float or decimal string conversion.
    return tuple(Fraction(*to_rational(bound)) for bound in magnitude._mpi_)


def right_line_magnitude(n, a, epsilon, *, offset=0, real_part=3, stats=None):
    """Enclose |D_n(s)| in rationals L,U with U ≤ (1+epsilon)L.

    The imaginary part is pi/log(a) + offset for integer a >= 2, or just
    offset for a=None. The offset is rational, real_part is 3 or 4, n is
    positive and odd, and 0 < epsilon <= 1.
    An optional mutable `stats` dictionary receives the trial-division work,
    the untreated cofactor and its tail bound, and the interval precision.
    """
    if not isinstance(n, int) or n < 1 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer")
    if a is not None and (not isinstance(a, int) or a < 2):
        raise ValueError("a must be None or an integer at least 2")
    if not isinstance(real_part, int) or real_part not in (3, 4):
        raise ValueError("the real part must be 3 or 4")
    offset = Fraction(offset)
    radial_power = real_part - 2
    eps = Fraction(epsilon)
    if not 0 < eps <= 1:
        raise ValueError("epsilon must be positive and at most one")
    cutoff = ceil(16 * n.bit_length() / eps)
    if radial_power == 2:
        cutoff = 1 + isqrt(cutoff - 1)
    factors, residual, trials = strip_small_prime_powers(n, cutoff)
    # Every prime left in the residual exceeds the cutoff. Each of its local
    # factors has distance at most (1+1/p)/(p^radial_power-1) from 1.
    prime_bound = (residual.bit_length() - 1) // ((cutoff + 1).bit_length() - 1)
    tail = Fraction(2 * prime_bound, cutoff ** radial_power)

    # At frequency zero all phases are +1. At pi/log(a) the p=a phase is -1.
    exact = Fraction(1)
    oscillatory = []
    for p, exponent in factors:
        if offset == 0 and (a is None or p == a):
            step = Fraction(1 if a is None else -1, p ** radial_power)
            exact *= 1 + (1 + Fraction(1, p)) * step * (1 - step ** exponent) / (1 - step)
        else:
            oscillatory.append((p, exponent))
    bits = 0
    low = high = Fraction(1)
    if oscillatory:
        offset_bits = max(0, abs(offset.numerator).bit_length() - offset.denominator.bit_length())
        bits = max(48, eps.denominator.bit_length() - eps.numerator.bit_length() + 32) + offset_bits
        while True:
            low, high = _known_magnitude(oscillatory, a, offset, radial_power, bits)
            if low > 0 and high <= (1 + eps / 4) * low:
                break
            bits *= 2

    # For total deviation T<1 the omitted product's magnitude is between
    # 1-T and 1/(1-T). Here T<=eps/8, so the final width meets the budget.
    low *= exact * (1 - tail)
    high *= exact / (1 - tail)
    if not 0 < low <= high <= (1 + eps) * low:
        raise ArithmeticError("the rational enclosure missed its width target")
    if stats is not None:
        stats.update(trial_candidates=trials, cutoff=cutoff,
                     known_prime_powers=tuple(factors), residual=residual,
                     unknown_prime_count_bound=prime_bound, tail_sum_bound=tail,
                     precision_bits=bits)
    return low, high


def right_line_squared_magnitude(n, a, absolute_error, *, offset=0, real_part=3, stats=None):
    """Approximate |D_n(s)|² by a rational with certified additive error.

    A coarse upper bound sets the relative precision for a second call if
    the first squared interval does not already meet the absolute error.
    """
    tolerance = Fraction(absolute_error)
    if tolerance <= 0:
        raise ValueError("absolute_error must be positive")
    options = dict(offset=offset, real_part=real_part)
    info = {} if stats is not None else None
    low, high = right_line_magnitude(n, a, 1, **options, stats=info)
    pilot_trials = 0
    if high * high - low * low > 2 * tolerance:
        if info is not None:
            pilot_trials = info['trial_candidates']
        # (1+eps)^2-1 <= 3eps for eps<=1. The coarse upper bound still
        # bounds the true magnitude and hence the refined lower bound.
        eps = min(Fraction(1), tolerance / (3 * high * high))
        low, high = right_line_magnitude(n, a, eps, **options, stats=info)
    lower_square, upper_square = low * low, high * high
    if upper_square - lower_square > 2 * tolerance:
        raise ArithmeticError("the squared enclosure missed its absolute error target")
    if stats is not None:
        stats.update(info, squared_interval=(lower_square, upper_square),
                     total_trial_candidates=info['trial_candidates'] + pilot_trials)
    return (lower_square + upper_square) / 2
