"""Certified magnitude enclosure for the four-square transform at Re(s)=3."""

from fractions import Fraction
from math import ceil

from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import to_rational

if __package__:
    from .four_square_prefix import strip_small_prime_powers
else:
    from four_square_prefix import strip_small_prime_powers


def _known_magnitude(factors, a, bits):
    # A private interval context avoids changing the precision of other users.
    ctx = MPIntervalContext()
    ctx.prec = bits
    magnitude = ctx.mpf(1)
    angle_scale = ctx.pi / ctx.log(ctx.mpf(a))
    for p, exponent in factors:
        r = ctx.mpf(p)
        angle = angle_scale * ctx.log(r)
        real_step = ctx.cos(angle) / r
        imag_step = -ctx.sin(angle) / r
        # Horner evaluation of z + z² + ... + z^e, z=p^(-1-it).
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


def right_line_magnitude(n, a, epsilon, *, stats=None):
    """Return rationals L,U enclosing |D_n(3+i*pi/log(a))| with U ≤ (1+eps)L.

    `n` is positive and odd, `a` an integer at least 2, and 0 < epsilon <= 1.
    An optional mutable `stats` dictionary receives the trial-division work,
    the untreated cofactor and its tail bound, and the interval precision.
    """
    if not isinstance(n, int) or n < 1 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer")
    if not isinstance(a, int) or a < 2:
        raise ValueError("a must be an integer at least 2")
    eps = Fraction(epsilon)
    if not 0 < eps <= 1:
        raise ValueError("epsilon must be positive and at most one")
    cutoff = ceil(16 * n.bit_length() / eps)
    factors, residual, trials = strip_small_prime_powers(n, cutoff)
    # Every prime left in the residual exceeds the cutoff. Each of its local
    # factors has distance at most (1+1/p)/(p-1) <= 2/cutoff from 1.
    prime_bound = (residual.bit_length() - 1) // ((cutoff + 1).bit_length() - 1)
    tail = Fraction(2 * prime_bound, cutoff)

    # At p=a, the phase is exactly -1 and the whole local factor is rational.
    exact = Fraction(1)
    oscillatory = []
    for p, exponent in factors:
        if p == a:
            exact *= 1 - Fraction(1, p) + Fraction((-1) ** exponent, p ** (exponent + 1))
        else:
            oscillatory.append((p, exponent))
    bits = 0
    low = high = Fraction(1)
    if oscillatory:
        bits = max(48, eps.denominator.bit_length() - eps.numerator.bit_length() + 32)
        while True:
            low, high = _known_magnitude(oscillatory, a, bits)
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
