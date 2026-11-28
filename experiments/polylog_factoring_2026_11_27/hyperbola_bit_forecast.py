"""Conditional upper-bit forecast from a product, its public factor size range and a period."""

from decimal import Decimal, localcontext
from fractions import Fraction


def _compare_intervals(z, r, step, h):
    """Intersect the two sorted lists of intervals where the selected bit is one."""
    count = h*r//step
    low = max(Fraction(r), Fraction(z, 2*r))
    high = min(Fraction(2*r), Fraction(z, r))
    cell = Fraction(step, 2*h)
    i, j = 2*count+1, 4*count-1
    while i < 4*count and j > 2*count:
        left_end = (i+1)*cell
        right_end = z/(j*cell)
        start = max(low, i*cell, z/((j+1)*cell))
        end = min(high, left_end, right_end)
        if start < end:
            yield start, end
        if left_end <= right_end:
            i += 2
        if right_end <= left_end:
            j -= 2


def upper_bit_forecast(n, r, period, h=1, *, digits=18):
    """Predict [hp mod period, hq mod period both lie in the upper half].

    The model draws the two factors continuously in [r,2r] with prime
    density 1/log(x), then conditions on their product n. Only n, r, period
    and h are used in evaluation. Work is proportional to h*r/period intervals.
    """
    if r < 2 or r & (r-1) or not r*r <= n <= 4*r*r:
        raise ValueError("expected a product in the given dyadic size range")
    if period < 2 or period & (period-1) or period > r or h < 1 or not h & 1:
        raise ValueError("expected a dyadic period at most r and a positive odd multiplier")
    if digits < 1:
        raise ValueError("expected positive precision")
    # Near each corner the whole feasible interval lies in the same half cell.
    if 2*h*(n-r*r) <= r*period:
        return Decimal(0)
    if h*(4*r*r-n) <= r*period:
        return Decimal(1)
    low = max(Fraction(r), Fraction(n, 2*r))
    high = min(Fraction(2*r), Fraction(n, r))
    count = h*r//period
    with localcontext() as ctx:
        # Guard cancellation in the mass of the whole feasible interval and
        # the sum over O(count) pieces. Unrounded work has logarithmic overhead.
        ctx.prec = digits + n.bit_length().bit_length() + count.bit_length() + 12
        log_n = Decimal(n).ln()

        def primitive(x):
            x_dec = Decimal(x.numerator)/Decimal(x.denominator)
            log_x = x_dec.ln()
            return (log_x/(log_n-log_x)).ln()

        total = primitive(high)-primitive(low)
        hit = sum((primitive(b)-primitive(a) for a,b in _compare_intervals(n,r,period,h)),
                  Decimal(0))
        ctx.prec = digits
        return +min(Decimal(1), max(Decimal(0), hit/total))
