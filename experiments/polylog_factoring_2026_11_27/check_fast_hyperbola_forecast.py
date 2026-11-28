"""Compare the scaled fast forecast with jump enumeration and a high-frequency limit."""

from fractions import Fraction
from math import sqrt

import numpy as np
from mpmath import mp

if __package__:
    from .fast_hyperbola_forecast import fast_scaled_forecast
else:
    from fast_hyperbola_forecast import fast_scaled_forecast


def direct_scaled(n, r, period, h=1):
    """Literal jump cells; integrate their log-prime weight by its antiderivative."""
    z = n/r**2
    lam = h*r/period
    low, high = max(1, z/2), min(2, z)
    jumps = np.arange(int(2*lam), int(4*lam)+2)/(2*lam)
    knots = np.unique(np.clip(np.r_[low, high, jumps, z/jumps], low, high))
    mids = (knots[1:]+knots[:-1])/2
    hit = (np.floor(2*lam*mids).astype(np.int64) &
           np.floor(2*lam*z/mids).astype(np.int64) & 1)
    log_ratios = np.log(np.log(r*knots)/np.log(r*z/knots))
    integrals = np.diff(log_ratios)*(np.log(r)**2/np.log(n))
    return sqrt(lam)*np.sum(integrals*(hit-.25))


def check_moderate():
    r = 1 << 16
    results = {}
    for lam in (16, 64, 256, 1024, 4096):
        errors = []
        for z_numerator in (19, 25, 31, 37, 44, 58):
            n, period = r*r*z_numerator//16, r//lam
            exact = direct_scaled(n, r, period)
            row = [float(fast_scaled_forecast(n, r, period, cutoff=k, digits=12))
                   for k in (3, 9, 15)]
            errors.append(np.subtract(row, exact))
        rms = np.sqrt(np.mean(np.square(errors), axis=0))
        assert rms[2] < rms[0] and max(abs(e[2]) for e in errors) < .05
        results[lam] = tuple(round(float(error), 8) for error in rms)
    assert results[4096][-1] < .001
    # The model includes an odd scale h independently from the dyadic period.
    n = r*r*37//16
    odd = fast_scaled_forecast(n, r, r//128, 3, cutoff=15)
    assert abs(float(odd)-direct_scaled(n, r, r//128, 3)) < .01
    assert fast_scaled_forecast(r*r, r, 16) == fast_scaled_forecast(4*r*r, r, 16) == 0
    return results


def infinite_endpoint_reference(n, r, period, cutoff):
    """Independent ordered-mode sum replacing distant Fresnel endpoints by limits.

    Exact rational comparisons detect a saddle at an endpoint and give it
    half weight. This reference is used only at extremely large lambda.
    """
    ctx = mp.clone()
    ctx.prec = 300+2*n.bit_length()
    z = Fraction(n, r*r)
    low, high = max(1, z/2), min(2, z)
    result = ctx.zero
    for m in range(1, cutoff+1, 2):
        for ell in range(1, cutoff+1, 2):
            center_sq = z*ell/m
            if not low**2 <= center_sq <= high**2:
                continue
            fraction = ctx.mpf(".5") if center_sq in (low**2, high**2) else ctx.one
            center = ctx.sqrt(ctx.mpf(center_sq.numerator)/center_sq.denominator)
            partner = ctx.mpf(n)/(r*r*center)
            amplitude = ctx.log(r)**2/(ctx.sqrt(m*center)*ctx.log(r*center)*ctx.log(r*partner))
            cycles = 2*ctx.sqrt(m*ell*n)/period
            angle = 2*ctx.pi*(cycles-ctx.floor(cycles))+ctx.pi/4
            # For a full saddle F(infty)-F(-infty) = exp(i*pi/4)/sqrt(2).
            result -= fraction*ctx.sqrt(2)*amplitude*ctx.cos(angle)/(ctx.pi**2*m*ell)
    return result


def check_huge():
    cutoff = 9
    moderate = {}
    fast_scaled_forecast(37*(1 << 20)**2//16, 1 << 20, 1 << 10, cutoff=cutoff, stats=moderate)
    results = []
    for exponent, z in ((80, Fraction(37, 16)), (80, Fraction(5, 3)), (256, Fraction(37, 16))):
        r = 3*(1 << (exponent+4)) if z.denominator == 3 else 1 << (exponent+4)
        period = r >> exponent
        n = int(z*r*r)
        work = {}
        answer = fast_scaled_forecast(n, r, period, cutoff=cutoff, digits=18, stats=work)
        precise = fast_scaled_forecast(n, r, period, cutoff=cutoff, digits=35)
        reference = infinite_endpoint_reference(n, r, period, cutoff)
        assert abs(answer-precise) < 1e-17 and abs(answer-reference) < 1e-10
        assert work["terms"] == moderate["terms"] and work["ordered_modes"] == moderate["ordered_modes"]
        assert work["fresnel_evaluations"]+work["limiting_values"] == 2*work["terms"]
        if exponent == 256:
            assert work["limiting_values"] == 2*work["terms"]
        if z.denominator == 3:
            # Move across an exact boundary saddle by one integer on either side.
            for other in (n-1, n+1):
                assert abs(answer-fast_scaled_forecast(other, r, period, cutoff=cutoff)) < 1e-17
        results.append((exponent, str(z), float(answer), work["terms"], work["ordered_modes"],
                        work["limiting_values"]))
    return results


if __name__ == "__main__":
    print("RMS scaled errors at cutoffs 3,9,15:", check_moderate())
    print("Large lambda (exponent, z, scaled result, modes evaluated, ordered modes, limits):", check_huge())
