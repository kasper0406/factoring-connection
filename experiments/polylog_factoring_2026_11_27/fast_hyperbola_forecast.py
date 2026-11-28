"""Stationary-phase approximation of the scaled, unnormalized hyperbola bit integral."""

from mpmath import mp


def fast_scaled_forecast(n, r, period, h=1, *, cutoff=15, digits=20, stats=None):
    """Approximate sqrt(lambda)*integral w(x,z)*(H(lambda*x)*H(lambda*z/x)-1/4) dx.

    z=n/r², lambda=h*r/period, H is the upper-half periodic indicator, and
    w=log(r)²/[x log(r*x) log(r*z/x)]. The Fourier cutoff is a separate
    asymptotic approximation. digits controls the absolute numerical precision.
    """
    if r < 4 or not r*r <= n <= 4*r*r or period < 1 or h < 1:
        raise ValueError("expected r>=4, r²<=n<=4r², and positive period and multiplier")
    if cutoff < 1 or digits < 1:
        raise ValueError("cutoff and digits must be positive integers")
    ctx = mp.clone()
    # Retain all input bits and guard the rapidly varying phase, endpoints,
    # and sum of O(cutoff²) terms. Precision grows with the input bit lengths.
    ctx.prec = 4*digits + n.bit_length() + h.bit_length() + period.bit_length() + 6*cutoff.bit_length() + 48
    z = ctx.mpf(n)/r**2
    lam = ctx.mpf(h)*r/period
    low, high = max(ctx.one, z/2), min(ctx.mpf(2), z)
    work = {} if stats is None else stats
    work.update(terms=0, ordered_modes=0, fresnel_evaluations=0, limiting_values=0, precision_bits=ctx.prec)
    if low == high:
        return ctx.zero
    log_r = ctx.log(r)
    limit = ctx.mpc(1, 1)/4
    # The exact Fresnel tail is at most 1/(2*pi*s); the extra factor covers
    # weighted accumulation without relying on cancellation among the modes.
    saturation = ctx.power(10, digits+4)*(cutoff+1)**2

    def fresnel(s):
        if abs(s) > saturation:
            work["limiting_values"] += 1
            return limit if s > 0 else -limit
        work["fresnel_evaluations"] += 1
        return ctx.mpc(ctx.fresnelc(2*s), ctx.fresnels(2*s))/2

    def endpoint(x, ratio, scale):
        # Rationalization retains a small signed distance near a boundary saddle.
        return scale*(ratio*x*x-z)/(x*(ctx.sqrt(ratio*x)+ctx.sqrt(z/x)))

    total = ctx.zero
    for m in range(1, cutoff+1, 2):
        # The substitution x -> z/x interchanges (m,ell) and (ell,m).
        # It reverses and negates the Fresnel endpoints; the integral is equal.
        for ell in range(m, min(cutoff, 4*m)+1, 2):
            ratio = ctx.mpf(m)/ell
            center = ctx.sqrt(z/ratio)
            amplitude = log_r**2/(ctx.sqrt(ratio*center)*ctx.log(r*center)*ctx.log(r*z/center))
            scale = ctx.sqrt(lam*ell)
            first = endpoint(low, ratio, scale)
            last = endpoint(high, ratio, scale)
            delta = fresnel(last)-fresnel(first)
            phase = ctx.expjpi(4*lam*ctx.sqrt(m*ell*z))
            multiplicity = 1 if m == ell else 2
            total -= multiplicity*amplitude*ctx.re(phase*delta)/(m*ell*ctx.sqrt(ell))
            work["terms"] += 1
            work["ordered_modes"] += multiplicity
    return 2*total/ctx.pi**2
