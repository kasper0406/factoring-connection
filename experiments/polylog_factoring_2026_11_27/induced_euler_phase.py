"""Decode a split prime from noisy imaginary-axis samples of the induced Euler quotient."""

from math import isqrt

from mpmath import mp

if __package__:
    from .affine_prime_window import affine_divisors
else:
    from affine_prime_window import affine_divisors


def _invert_near(center, radius, observed, multiple):
    """Fit sin(multiple*angle) on the certified monotone interval about center."""
    lo, hi = center-radius, center+radius
    increasing = mp.cos(multiple*center) > 0
    for _ in range(18):
        mid = (lo+hi)/2
        if (mp.sin(multiple*mid) < observed) == increasing:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def factor_from_phase(n, euler_ratio, *, stats=None, round_only=False):
    """Use ratios with complex absolute error <=1/256 on distinct balanced odd primes.

    euler_ratio(t) approximates L(it,psi)/L(it,chi); chi is positive at the
    positive-character prime. Rounding returns that prime, while the lattice
    variant returns a verified factor if both primes are inside one interval.
    """
    if n < 15 or n % 2 == 0:
        raise ValueError("expected an odd balanced semiprime")
    work = {} if stats is None else stats
    with mp.workprec(n.bit_length()+75):
        half_log = mp.log(n)/2
        sine = {}

        def sample(scale):
            if scale not in sine:
                rotation = mp.expj(scale*half_log)
                # rotation*(R-1+exp(-it log n)) has imaginary part 2 sin(t h).
                sine[scale] = mp.im(rotation*(euler_ratio(scale)-1) + rotation.conjugate())/2
            return sine[scale]

        scale = 1
        error = mp.mpf(1)/32
        angle = _invert_near(mp.mpf(0), mp.log(2)/2, sample(1), 1)
        while True:
            estimate = mp.exp(half_log+angle/scale)
            center = int(mp.nint(estimate))
            if round_only:
                if scale >= isqrt(2*n)+1:
                    if 1 < center < n and n % center == 0:
                        work.update(samples=len(sine), max_height=max(sine), radius=0)
                        return center
                    raise ArithmeticError("the rounded estimate is not a factor")
            else:
                # |r-estimate| <= sqrt(2n)*(exp(error/scale)-1).
                upper = isqrt(2*n)+1
                radius = (upper+16*scale-1)//(16*scale) + 1
                # Extra margin permits smaller practical lattice dimensions.
                if radius <= max(2, isqrt(isqrt(n))//128):
                    lower = isqrt((n-1)//2)+1
                    possible = affine_divisors(n, center, radius, lower)
                    factors = [v for v in possible if v < n]
                    work.update(samples=len(sine), max_height=max(sine), radius=radius)
                    if factors:
                        return min(factors, key=lambda v: abs(v-center))
                    raise ArithmeticError("the decoded interval contains no factor")
            scale *= 2
            predicted = 2*angle
            multiple = 1 if abs(mp.cos(predicted)) >= mp.mpf("0.5") else 2
            angle = _invert_near(predicted, 2*error, sample(multiple*scale), multiple)
