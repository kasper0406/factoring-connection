"""Native gcd-prefix estimates on a promised product of two primes."""

from fractions import Fraction
from math import isqrt


def exact_two_prime_prefix(p, q, t):
    if p == q:
        return t + (p-1)*(t//p) + (p*p-p)*(t//(p*p))
    return t + (p-1)*(t//p) + (q-1)*(t//q) + (p-1)*(q-1)*(t//(p*q))


def semiprime_prefix(n, t):
    """One-third relative accuracy for all t >= 0, under the semiprime promise."""
    if n < 4 or t < 0:
        raise ValueError("Expected semiprime n and nonnegative endpoint")
    root = isqrt(n)
    if n % 2 == 0:
        return exact_two_prime_prefix(2, n//2, t)
    if root*root == n:
        return exact_two_prime_prefix(root, root, t)
    periods, remainder = divmod(t, n)
    if remainder <= root:
        tail = Fraction(4*remainder, 3)
    else:
        tail = 2*remainder-1
    return periods*4*n + tail


def balanced_short_prefix(n, t):
    """Return estimate and certified relative error for p<q<2p and 0<=t<=sqrt(n)."""
    from sympy import integer_nthroot
    if __package__:
        from .quartic_collision_coarse_split import interval_divisor
    else:
        from quartic_collision_coarse_split import interval_divisor

    if not 0 <= t <= isqrt(n):
        raise ValueError("Expected a cutoff at most the integer square root")
    if 2*t*t <= n:
        return t, Fraction(0)

    width = max(1, int(integer_nthroot(n, 6)[0]) // 128)
    first = max(t-width+1, isqrt(n//2)+1)
    q_low = (n+t-1)//t
    q_high = n//first+1
    if q_low < q_high:
        divisor = interval_divisor(n, q_low, q_high, 4)
        if divisor is not None:
            return exact_two_prime_prefix(divisor, n//divisor, t), Fraction(0)

    # The largest possible p <= t is now t-width. If first was clipped,
    # the entire feasible p-range was checked.
    if first == isqrt(n//2)+1:
        return t, Fraction(0)
    upper = 2*t-width-1
    return Fraction(2*t*upper, t+upper), Fraction(upper-t, upper+t)


def degree_count(j, m):
    """Number of (a,b) in [0,j]^2 with a+b <= m."""
    def triangle(s):
        return (s+1)*(s+2)//2 if s >= 0 else 0
    return triangle(m) - 2*triangle(m-j-1) + triangle(m-2*j-2)


def degree_limit(base, cap, threshold):
    value, degree = 1, 0
    if threshold < 1:
        return -1
    while degree < cap and value <= threshold // base:
        value *= base
        degree += 1
    return degree


def balanced_power_bounds(n, j, t):
    """Bound the gcd prefix; outside the critical band the upper/lower ratio is <= 1.9."""
    if n < 15 or n % 2 == 0 or j < 1 or t < 0:
        raise ValueError("Expected a balanced odd semiprime, positive power and nonnegative endpoint")
    base = isqrt(n//2)
    ceiling = isqrt(2*n) + 1
    onset = base + 1
    if t < onset:
        return t, t
    if t < 2*onset:
        if t <= isqrt(n):
            return t, 2*t
        return Fraction(3*t-1, 2), 3*t-2
    if t < base*base:
        return 2*t-1, 3*t

    m = degree_limit(base, 2*j, t)
    v = degree_limit(ceiling, 2*j, t//32)
    possible = degree_count(j, m)
    if m == 2 and t < n:
        possible -= 1 + int(j >= 2)  # Exclude pq and, when present, q**2.
    upper = possible*t
    lower = Fraction((base-1)**2, base**2) * (
        Fraction(31*t, 32)*degree_count(j, v)
        + (Fraction(t, 2) if v < 2 and t >= n else 0)
    )
    lower = max(t, lower)
    if 10*upper <= 19*lower:
        return lower, upper

    # The degree-count bounds certify all base >= 32*3**19 (see the proof).
    # Only fixed-size bases reach this exact fallback.
    p = next(d for d in range(3, isqrt(n)+1, 2) if n % d == 0)
    q = n//p
    total = 0
    pa = 1
    for a in range(degree_limit(p, j, t)+1):
        qb = 1
        for b in range(degree_limit(q, j, t//pa)+1):
            phi_p = pa-pa//p if a else 1
            phi_q = qb-qb//q if b else 1
            total += phi_p*phi_q*(t//(pa*qb))
            qb *= q
        pa *= p
    return total, total


def balanced_power_prefix(n, j, t):
    """Approximate the gcd prefix on n**j within 1/3, given n=pq with odd p<q<2p."""
    lower, upper = balanced_power_bounds(n, j, t)
    return Fraction(2*upper*lower, upper+lower) if upper else 0


def padded_power_prefix(n, j, padding_factors, t):
    """One-third approximation on s*n**j with known prime powers of s coprime to n."""
    if __package__:
        from .gcd_prefix_compression import divisor_mass_jumps
    else:
        from gcd_prefix_compression import divisor_mass_jumps

    if any(n % p == 0 for p in padding_factors):
        raise ValueError("Padding must be coprime to the base")
    base_low, base_high = balanced_power_bounds(n, j, t)
    if not padding_factors or t == 0:
        return Fraction(2*base_low*base_high, base_low+base_high) if t else 0

    onset = isqrt(n//2) + 1
    precision = 8
    while precision <= 512*(len(padding_factors)+1):
        jumps = divisor_mass_jumps(padding_factors, precision)
        if t < 2*onset:
            tail = sum(w*(t//d) for d,w in jumps if d > 1)
            lower = base_low + tail
            upper = base_high + Fraction(precision+4, precision)*tail
        else:
            lower = upper = 0
            for d, w in jumps:
                if d > t:
                    break
                lo, hi = balanced_power_bounds(n, j, t//d)
                lower += w*lo
                upper += w*hi
            upper *= Fraction(precision+1, precision)
        if upper <= 2*lower:
            return Fraction(2*upper*lower, upper+lower)
        precision *= 2
    raise ArithmeticError("Compression did not reach the proved enclosure")
