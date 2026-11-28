"""Turn a supplied partial multiple of r−1 or r+1 into an affine divisor hint.

The public input promises a balanced N=pq. If g=gcd(e,r±1), q'=N/r, the
divisor g*q' of e*N divides N±q'. The ordinary affine reader on e*N recovers
the gcd at that offset, whose intersection with N is exactly q'. Neither the
order sign, r, nor the missing cofactor (r±1)/g is an algorithm input.
"""

from math import gcd, isqrt

from sympy import integer_nthroot

if __package__:
    from .affine_prime_window import _parameters, affine_divisors
else:
    from affine_prime_window import _parameters, affine_divisors


def factor_from_affine_order_hint(n, e, *, a=1, b=1, stats=None):
    """Return a verified proper divisor or None. Assumes distinct p<q<=2p.

    This decoder targets 1<e<n. For each fixed epsilon>0, e<=n^(1−epsilon)
    and missing cofactor n^o(1) at either hidden prime suffice asymptotically.
    With nonzero a,b of size n^o(1), the local orders become |a*r±b|. An
    inverse for b modulo n is unnecessary: gcd(b,n) is checked directly.
    An individual call reports the attained exact divisor threshold in stats.
    """
    if n < 6 or e < 1 or not a or not b:
        raise ValueError('expected n>=6, e>=1 and nonzero a,b')
    if stats is None:
        stats = {}
    stats.update(method='range', dimension=0, radius=0, lower=None, offsets_tested=0)
    for value in (e, a, b):
        common = gcd(n, value)
        if 1 < common < n:
            stats['method'] = 'shared factor'
            return common
    if e == 1 or e >= n or b % n == 0:
        return None

    modulus = n*e
    max_prime = isqrt(2*n-1)+1
    radius = abs(b)*max_prime
    stats['radius'] = radius
    # With s=log_n(n/e)>=epsilon, the target has exponent 1−s/32+o(1).
    # Its square strictly exceeds (2−s)/2, the affine threshold for this radius.
    shrink = int(integer_nthroot(n//e, 32)[0])
    target = max(2, (n-max_prime+shrink-1)//shrink)
    if e*max_prime < target:
        stats['method'] = 'small exponent'
        return None
    # Every possible d,m requires log_M(sqrt(2)*radius) < log_M(target)^2.
    # After omitting sqrt(d), this follows by completing the square in
    # beta=log_M(target): d(d-1)*beta^2 - 2*d*m*beta + m(m+1) >= 0.
    # One-sixteenth-bit enclosures can establish impossibility without an LLL.
    log_modulus = (modulus**16).bit_length()-1
    log_radius = (radius**16).bit_length()-1 + 8
    log_target = (target**16-1).bit_length()
    if log_radius*log_modulus >= log_target**2:
        stats['method'] = 'no certificate'
        return None
    try:
        d, m = _parameters(modulus, radius, target)
    except ValueError:
        stats['method'] = 'no certificate'
        return None
    # At that same dimension, use the least lower bound certified by the
    # determinant. This can only increase the number of covered cofactors.
    det4 = modulus**(2*m*(m+1)) * radius**(2*d*(d-1))
    lhs = (1 << (d*(d-1))) * det4 * d**(2*d)
    lower = int(integer_nthroot(lhs, 4*m*d)[0])+1
    stats['lower'] = lower
    for divisor in affine_divisors(modulus, a*n, radius, lower, stats=stats):
        common = gcd(divisor, n)
        if 1 < common < n:
            return common
    return None
