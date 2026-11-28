"""Cubic Euler and full 3-power tests on norm-one residues of an odd semiprime."""

from math import gcd, lcm

if __package__:
    from .eisenstein_torus import ONE, OMEGA, norm, torus_power
else:
    from eisenstein_torus import ONE, OMEGA, norm, torus_power


def _compare_roots(n, z):
    """Return a proper component gcd if present, and whether z is a global cube root of 1."""
    for a, b in (ONE, OMEGA, (-1, -1)):
        common = gcd(n, z[0]-a, z[1]-b)
        if common > 1:
            return (None, True) if common == n else (common, False)
    return None, False


def torus_cubic_split(n, z, *, full_ladder=True, bound=1):
    """Scan one cubic Euler power, or remove small cofactors and scan the full 3-power ladder.

    Returns a verified divisor or None. The probability analysis assumes a
    uniform norm-one input, optionally conditioned on its public Jacobi mark.
    """
    if n <= 6 or gcd(n, 6) != 1 or norm(z) % n != 1:
        raise ValueError("expected an odd modulus coprime to 3 and a norm-one element")
    if bound < 1 or (not full_ladder and bound != 1):
        raise ValueError("the optional cofactor bound belongs to the full ladder")
    exponent = n - (1 if n % 3 == 1 else -1)
    if not full_ladder:
        return _compare_roots(n, torus_power(n, z, exponent//3))[0]
    while exponent % 3 == 0:
        exponent //= 3
    multiplier = 1
    for t in range(2, bound+1):
        if t % 3:
            multiplier = lcm(multiplier, t)
    value = torus_power(n, z, exponent*multiplier)
    for _ in range(n.bit_length()+1):
        factor, global_root = _compare_roots(n, value)
        if factor is not None or global_root:
            return factor
        following = torus_power(n, value, 3)
        if following == value:
            break
        value = following
    return None
