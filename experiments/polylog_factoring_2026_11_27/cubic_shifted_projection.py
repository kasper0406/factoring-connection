"""Unpadded phase extraction from shifted critical cubic projections."""

from fractions import Fraction
from math import gcd
from secrets import randbelow
from sympy import jacobi_symbol


def aligned_projection(n, character_argument, scale, evaluate):
    """Set -3v to the character of the argument, and the CRT phase to scale."""
    modulus = n*n
    v = -3*character_argument
    a = -v*v*pow(4*scale, -1, modulus) % modulus
    return evaluate(modulus, a, 2 + modulus*v)


def quadratic_residue_semiprime(n, unit, evaluate):
    """Decide quadratic residuosity on distinct semiprimes coprime to 6 (unit input)."""
    if gcd(n, 6*unit) != 1:
        raise ValueError("Expected a unit and a semiprime coprime to six")
    if jacobi_symbol(unit, n) == -1:
        return False
    for scale in (1, 2):
        difference = (aligned_projection(n, 1, scale, evaluate)
                      - aligned_projection(n, unit, scale, evaluate))
        if difference.real > 1:
            return False
    return True


def circular_distance(left, right):
    offset = (left-right) % 1
    return min(offset, 1-offset)


def try_shifted_split(n, evaluate, draw=randbelow):
    """One trial, for exactly two distinct prime-power components coprime to six.

    evaluate(M,a,B) estimates the complex projected cubic to absolute error <=1/32.
    """
    if n < 5 or gcd(n, 6) != 1:
        raise ValueError("Expected base coprime to six")
    mask = 1 + draw(n-1)
    divisor = gcd(mask, n)
    if divisor != 1:
        return divisor
    def point_at(scale):
        base = aligned_projection(n, 1, scale, evaluate)
        masked = aligned_projection(n, mask, scale, evaluate)
        return (base-masked)/2 - 2

    return divisor_from_phases(n, point_at)


def divisor_from_phases(n, point_at):
    modulus = n*n
    observations = []
    scale = 1
    for _ in range(modulus.bit_length()+1):
        point = point_at(scale)
        coordinates = (point.real, point.imag, -point.real, -point.imag)
        observations.append(Fraction(max(range(4), key=coordinates.__getitem__), 4))
        scale = 2*scale % modulus

    angle = observations[-1]
    for observation in reversed(observations[:-1]):
        half = angle/2
        angle = min((half, half+Fraction(1, 2)), key=lambda x: circular_distance(x, observation))
    projector = round(modulus*angle) % modulus
    divisor = gcd(projector, n)
    return divisor if 1 < divisor < n else None


def divisor_from_cosines(n, cosine_at):
    """Decode a signed CRT projector from 0.01-accurate cosines at unit scales."""
    from mpmath import mp

    modulus = n*n
    dyadic = {}
    scale = 1
    for _ in range(modulus.bit_length()+1):
        dyadic[scale] = max(-1, min(1, cosine_at(scale)))
        scale = 2*scale % modulus
    reference = next((s for s,c in dyadic.items() if abs(c) <= mp.mpf("0.9")), None)
    if reference is None:
        return None
    ref_cos = dyadic[reference]
    ref_sin = mp.sqrt(1-ref_cos**2)
    points = {}
    for scale, cosine in dyadic.items():
        direction = 1
        shifted = (scale+reference) % modulus
        divisor = gcd(shifted, n)
        if divisor == n:
            direction = -1
            shifted = (scale-reference) % modulus
            divisor = gcd(shifted, n)
        if 1 < divisor < n:
            return divisor
        shifted_cos = dyadic.get(shifted)
        if shifted_cos is None:
            shifted_cos = max(-1, min(1, cosine_at(shifted)))
        sine = direction*(cosine*ref_cos-shifted_cos)/ref_sin
        points[scale] = mp.mpc(cosine, sine)
    return divisor_from_phases(n, points.__getitem__)
