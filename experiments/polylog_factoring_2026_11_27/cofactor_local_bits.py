"""Which balanced cofactor-bit answers follow just from a product modulo a power of two?"""

from math import isqrt
from pathlib import Path
import sys

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.ellipse_disk_mod8.native_dyadic import floor_sum


# Exactly the singleton high bits that cannot occur on both odd units.
_FORCED_NO = {1: (3,), 2: (5, 7), 3: (11, 13), 4: (19, 29)}


def _modulus(n, mask):
    if n < 1 or not n & 1 or mask < 0 or mask & 1:
        raise ValueError("expected an odd positive product and a mask using only positions >=1")
    return 1 << max(1, mask.bit_length())


def local_outcomes(n, mask):
    """Boolean outcomes on all odd unit pairs a*b = n modulo the full bit period."""
    modulus = _modulus(n, mask)
    if not mask:
        return frozenset((False,))
    j = mask.bit_length()-1
    if mask == 1 << j and n % modulus in _FORCED_NO.get(j, ()):
        return frozenset((False,))
    return frozenset((False, True))


def local_witness(n, mask, answer):
    """Construct a compatible pair of residues, or return None for an impossible answer."""
    modulus = _modulus(n, mask)
    if not answer:
        return 1, n % modulus
    if modulus <= 128:
        for a in range(1, modulus, 2):
            b = n*pow(a, -1, modulus) % modulus
            if a & mask and b & mask:
                return a, b
        return None

    # In each of these two fibers the inverse map on the high digits is affine.
    step = 1 << ((modulus.bit_length())//2)
    length = modulus//step
    half = length//2
    for base in (1, 3):
        slope = -n*pow(base, -2, length) % length
        offset = (n*pow(base, -1, modulus) % modulus)//step + slope*half

        def hits(count):
            return (floor_sum(count, length, slope, offset+half)
                    - floor_sum(count, length, slope, offset))

        if hits(half) == 0:
            continue
        lo, hi = 1, half
        while lo < hi:
            mid = (lo+hi)//2
            if hits(mid):
                hi = mid
            else:
                lo = mid+1
        a = base+step*(half+lo-1)
        return a, n*pow(a, -1, modulus) % modulus
    raise ArithmeticError("the high-bit fibers did not produce the guaranteed witness")


def screened_outcomes(n, mask):
    """Sound answers from local and real relaxations on a promised balanced semiprime."""
    if __package__:
        from .cofactor_arch_bits import arch_outcomes
    else:
        from cofactor_arch_bits import arch_outcomes
    _modulus(n, mask)
    possible_bits = (1 << isqrt(2*n).bit_length())-1
    effective = mask & possible_bits
    return local_outcomes(n, effective) & arch_outcomes(n, effective)
