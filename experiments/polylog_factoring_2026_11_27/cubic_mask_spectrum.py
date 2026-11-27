"""Split using log-magnitude masks on the single unpadded square modulus."""

from math import gcd
from secrets import randbelow
from mpmath import mp
from sympy import primerange

if __package__:
    from .cubic_shifted_projection import aligned_projection, divisor_from_phases, divisor_from_cosines
else:
    from cubic_shifted_projection import aligned_projection, divisor_from_phases, divisor_from_cosines


def mask_spectrum(n, products, multiplier, scale, evaluate):
    values = [mp.log(abs(aligned_projection(n, multiplier*u % n, scale, evaluate)))
              for u in products]
    half = 1
    while half < len(values):
        for start in range(0, len(values), 2*half):
            for i in range(start, start+half):
                left, right = values[i], values[i+half]
                values[i], values[i+half] = (left+right)/2, (left-right)/2
        half *= 2
    return values


def try_multicomponent_split(n, evaluate, draw=randbelow, *, magnitude_only=False):
    """One trial on a composite with at least two distinct odd prime components.

    Complex mode takes absolute error <=1/256. Magnitude mode takes positive
    estimates at relative error <=1/4096, up to a common gain at each scale.
    """
    if n < 9 or n % 2 == 0:
        raise ValueError("Expected an odd composite")
    bits = n.bit_length()
    cutoff = max(2*bits, 4096 if magnitude_only else 0)
    for p in primerange(3, cutoff+1):
        if n % p == 0:
            return p
    bound = bits // ((2*bits).bit_length()-1)
    products = [1]
    while len(products) < 8*(bound+1)**2:
        unit = 1+draw(n-1)
        divisor = gcd(unit, n)
        if divisor > 1:
            return divisor
        products += [u*unit % n for u in products]

    with mp.workprec(3*bits+64):
        spectra = [mask_spectrum(n, products, 1, scale, evaluate) for scale in (1, 2)]
        channels = []
        for index in range(1, len(products)):
            scale = 1 if spectra[0][index] >= spectra[1][index] else 2
            if spectra[scale-1][index] > mp.mpf(1)/8:
                channels.append((index, scale))
        if not 2 <= len(channels) <= bound:
            return None

        # Each row carries the actual unit realizing its vector of Legendre signs.
        basis = {}
        for _ in range(len(channels)+4):
            unit = 1+draw(n-1)
            divisor = gcd(unit, n)
            if divisor > 1:
                return divisor
            shifted = [mask_spectrum(n, products, unit, scale, evaluate) for scale in (1, 2)]
            signs = sum(1 << i for i,(index, scale) in enumerate(channels)
                        if shifted[scale-1][index] < 0)
            while signs:
                pivot = signs.bit_length()-1
                if pivot not in basis:
                    basis[pivot] = (signs, unit)
                    break
                row, row_unit = basis[pivot]
                signs ^= row
                unit = unit*row_unit % n
            if 0 in basis:
                break
        if 0 not in basis:
            return None
        mask = basis[0][1]

        if magnitude_only:
            def cosine_at(scale):
                positive = aligned_projection(n, 1, scale, evaluate)
                masked = aligned_projection(n, mask, scale, evaluate)
                return ((positive/masked)**2-5)/4

            return divisor_from_cosines(n, cosine_at)

        def point_at(scale):
            positive = aligned_projection(n, 1, scale, evaluate)
            masked = aligned_projection(n, mask, scale, evaluate)
            difference = positive - masked
            if abs(difference) < 1:  # May occur on a signature collision.
                return 0j
            return 2*masked/difference

        return divisor_from_phases(n, point_at)
