"""Vectorized finite-size experiment for cyclic Fourier interval-energy retrieval."""
from math import gcd, isqrt
import numpy as np


def split_by_interval_energy(n, log_magnitude, rng, samples=32768):
    """Return (frequency, divisor) or None. Callback accepts arrays of unit scales.

    For this experiment missing scales receive the pilot mean. Only the final
    frequency gcd can return a factor. The sample budget is empirical.
    """
    modulus = n*n
    if n < 5 or modulus > isqrt(2**63-1) or samples < 1:
        raise ValueError("Expected a positive sample budget and an int64-safe modulus")
    pilot = rng.integers(1, modulus, 1024)
    units = pilot[np.gcd(pilot, n) == 1]
    if not len(units):
        return None
    center = np.mean(log_magnitude(units))

    def centered(scales):
        valid = np.gcd(scales, n) == 1
        result = np.zeros(len(scales))
        result[valid] = log_magnitude(scales[valid])-center
        return result

    return _frequency_search(n, modulus, centered, rng, samples)


def _frequency_search(n, modulus, centered, rng, samples):
    starts = np.array([1 << j for j in range((modulus//2).bit_length())])
    stops = np.minimum(2*starts, modulus//2+1)
    masses = np.r_[1., 2*(stops-starts)/starts]
    if n % 2 == 0:
        raise ValueError("Expected an odd base")
    intervals = [(1, modulus)]  # DC is irrelevant
    while intervals:
        shells = rng.choice(len(masses), samples, p=masses/masses.sum())
        indices = np.maximum(shells-1, 0)
        distance = rng.integers(starts[indices], stops[indices])
        distance[shells == 0] = 0
        y = (distance*rng.choice([-1, 1], samples)) % modulus
        t = rng.integers(modulus, size=samples)
        weighted = centered(t)*centered((t+y) % modulus)
        weighted *= masses.sum()*np.where(shells == 0, 1, starts[indices])/modulus
        nonzero = y != 0
        angle = 2j*np.pi/modulus
        base = np.expm1(angle*y[nonzero])
        def prefix(end):
            kernel = np.full(samples, float(end))
            kernel[nonzero] = (np.expm1(angle*((end*y[nonzero]) % modulus))/base).real
            return np.mean(weighted*kernel)

        children = []
        for lo, hi in intervals:
            mid = (lo+hi)//2
            middle = prefix(mid)
            for left, right, energy in [(lo, mid, middle-prefix(lo)),
                                        (mid, hi, prefix(hi)-middle)]:
                if energy < 0.025 or left == right:
                    continue
                if right-left == 1:
                    divisor = gcd(left, n)
                    if 1 < divisor < n:
                        return left, divisor
                else:
                    children.append((left, right))
        if len(children) > 50*(2*n.bit_length()+1)**2:
            return None
        intervals = children
    return None


def split_by_scale_comparisons(n, compare, rng, samples=32768):
    """Finite-size unit-only test: compare(left,right) returns a vector of booleans."""
    def score(scales):
        references = rng.integers(1, n*n, len(scales))
        invalid = np.gcd(references, n) != 1
        while np.any(invalid):
            references[invalid] = rng.integers(1, n*n, np.count_nonzero(invalid))
            invalid = np.gcd(references, n) != 1
        reverse = rng.integers(2, size=len(scales)).astype(bool)
        answer = compare(np.where(reverse, references, scales),
                         np.where(reverse, scales, references))
        return np.where(reverse, -1, 1)*np.where(answer, 1, -1)

    return split_by_interval_energy(n, score, rng, samples)


def split_by_doubling_bit(n, greater_at_t, rng, samples=32768):
    """greater_at_t(t) compares fixed-mask magnitudes at t and 2t modulo n²."""
    def score(scales):
        return np.where(greater_at_t(scales), 1, -1)

    return split_by_interval_energy(n, score, rng, samples)

def split_by_lift_comparison(n, compare, rng, samples=32768):
    """Compare mask-one magnitudes only at scales congruent to 1 modulo n."""
    if n < 5 or n % 2 == 0 or n > isqrt(2**63-1) or samples < 1:
        raise ValueError("Expected an odd base and int64-safe scale arithmetic")

    def score(k):
        left = 1+n*k
        right = 1+n*((2*k) % n)
        return np.where(compare(left, right), 1.0, -1.0)

    return _frequency_search(n, n, score, rng, samples)


def split_by_power_threshold(n, above_power, exponent, rng, samples=32768):
    """above_power(scales,j) compares a single mask-one magnitude to 2**j."""
    if n < 5 or n % 2 == 0 or n > isqrt(2**63-1) or samples < 1:
        raise ValueError("Expected an odd base and int64-safe scale arithmetic")
    if not 0 <= exponent <= n.bit_length():
        raise ValueError("Threshold exponent outside the search range")

    def score(k):
        return np.where(above_power(1+n*k, exponent), 1.0, -1.0)

    return _frequency_search(n, n, score, rng, samples)


def split_by_average_case_windows(n, length, window_mean, rng, samples=32768, padding=1):
    """Search on n using window means modulo n*padding; only the padding value is needed."""
    base = n*padding
    if n < 5 or padding < 1 or gcd(base, 6) != 1:
        raise ValueError("Expected positive target and padding, both coprime to six")
    if base > isqrt(2**63-1) or samples < 1:
        raise ValueError("Expected int64-safe scale arithmetic and a positive sample budget")
    if length < 1 or gcd(length, base) != 1:
        raise ValueError("Window length must be positive and invertible")
    inverse = pow(length, -1, base)
    coprime = gcd(padding, n) == 1
    if coprime:
        padding_inverse = pow(padding, -1, n)
    else:
        offset = int(rng.integers(padding))

    def score(k):
        if coprime:
            left = k if padding == 1 else k+n*rng.integers(padding, size=len(k))
            right = (left+padding*((k*padding_inverse) % n)) % base
        else:
            left = offset+padding*k
            right = offset+padding*((2*k) % n)
        anchor = rng.integers(base, size=len(k))
        def relative(start):
            step = ((anchor-start)*inverse) % base
            return window_mean(start, step)-window_mean((start+step) % base, step)
        return np.where(relative(left) > relative(right), 1., -1.)

    return _frequency_search(n, n, score, rng, samples)
