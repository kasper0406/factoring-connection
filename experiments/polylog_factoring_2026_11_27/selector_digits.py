"""Synthetic reciprocal-digit checks for bounded idempotent box counts.

The parity version keeps one selector in a fixed lower half. The ternary
version keeps both mixed selectors and recovers the sum of their next
binary digits. Literal counts are checked on small grids. Larger fixtures
provide the factors only to synthesize the oracle answers.
"""

from math import gcd, isqrt, pi

from stationary_log import spectrum_at


def bounded_count(modulus, bound, half):
    """Count the pairs literally, omitting the two public selectors."""
    stop = (modulus + 1) // 2 if half else modulus
    return sum(
        selector * value % modulus == 0
        for selector in range(2, stop)
        if selector * (selector - 1) % modulus == 0
        for value in range(bound + 1)
    )


def scaled_candidates(modulus, numerator, denominator, uncertainty):
    """Integers in [N*numerator/denominator, N*(numerator+uncertainty)/denominator)."""
    low = (modulus * numerator + denominator - 1) // denominator
    stop = (modulus * (numerator + uncertainty) + denominator - 1) // denominator
    return range(low, stop)


def half_parity_candidates(modulus, parity):
    """Recover the selected annihilator size to a window of two integers."""
    quotient = 0
    scale = 1
    while 2 * scale <= modulus - 1:
        scale *= 2
        quotient = 2 * quotient + ((parity(scale) - 1) % 2)
    candidates = scaled_candidates(modulus, quotient, scale, 1)
    assert len(candidates) <= 2
    return candidates


def full_ternary_candidates(modulus, residue):
    """Recover the component sum to a window of four integers."""
    quotient_sum = 0
    scale = 1
    while 2 * scale <= modulus - 1:
        scale *= 2
        next_digit_sum = (residue(scale) - 2 - 2 * quotient_sum) % 3
        quotient_sum = 2 * quotient_sum + next_digit_sum
    candidates = scaled_candidates(modulus, quotient_sum, scale, 2)
    assert len(candidates) <= 4
    return candidates


def split_from_total(modulus, total):
    discriminant = total**2 - 4 * modulus
    if discriminant < 0:
        return None
    difference = isqrt(discriminant)
    if difference**2 != discriminant or (total - difference) % 2:
        return None
    low = (total - difference) // 2
    high = (total + difference) // 2
    if low > 1 and low * high == modulus and gcd(low, high) == 1:
        return low, high
    return None


def check_full_parity_spectrum(p, q):
    modulus = p * q
    signal = [(-1)**((bound // p + bound // q) % 2) for bound in range(modulus)]
    positions = [(p + q) // 2, abs(p - q) // 2]
    peaks = spectrum_at(signal, positions)
    assert min(peaks) >= (2 / pi)**2 - 1e-12
    for j in range(1, modulus.bit_length()):
        bound = modulus >> j
        assert int(signal[bound] < 0) == ((p ^ q) >> j) & 1
    return min(peaks)


def check_fixture(p, q):
    modulus = p * q
    first = p * pow(p, -1, q)
    selected = min(first, modulus + 1 - first)
    kept = gcd(selected, modulus)
    assert kept in (p, q)

    half_calls, full_calls = [], []

    def parity(bound):
        half_calls.append(bound)
        count = 1 + bound // (modulus // kept)
        if modulus <= 1_225:
            assert bounded_count(modulus, bound, half=True) == count
        return count % 2

    def ternary(bound):
        full_calls.append(bound)
        count = 2 + bound // p + bound // q
        if modulus <= 1_225:
            assert bounded_count(modulus, bound, half=False) == count
        return count % 3

    half = half_parity_candidates(modulus, parity)
    full = full_ternary_candidates(modulus, ternary)
    assert kept in half and p + q in full
    assert [d for d in half if 1 < d < modulus and modulus % d == 0] == [kept]
    assert [parts for total in full if (parts := split_from_total(modulus, total))] == [
        tuple(sorted((p, q)))
    ]
    assert half_calls == full_calls
    assert half_calls[-1] <= modulus - 1 < 2 * half_calls[-1]
    print(
        f"N={modulus}: {len(half_calls)} queries per decoder, maximum bound "
        f"{half_calls[-1]}, parity candidates {len(half)}, ternary candidates {len(full)}"
    )


if __name__ == "__main__":
    peaks = []
    for components in (
        (3, 5), (9, 5), (9, 7), (9, 25), (25, 49), (101, 109), (101**2, 109),
        (153247, 254119), (1000000007, 1000000009), (3**90, 5**46),
    ):
        check_fixture(*components)
        if components[0] * components[1] <= 1225:
            peaks.append(check_full_parity_spectrum(*components))
    print(f'{len(peaks)} unrestricted parity spectra checked; smallest predicted peak {min(peaks):.6f}')
