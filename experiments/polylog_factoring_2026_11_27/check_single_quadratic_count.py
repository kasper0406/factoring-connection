"""Check the rank-two count against four independently counted slices."""

from math import gcd
from random import Random


def packed_count(n, c, j, width, offset, alpha, beta, polarity):
    """Enumerate the residual-quotient formulation over its free coordinates."""
    modulus = n * n
    half = (n - 1) // 2
    delta = (c + n * j) ** 2 - c**2
    count = 0
    for s in (0, 1):
        constant = c**2 + s * delta
        for lower in range(width):
            low_value = lower - offset
            eta = int((alpha * low_value + beta) % n > half)
            exclusive = s ^ eta
            complement = 1 - exclusive if polarity == 1 else exclusive
            low = complement
            high = (modulus - 1) * complement
            for y in range(n):
                x = low_value + n * y
                base = -x * x - constant
                for t in range(n * (half + 1)):
                    quotient, residual = divmod(x * t + base, modulus)
                    assert abs(quotient) <= 8 * modulus
                    count += low <= residual <= high
    return count


def direct_contrast(n, c, j, width, offset, alpha, beta):
    """Use direct inversion for the four positive congruence counts."""
    modulus = n * n
    half = (n - 1) // 2
    by_digit = [0, 0]
    for x in range(modulus):
        x0 = x % n
        if (x0 + offset) % n >= width or gcd(x, n) != 1:
            continue
        sign = 1 if (alpha * x0 + beta) % n <= half else -1
        inverse = pow(x, -1, modulus)
        for s, scale in enumerate((c, c + n * j)):
            phase = (x + scale * scale * inverse) % modulus
            if phase // n <= half:
                by_digit[s] += sign
    return by_digit[1] - by_digit[0]


def stationary_correction(n, c, j, width, offset, alpha, beta):
    """Count the two public stationary blocks exactly."""
    modulus = n * n
    half = (n - 1) // 2
    by_digit = [0, 0]
    for root_sign in (-1, 1):
        lower = root_sign * c % n
        if (lower + offset) % n >= width:
            continue
        sign = 1 if (alpha * lower + beta) % n <= half else -1
        for s, scale in enumerate((c, c + n * j)):
            phase = 2 * root_sign * scale % modulus
            if phase // n <= half:
                by_digit[s] += n * sign
    return by_digit[1] - by_digit[0]


def check():
    rng = Random(261136)
    tested = 0
    for n in (15, 21, 45):
        half = (n - 1) // 2
        for _ in range(5):
            c = next(value for value in range(rng.randrange(1, n), n) if gcd(value, n) == 1)
            j, offset, alpha, beta = [rng.randrange(n) for _ in range(4)]
            width = rng.randrange(1, 5)
            baseline = width * n * n * (half + 1)
            arguments = n, c, j, width, offset, alpha, beta
            contrast = direct_contrast(*arguments)
            correction = stationary_correction(*arguments)
            for polarity in (-1, 1):
                packed = packed_count(*arguments, polarity)
                assert packed == baseline + polarity * contrast
                centered = packed - baseline - polarity * correction
                assert centered == polarity * (contrast - correction)
                assert (4 * centered <= n) == (4 * polarity * (contrast - correction) <= n)
                assert (4 * centered >= 3 * n) == (4 * polarity * (contrast - correction) >= 3 * n)
                tested += 1
    print(f"{tested} packed counts and threshold comparisons verified")


if __name__ == "__main__":
    check()
