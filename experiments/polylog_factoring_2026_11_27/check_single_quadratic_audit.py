"""Independent checks of the ten-coordinate quadratic count encoding."""

from math import gcd
from random import Random

from check_single_quadratic_count import direct_contrast, stationary_correction


def encoded_count(n, c0, j, width, offset, alpha, beta, polarity):
    """Enumerate the free coordinates and check every auxiliary constraint."""
    modulus = n * n
    half = (n - 1) // 2
    c1 = c0 + n * j
    delta = c1 * c1 - c0 * c0
    total = 0
    max_k = max_q = 0
    hyperbola_checks = 0
    for lower in range(width):
        low_value = lower - offset
        for y in range(n):
            x = low_value + n * y
            for t0 in range(n):
                for z in range(half + 1):
                    t = t0 + n * z
                    for selector in (0, 1):
                        signed = alpha * low_value + beta
                        k = signed // n
                        residue = signed - n * k
                        eta = int(residue > half)
                        assert (half + 1) * eta <= residue
                        assert residue <= half + (n - 1 - half) * eta
                        xor = selector ^ eta
                        assert xor >= selector - eta
                        assert xor >= eta - selector
                        assert xor <= selector + eta
                        assert xor <= 2 - selector - eta
                        complement = 1 - xor if polarity == 1 else xor
                        constant = c0 * c0 + selector * delta
                        quotient, residual = divmod(x * t - x * x - constant, modulus)
                        if residual == 0:
                            assert (x * (t - x) - constant) % modulus == 0
                        assert abs(k) <= n + 1
                        assert abs(quotient) <= 8 * modulus
                        max_k = max(max_k, abs(k))
                        max_q = max(max_q, abs(quotient))
                        if complement == 0:
                            accepted = residual == 0
                        else:
                            accepted = 1 <= residual <= modulus - 1
                        total += accepted
                        if accepted and complement == 0:
                            assert x * (t - x) == constant + modulus * quotient
                            hyperbola_checks += 1
    return total, max_k, max_q, hyperbola_checks


def audit():
    rng = Random(271127)
    packed = hyperbolas = 0
    max_k_ratio = max_q_ratio = 0.0
    for n in (9, 15, 21, 25, 45):
        for _ in range(3):
            c0 = next(x for x in range(rng.randrange(1, n), n) if gcd(x, n) == 1)
            j, offset, alpha, beta = (rng.randrange(n) for _ in range(4))
            width = rng.randrange(1, min(5, n + 1))
            baseline = width * n * n * ((n + 1) // 2)
            args = n, c0, j, width, offset, alpha, beta
            contrast = direct_contrast(*args)
            correction = stationary_correction(*args)
            for polarity in (-1, 1):
                value, max_k, max_q, checked = encoded_count(*args, polarity)
                assert value == baseline + polarity * contrast
                centered = value - baseline - polarity * correction
                assert centered == polarity * (contrast - correction)
                max_k_ratio = max(max_k_ratio, max_k / (n + 1))
                max_q_ratio = max(max_q_ratio, max_q / (n * n))
                packed += 1
                hyperbolas += checked
    print(
        f"{packed} encodings, {hyperbolas} accepted hyperbola points, "
        f"max |k|/(N+1)={max_k_ratio:.3f}, max |q|/N^2={max_q_ratio:.3f}"
    )


if __name__ == "__main__":
    audit()
