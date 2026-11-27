"""Check half-period witnesses in the odd class powers."""

from math import isqrt

from cyclic_class_queries import class_exponent
from fixed_discriminant_split import query_form
from fixed_window_split import mixed_root, reference_count


def mixed_disk_parity(n, scale, diagonal, cutoff, shear):
    form = (scale**2, 2*scale*shear, diagonal + shear**2)
    radius = isqrt(n*cutoff - diagonal)
    first = -((radius + shear) // scale)
    last = (radius - shear) // scale
    offset = shear * pow(scale, -1, n) % n
    trivial = sum((x + offset) % n in (1, n - 1) for x in range(first, last + 1))
    axis = 2 * sum(x*x % n == 0 for x in range(1, isqrt(n*cutoff) // scale + 1))
    return ((reference_count(n, form, cutoff) - axis - 2*trivial) % 4) // 2


def check(n):
    (leading, cross, last), cutoff = query_form(n, 0)
    scale = isqrt(leading)
    diagonal = last - (cross // (2*scale))**2
    half = (n - 1) // 2
    k = half // 2
    omitted = (half,) if half % 2 == 0 else (half - 1, half, 2*half)
    for start in range(n):
        expected = sum(mixed_root(n, start + i) for i in omitted) % 2
        if not expected:
            continue
        lower = scale // 8
        shear = lower + (scale*(start + k) - lower) % n
        if shear % 2 == 0:
            shear += n
        shifted = shear - scale // 2
        assert lower <= shear < scale//4
        exponents = [
            class_exponent(pow(value, -1, scale), diagonal, scale)
            for value in (shear, shifted)
        ]
        assert (exponents[1] - exponents[0]) % scale == scale//2
        first = mixed_disk_parity(n, scale, diagonal, cutoff, shear)
        second = mixed_disk_parity(n, scale, diagonal, cutoff, shifted)
        assert first ^ second == expected
        return True
    assert all(not mixed_root(n, value) for value in range(n))
    return False


if __name__ == "__main__":
    witnesses = sum(check(n) for n in range(5, 321, 2))
    assert witnesses == 84
    print(f"{witnesses} half-period witnesses among 158 odd moduli")
