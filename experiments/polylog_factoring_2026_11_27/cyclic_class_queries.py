"""Encode the fixed-discriminant queries as powers of one explicit class."""

from math import gcd, isqrt

from fixed_discriminant_split import query_form


def compose(left, right, diagonal, scale):
    return (left + right) * pow(1 - diagonal * left * right, -1, scale) % scale


def class_exponent(target, diagonal, scale):
    current, generator_power, exponent = 0, 1, 0
    for bit in range(scale.bit_length() - 1):
        if (current ^ target) & (1 << bit):
            current = compose(current, generator_power, diagonal, scale)
            exponent |= 1 << bit
        generator_power = compose(generator_power, generator_power, diagonal, scale)
    assert current == target
    return exponent


def parameter(exponent, diagonal, scale):
    real, imag = 1, 0
    base_real, base_imag = 1, 1
    while exponent:
        if exponent & 1:
            real, imag = (
                (real * base_real - diagonal * imag * base_imag) % scale,
                (real * base_imag + imag * base_real) % scale,
            )
        base_real, base_imag = (
            (base_real**2 - diagonal * base_imag**2) % scale,
            2 * base_real * base_imag % scale,
        )
        exponent //= 2
    return imag * pow(real, -1, scale) % scale


def encode(n, start):
    (a, b, c), cutoff = query_form(n, start)
    scale = isqrt(a)
    shear = b // (2 * scale)
    diagonal = c - shear**2
    exponent = class_exponent(pow(shear, -1, scale), diagonal, scale)
    return (scale**2, 2 * scale, diagonal + 1), exponent, cutoff


def check_product(left, right, diagonal, scale):
    """Compare product lattices by their generators and integer index."""
    expected = compose(left, right, diagonal, scale)
    vectors = [
        (1 - diagonal * left * right, left + right),
        (-diagonal * left * scale, scale),
        (-diagonal * right * scale, scale),
        (-diagonal * scale**2, 0),
    ]
    assert all((y - expected * x) % scale == 0 for x, y in vectors)
    index = 0
    for i, (x, y) in enumerate(vectors):
        for z, w in vectors[:i]:
            index = gcd(index, x*w-y*z)
    assert index == scale


if __name__ == "__main__":
    products = 0
    for scale in (2, 4, 8, 16, 32):
        for diagonal in (2, 4, 10, 22):
            for left in range(scale):
                for right in range(scale):
                    check_product(left, right, diagonal, scale)
                    products += 1
    queries = 0
    for n in range(5, 321, 2):
        for start in range((n - 1) // 4 + 1):
            generator, exponent, cutoff = encode(n, start)
            leading, cross, constant = generator
            scale = isqrt(leading)
            assert cross == 2 * scale
            diagonal = constant - 1
            (a, b, c), reference_cutoff = query_form(n, start)
            assert cutoff == reference_cutoff
            shear = b // (2 * isqrt(a))
            assert c - shear**2 == diagonal
            assert parameter(exponent, diagonal, scale) * shear % scale == 1
            assert exponent % 2 == 1
            queries += 1
    print(f"{products} lattice products and {queries} encoded disk queries")
