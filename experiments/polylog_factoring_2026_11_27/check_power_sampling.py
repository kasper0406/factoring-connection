"""Compare the two decoders on uniform multiplicative power cosets."""

from fractions import Fraction
from math import gcd, isqrt, lcm


def square_splits(z, n):
    residue = z * z % n
    t = isqrt(residue)
    return t * t == residue and 1 < gcd(z - t, n) < n


def main():
    comparisons = 0
    maximum = Fraction(0)
    for p, q in [(3, 5), (5, 7), (7, 11), (11, 13), (13, 19), (11, 37)]:
        n = p * q
        units = [u for u in range(1, n) if gcd(u, n) == 1]
        for exponent in range(lcm(p - 1, q - 1)):
            image = {pow(u, exponent, n) for u in units}
            a = len({u % p for u in image})
            b = len({u % q for u in image})
            assert len(image) == a * b
            direct = sum(1 < gcd(u - 1, n) < n for u in image)
            assert direct == a + b - 2
            for scale in units:
                count = sum(square_splits(scale * u % n, n) for u in image)
                if direct:
                    assert count <= 2 * direct
                    maximum = max(maximum, Fraction(count, direct))
                    comparisons += 1
                else:
                    assert count == square_splits(scale, n)
    print(f"{comparisons} nonconstant comparisons, largest ratio {maximum}")


if __name__ == "__main__":
    main()
