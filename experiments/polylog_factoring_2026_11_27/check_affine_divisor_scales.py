"""Check all public centers on small moduli and certificates at larger scales."""

from math import gcd
from random import Random

if __package__:
    from .affine_divisor_scales import divisor_scales, try_affine_scales
    from .affine_prime_window import rigorous_ok
else:
    from affine_divisor_scales import divisor_scales, try_affine_scales
    from affine_prime_window import rigorous_ok


def run():
    centers = 0
    for n in (4, 8, 9, 15, 25, 27, 30, 35, 49, 60, 77, 91, 121, 143):
        bands = list(divisor_scales(n))
        b = max(t for t in range(2,n) if n % t == 0)
        lower, radius = next((lower, radius) for lower, radius in bands if lower <= b < 2*lower)
        hits = 0
        guaranteed = 0
        for x in range(n):
            expected = {g for lim,h in bands for y in range(-h,h+1)
                        if lim <= (g := gcd(x+y,n)) < n and g > 1}
            result = try_affine_scales(n,x)
            assert result in expected if expected else result is None
            hits += result is not None
            guaranteed += any(b <= gcd(x+y,n) < n and (x+y) % b == 0
                              for y in range(-radius,radius+1))
        assert guaranteed == (2*radius+1)*(n//b-1) and hits >= guaranteed
        centers += n

    rng = Random(13)
    certificates = 0
    for bits in (12, 24, 36, 54, 80):
        for n in (1 << bits, (1 << bits)-rng.randrange(1,100)):
            for lower, radius in divisor_scales(n):
                assert 2*radius < lower
                if radius:
                    assert rigorous_ok(n, lower, radius, (n-1).bit_length(), lower.bit_length()-1)
                    certificates += 1

    # Prime powers and multiple prime factors are valid hidden divisors for
    # the affine lattice. These cases use a nonzero error in the chosen band.
    for n, b in ((3*1009, 1009), (7**4 * 11, 7**4), (101*103*107, 101*103)):
        radius = next(h for lower,h in divisor_scales(n) if lower <= b < 2*lower)
        assert radius >= 1
        f = try_affine_scales(n, b + radius)
        assert f is not None and 1 < f < n and n % f == 0
    print(f"{centers} complete center checks, {certificates} dyadic lattice certificates, "
          "three planted unbalanced or repeated-factor cases")


if __name__ == "__main__":
    run()
