"""Exact controls for weighted_order_hint.py. Run directly with Python."""

from math import gcd, isqrt, prod
from pathlib import Path
import sys

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.weighted_order_hint import (
    _box, _certifies, _parameters, _weighted_basis, minimum_hint_bound, weighted_order_divisors,
)
from src.trackE.coppersmith import rigorous_ok


def check_basis():
    modulus, exponent = 23*31, 110
    for sign, r in [(1, 23), (-1, 31)]:
        g = gcd(exponent, r-sign)
        assert g > 1
        for center in (0, 27):
            radius, d = 40, 6
            t = r-center
            for m, n in [(0, 3), (1, 4), (2, 2), (4, 1), (3, 0)]:
                rows = _weighted_basis(modulus, exponent, sign, center, radius, d, m, n)
                determinant = (modulus**(m*(m+1)//2) * exponent**(n*(n+1)//2)
                               * radius**(d*(d-1)//2))
                assert prod(row[i] for i, row in enumerate(rows)) == determinant
                for i, row in enumerate(rows):
                    assert not any(row[i+1:])
                    assert all(v % radius**j == 0 for j, v in enumerate(row))
                    value = sum(v//radius**j * t**j for j, v in enumerate(row))
                    assert value % (r**m * g**n) == 0


def check_signs_and_bounds():
    for p, q, exponent, bound, sign, threshold in [
        (65_537, 98_317, 196_608, 65_536, 1, 42_336),
        (8_191, 12_289, 24_576, 8_192, -1, 7_484),
    ]:
        modulus = p*q
        lower, _, radius = _box(modulus, exponent, 0, None)
        assert gcd(exponent, p-sign) == bound < exponent
        assert minimum_hint_bound(modulus, exponent, 3) == threshold
        for trial_bound in (threshold-1, threshold, bound):
            brute = [(m, n) for m in range(3) for n in range(3)
                     if _certifies(modulus, exponent, lower, trial_bound, radius, 3, m, n)]
            assert bool(brute) == (trial_bound >= threshold)
        stats = {}
        assert weighted_order_divisors(modulus, exponent, threshold, sign=sign,
                                       max_dimension=3, stats=stats) == (p,)
        assert stats == dict(dimension=3, modulus_multiplicity=1, exponent_multiplicity=2)
        opposite = tuple(r for r in (p, q) if gcd(exponent, r+sign) >= bound)
        assert weighted_order_divisors(modulus, exponent, bound, sign=-sign,
                                       max_dimension=3) == opposite
        assert not any(rigorous_ok(modulus*exponent, lower*bound, radius, d, m)
                       for d in (2, 3) for m in range(1, d+1))
    assert weighted_order_divisors(65_537*114_689, 3*16_384, 16_384,
                                   max_dimension=5) == (65_537, 114_689)
    assert minimum_hint_bound(143, 2, 2) is None
    try:
        weighted_order_divisors(143, 2, 2, max_dimension=2)
    except ValueError as error:
        assert 'certificate' in str(error)
    else:
        raise AssertionError('expected an uncertified box')


def check_all_dimension_separation():
    p = 3_138_550_867_693_340_381_918_434_710_477_120_830_124_544_893_850_256_146_433
    q = 4_707_826_301_540_010_572_877_652_065_715_681_245_186_817_340_775_384_219_871
    modulus, exponent, bound = p*q, 3 << 113, 1 << 113
    lower, _, radius = _box(modulus, exponent, 0, None)
    assert gcd(exponent, p-1) == bound < exponent
    assert (lower*bound)**13 <= (modulus*exponent)**8
    assert radius**169 >= (modulus*exponent)**64
    # Put u=8/13. Even dropping the LLL/HG penalties, equality of the
    # multiplicities would require X^(d(d-1)) (Ne)^(m(m+1)) < (LG)^(2md).
    # The two integer inequalities above show that the exponent gap is at least
    # (m-du+1/2)^2 + du(1-u)-1/4 >= 0 for every integer d>=2 and any m.
    assert _parameters(modulus, exponent, lower, bound, radius, 7) == (7, 3, 6)
    assert weighted_order_divisors(modulus, exponent, bound, max_dimension=7) == (p,)
    print(f'{modulus.bit_length()}-bit uncentered example recovered in dimension 7; '
          'equal multiplicities fail the determinant bound in every dimension')


def check_centered():
    p, q, exponent, bound = 65_537, 66_541, 2**33, 2**16
    modulus = p*q
    center, radius = isqrt(modulus), 512
    assert exponent > modulus and gcd(exponent, p-1) == bound
    assert isqrt(isqrt(modulus)) < abs(p-center) <= radius
    assert weighted_order_divisors(modulus, exponent, bound, center=center, radius=radius,
                                   max_dimension=3) == (p,)
    assert weighted_order_divisors(modulus, exponent, bound, center=p, radius=0) == (p,)
    assert weighted_order_divisors(modulus, exponent, bound, center=-p, radius=0) == ()


if __name__ == '__main__':
    check_basis()
    check_signs_and_bounds()
    check_all_dimension_separation()
    check_centered()
    print('weighted order hint checks passed')
