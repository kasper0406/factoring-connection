"""Check the line permutation, Fourier discrepancy ingredients, and actual lifted witnesses."""

from math import gcd
from random import Random

from sympy import nextprime

if __package__:
    from .cofactor_tangent_lifts import inverse_line, quotient_witnesses
else:
    from cofactor_tangent_lifts import inverse_line, quotient_witnesses


def verify_permutation_and_windows():
    controls = 0
    for bits in range(3, 10):
        period = 1 << bits
        for n in range(1, min(period, 41), 2):
            actual = {((3*x+2*n*pow(x, -1, period)) % period):
                      (x, n*pow(x, -1, period) % period)
                      for x in range(1, period, 2)}
            assert len(actual) == period//2
            for line, pair in actual.items():
                assert inverse_line(n, line, period) == pair
                controls += 1
            # Exact second moment of hits in all shifted windows; independently
            # check its autocorrelation expansion on half masks of both factors.
            half = period//2
            yes = [int(bool(line in actual and actual[line][0] & half and actual[line][1] & half))
                   for line in range(period)]
            total = sum(yes)
            length = period//8
            cur = sum(yes[:length])
            squares = 0
            for line in range(period):
                squares += cur*cur
                cur += yes[(line+length) % period]-yes[line]
            other = length*total
            for gap in range(1, length):
                pairs = sum(a*yes[(j+gap) % period] for j, a in enumerate(yes))
                other += 2*(length-gap)*pairs
            assert squares == other
    return controls


def verify_witness(n, mask, wheel, answer, point):
    u, v, line, remainder = point
    assert u % 2 and v % 2 and u < v < 2*u and v == n//u
    assert u*v+remainder == n and 0 <= remainder < u
    assert 3*u+2*v == line
    assert remainder % ((1 << mask.bit_length())*wheel) == 0
    assert gcd(u*v, wheel) == 1 and (u-1) % wheel == (v-n) % wheel == 0
    assert bool(u & mask and v & mask) == answer


def verify_large():
    rng = Random(47708)
    summary = {}
    for side_bits in (32, 48, 64, 96, 128):
        counts = {label: [0, 0, 0, 0] for label in ("below", "quarter", "wheel")}
        for _ in range(16):
            while True:
                p, q = [int(nextprime(2**(side_bits-1)+rng.randrange(2**(side_bits-1))))
                        for _ in range(2)]
                if p != q:
                    break
            n = p*q
            quarter = n.bit_length()//4
            for label, j, wheel in (("below", max(5, quarter-5), 1),
                                    ("quarter", quarter, 1),
                                    ("wheel", max(6, quarter-12), 105)):
                for mask in (1 << j, (1 << j) | (1 << (j//2))):
                    found, cost = quotient_witnesses(n, mask, 512, wheel=wheel)
                    stats = counts[label]
                    stats[0] += 1
                    stats[1] += len(found) == 2
                    if len(found) == 2:
                        stats[2] = max(stats[2], cost)
                    for answer, point in found.items():
                        verify_witness(n, mask, wheel, answer, point)
                        stats[3] += point.remainder == 0
        summary[side_bits] = counts
    return summary


if __name__ == "__main__":
    print("Direct inverse-line comparisons", verify_permutation_and_windows())
    print("By bits per prime: queries, both outcomes found, maximum successful lines, exact products", verify_large())
