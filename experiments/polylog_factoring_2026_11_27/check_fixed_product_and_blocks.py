"""Exact high-word histograms on dyadic inverse graphs, including a prescribed low residue."""

from cmath import exp
from collections import Counter
from fractions import Fraction
from math import gcd, tau


def block_histogram(n, k, b, low_bits=0, residue=1):
    modulus, bins = 1 << k, 1 << b
    assert n & 1 and residue & 1 and b+max(1, low_bits) <= k
    step = 1 << max(1, low_bits)
    counts = Counter()
    for x in range(residue % step, modulus, step):
        y = n*pow(x, -1, modulus) % modulus
        counts[x >> (k-b), y >> (k-b)] += 1
    size = modulus//step
    assert sum(counts.values()) == size
    for u in range(bins):
        assert sum(counts[u, v] for v in range(bins)) == size//bins
        assert sum(counts[v, u] for v in range(bins)) == size//bins
    tv = sum(abs(counts[u, v]*bins*bins-size) for u in range(bins) for v in range(bins))
    return Fraction(tv, 2*size*bins*bins)


def check_characters():
    probes = 0
    for k in range(3, 8):
        modulus = 1 << k
        phases = [exp(1j*tau*j/modulus) for j in range(modulus)]
        for b in range(1, k):
            bins = 1 << b
            length = modulus//bins
            for r in range(1, bins):
                coeffs = [sum(phases[(r*(x//length)*length-j*x) % modulus]
                              for x in range(modulus))/modulus for j in range(modulus)]
                assert all(abs(value) < 1e-10 for j, value in enumerate(coeffs) if (j-r) % bins)
                assert sum(map(abs, coeffs)) <= 4+2*k
                probes += 1
        inv = [(x, pow(x, -1, modulus)) for x in range(1, modulus, 2)]
        for a in range(1, min(21, modulus)):
            for c in range(1, min(21, modulus)):
                value = sum(phases[(a*x+c*y) % modulus] for x, y in inv)/len(inv)
                divisor = gcd(a, c, modulus)
                assert abs(value)**2 <= 128*divisor/modulus+1e-10
    return probes


def fiber_census(n, k, h, pattern):
    modulus = 1 << k
    mask = sum(pattern[j % len(pattern)] << j for j in range(k))
    cells = Counter()
    for x in range(1, modulus, 2):
        y = n*pow(x, -1, modulus) % modulus
        cells[x & y, (h*x) & (h*y) & mask] += 1
    sums = 0
    for size in sorted(cells.values()):
        sums += size
        if 2*sums >= modulus//2:
            return size, max(cells.values())


if __name__ == "__main__":
    print("Literal coarse character probes:", check_characters())
    for k, b in ((12, 2), (16, 3), (20, 4)):
        values = [block_histogram(n, k, b, ell, a)
                  for n in (1, -1, 82791) for ell, a in ((0, 1), (b, 1), (b, 3))]
        print("k,b,max TV unconditioned, max TV with prescribed low word:",
              k, b, float(max(values[::3])), float(max(values[1::3]+values[2::3])))
    for n in (1, -1, 82791):
        print("n, 16-bit median and maximum fiber for the two sparse schedules:",
              n, fiber_census(n, 16, 3, (0, 1, 1)), fiber_census(n, 16, 5, (0, 1, 0, 1, 1)))
