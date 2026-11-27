"""Recover roots in short public intervals via a five-dimensional lattice."""

from math import gcd, isqrt
from random import Random

from sympy import Matrix, Poly, Symbol, integer_nthroot

from fixed_window_split import reference_count
from nonadaptive_quarter import queries


def interval_limit(n):
    return max(1, int(integer_nthroot(n*n, 5)[0]) // 4)


def quadratic_roots(n, linear, constant, limit):
    """Find all 0 <= x <= limit with x^2 + linear*x + constant = 0 mod n."""
    if limit <= 1:
        return [x for x in range(limit+1) if (x*x+linear*x+constant) % n == 0]
    if limit > interval_limit(n):
        raise ValueError("Interval exceeds the proven lattice bound")
    linear %= n
    constant %= n
    rows = [
        [n*n, 0, 0, 0, 0],
        [0, n*n, 0, 0, 0],
        [n*constant, n*linear, n, 0, 0],
        [0, n*constant, n*linear, n, 0],
        [constant**2, 2*constant*linear, linear**2+2*constant, 2*linear, 1],
    ]
    scales = [limit**i for i in range(5)]
    lattice = Matrix([[value*scale for value, scale in zip(row, scales)] for row in rows])
    short = list(lattice.lll()[0, :])
    assert sum(abs(entry) for entry in short) < n*n
    x = Symbol("x")
    auxiliary = Poly(sum((entry//scale)*x**i for i, (entry, scale) in
                         enumerate(zip(short, scales))), x)
    return sorted(
        int(root) for root in auxiliary.ground_roots()
        if root.is_Integer and 0 <= root <= limit and
        (root*root+linear*root+constant) % n == 0
    )


def roots_from_quotient_interval(n, y, low, high):
    """Find roots in [ceil(n/2), n) whose square is y+k*n, low <= k <= high."""
    if not (0 <= y < n and 0 <= low <= high and high-low <= interval_limit(n)):
        raise ValueError("Expected a canonical target and a short quotient interval")
    start = max((n+1)//2, isqrt(max(0, y+low*n-1))+1)
    end = min(n-1, isqrt(y+high*n))
    if end < start:
        return []
    offsets = quadratic_roots(n, 2*start, start**2-y, end-start)
    return [start+t for t in offsets]


def progression_roots(n, origin, step, limit):
    """Recover roots of 1 represented by origin+step*t, 0 <= t <= limit."""
    inverse = pow(step, -1, n)
    roots = quadratic_roots(n, 2*origin*inverse, (origin**2-1)*inverse**2, limit)
    return [(origin+step*t) % n for t in roots]


def split_shortened_quarter(n, disk_count):
    cutoff, reference, forms = queries(n)
    baseline = disk_count(n, reference, cutoff)
    remaining = interval_limit(n) + 1
    bits = max(1, ((n + remaining - 1)//remaining - 1).bit_length())
    bit, prefix = 0, 0
    for scale, form in forms[:bits-1]:
        residue = (disk_count(n, form, cutoff) - baseline) % 8
        assert residue in (0, 4)
        known = int(min(scale, n-scale) <= n//4)
        bit ^= 1 ^ (residue//4) ^ known
        prefix = 2*prefix + bit
    denominator = 1 << bits
    low = (prefix*n + denominator - 1) // denominator
    high = ((prefix+1)*n - 1) // denominator
    for root in progression_roots(n, low, 1, high-low):
        divisor = gcd(root-1, n)
        if 1 < divisor < n:
            return divisor, bits-1
    raise AssertionError("Short interval contained no mixed root")


if __name__ == "__main__":
    rng = Random(261127)
    checked = 0
    for p, q in [(101, 109), (503, 937), (1009, 1877), (541, 1201), (101**2, 109), (10007, 10009)]:
        n = p*q
        limit = interval_limit(n)
        actual = 1+p*(-2*pow(p, -1, q) % q)
        for _ in range(50):
            step = rng.randrange(1, n)
            while gcd(step, n) != 1:
                step = rng.randrange(1, n)
            origin = (actual - step*rng.randrange(limit+1)) % n
            observed = progression_roots(n, origin, step, limit)
            expected = [(origin+step*t) % n for t in range(limit+1)
                        if ((origin+step*t)**2-1) % n == 0]
            assert observed == expected
            checked += 1
    for p, q in [(31, 37), (101, 109), (503, 937)]:
        divisor, count = split_shortened_quarter(p*q, reference_count)
        assert divisor in (p, q)
        print(f"N={p*q}, {count} scaled disk queries, divisor {divisor}")
    print(f"{checked} bounded lattice checks")
