"""Check Fourier sampling of the normalized nonunit indicator of a semiprime."""

from cmath import exp
from fractions import Fraction
from math import gcd, isqrt, tau
from random import Random

from sympy import nextprime


def weight(p, q, k):
    n = p*q
    return 1+p*(p-2)*(k % p == 0)+q*(q-2)*(k % q == 0)+2*n*(k == 0)


def interval_mass(p, q, lo, hi):
    """Exact mass on [lo,hi) using factor-supplied progression counts."""
    n = p*q
    assert 0 <= lo <= hi <= n
    multiples = lambda s: (hi-1)//s-(lo-1)//s
    value = hi-lo+p*(p-2)*multiples(p)+q*(q-2)*multiples(q)+2*n*multiples(n)
    return Fraction(value, n*(p+q-1))


def exact_metrics(p, q):
    n, support = p*q, p+q-1
    phi = (p-1)*(q-1)
    split = Fraction(phi*(support-1), n*support)
    categories = ((phi, 1), (q-1, (p-1)**2), (p-1, (q-1)**2), (1, support**2))
    tv = Fraction(sum(count*abs(w-support) for count, w in categories), 2*n*support)
    assert sum(count*w for count, w in categories) == n*support
    assert split == 1-Fraction(support**2+phi, n*support)
    return Fraction(support-1, n), tv, split


def factor_from_prefix_estimates(n, estimate):
    """Find the smaller balanced prime using estimates of the mass on [1,t].

    The callable estimate(n,t) returns a rational with absolute error at most
    1/(40*sqrt(n)). The decoder compares against 1/(16*sqrt(n)) exactly.
    """
    if n < 15 or n % 2 == 0:
        raise ValueError("expected an odd balanced semiprime")
    upper = isqrt(n)
    if n < 1024:
        return next(x for x in range(3, upper+1, 2) if n % x == 0)
    lo, hi = 1, upper
    while lo < hi:
        mid = (lo+hi)//2
        score = estimate(n, mid)
        if score > 0 and 256*n*score*score > 1:
            hi = mid
        else:
            lo = mid+1
    if n % lo or lo <= 1:
        raise ValueError("inconsistent prefix estimates or balanced semiprime promise")
    return lo


def check_bisection(p, q):
    n = p*q
    shifts = Fraction(1, 40*(isqrt(n-1)+1))
    calls = 0
    def least_favorable(modulus, t):
        nonlocal calls
        calls += 1
        assert modulus == n and 1 <= t <= isqrt(n)
        exact = interval_mass(p, q, 1, t+1)
        assert exact == Fraction(t+p*(p-2)*(t >= p), n*(p+q-1))
        return exact+shifts if t < p else exact-shifts
    assert factor_from_prefix_estimates(n, least_favorable) == p
    assert calls <= isqrt(n).bit_length()
    return calls


def check_direct(p, q):
    n, support = p*q, p+q-1
    nonunits = [x for x in range(n) if gcd(x, n) > 1]
    assert len(nonunits) == support
    weights = []
    for k in range(n):
        transform = sum(exp(-1j*tau*((k*x) % n)/n) for x in nonunits)
        assert abs(abs(transform)**2-weight(p, q, k)) < 1e-7
        weights.append(weight(p, q, k))
    discrepancy, tv, split = exact_metrics(p, q)
    prefixes = [0]
    for w in weights:
        prefixes.append(prefixes[-1]+w-support)
    denominator = n*support
    assert prefixes[-1] == 0
    assert Fraction(max(prefixes)-min(prefixes), denominator) == discrepancy
    assert Fraction(sum(abs(w-support) for w in weights), 2*denominator) == tv
    assert Fraction(sum(w for k, w in enumerate(weights) if 1 < gcd(k, n) < n), denominator) == split


def check_large(p, q, rng):
    n = p*q
    discrepancy, tv, split = exact_metrics(p, q)
    assert tv == split
    assert interval_mass(p, q, 0, 1)-Fraction(1, n) == discrepancy
    for _ in range(500):
        lo, hi = sorted((rng.randrange(n+1), rng.randrange(n+1)))
        assert abs(interval_mass(p, q, lo, hi)-Fraction(hi-lo, n)) <= discrepancy
        pivot = rng.randrange(lo, hi+1)
        assert interval_mass(p, q, lo, pivot)+interval_mass(p, q, pivot, hi) == interval_mass(p, q, lo, hi)
    return float(discrepancy), float(tv), float(split)


if __name__ == "__main__":
    for primes in ((3, 5), (5, 7), (7, 11), (11, 13), (13, 17)):
        check_direct(*primes)
    for p, q in ((29, 37), (31, 53), (37, 43)):
        check_direct(p, q)
        check_bisection(p, q)
    rng = Random(410322)
    larger = sorted(int(nextprime((1 << 29)+rng.randrange(1 << 29))) for _ in range(2))
    for p, q in ((5, 7), (11, 13), (42697, 70019), larger):
        print(p, q, p*q, *check_large(p, q, rng), "bisection calls", check_bisection(p, q))
