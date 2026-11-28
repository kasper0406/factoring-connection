"""Check fixed first-vector orbits, public character menus, and the binary norm search."""

from collections import Counter
from itertools import combinations, product
from math import gcd, isqrt, prod

from sympy import jacobi_symbol, primerange


# Columns give (a,b,c,d) with ad-bc=1. The public frame is (a,b,0), N(c,d,0).
ROWS = [(0, 1, -1, 0), (1, 1, -1, 0), (2, 1, -1, 0), (6, 1, -1, 0),
        (2, 3, -1, -1), (3, 7, -1, -2), (3, 1, -1, 0)]


def slope_partners(n, row):
    a, b, c, d = row
    norm = a*a + b*b
    result = {(n*c, n*d, 0)}
    for ell in range(1, isqrt(2*n // norm) + 1):
        bound = isqrt(2*n - norm*ell*ell)
        for k in range(-bound, bound + 1):
            if gcd(k, ell) != 1:
                continue
            denominator = k*k + norm*ell*ell
            scale, rem = divmod(2*n, denominator)
            if rem == 0:
                t = scale*ell*ell
                v = (n*c + b*t, n*d - a*t, scale*k*ell)
                assert gcd(*v) == n // gcd(n, denominator)
                result.add(v)
    return result


def direct_partners(n, row):
    a, b, c, d = row
    norm = a*a + b*b
    result = set()
    for t in range(2*n // norm + 1):
        squared = t*(2*n-norm*t)
        z = isqrt(squared)
        if z*z == squared:
            result.update((n*c+b*t, n*d-a*t, signed) for signed in {z, -z})
    return result


def orbit(row):
    a, b, c, d = row
    assert a*d - b*c == 1
    norm = a*a + b*b
    interval = range(-isqrt(norm), isqrt(norm) + 1)
    orbit_size = 0
    for u in product(interval, repeat=3):
        if sum(t*t for t in u) == norm:
            assert sorted(map(abs, u)) == sorted((a, b, 0))
            orbit_size += 1
    return orbit_size


def check_residue_signatures():
    # The CRT slots give n mod 8 and the Legendre symbols of n at odd primes.
    # All 5,13,29,37 are 1 mod 4, so quadratic reciprocity has positive sign.
    def signatures(extended):
        primes = [5, 37, 13, 29] if extended else [5, 37]
        local = [[(1 if pow(a, (r-1)//2, r) == 1 else -1) for a in range(1, r)]
                 for r in primes]
        hist = Counter()
        for residue in (1, 3, 5, 7):
            minus = 1 if residue % 4 == 1 else -1
            two = 1 if residue in (1, 7) else -1
            for values in product(*local):
                signs = (minus, minus*two, minus*values[0], minus*values[1])
                if extended:
                    signs += (minus*values[2], minus*two*values[3])
                hist[signs] += 1
        assert len(hist) == (64 if extended else 16) and len(set(hist.values())) == 1
        return len(hist), sum(hist.values()), 8*prod(primes)
    return signatures(False), signatures(True)


def check_frames():
    primes = list(primerange(7, 48))
    cases, partners = 0, 0
    for row in ROWS:
        a, b, c, d = row
        norm = a*a + b*b
        orbit_size = orbit(row)
        for p, q in list(combinations(primes, 2)) + [(101, 109)]:
            n = p*q
            if jacobi_symbol(-norm, n) != -1:
                continue
            direct = direct_partners(n, row)
            assert direct == slope_partners(n, row)
            r = p if jacobi_symbol(-norm, p) == 1 else q
            contents = Counter(gcd(*v) for v in direct)
            assert contents == {n: 24 // orbit_size, n // r: 48 // orbit_size}
            for v in direct:
                assert a*v[0]+b*v[1] == n*(a*c+b*d)
                assert sum(x*x for x in v) == n*n*(c*c+d*d)
            cases += 1
            partners += len(direct)
    assert [jacobi_symbol(-a*a-b*b, 11009) for a,b,_,_ in ROWS[:4]] == [1, 1, 1, -1]
    assert 8**2 + 37 == 101 and 11009 % 101 == 0
    return cases, partners


if __name__ == "__main__":
    print("Signatures, reduced residue classes, moduli:", check_residue_signatures())
    print("Mixed Grams and fixed-first-vector partners:", check_frames())
