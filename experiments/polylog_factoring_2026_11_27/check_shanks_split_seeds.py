"""Check input-generated cyclic cubic splitting certificates and their seed counts."""

from collections import Counter
from itertools import combinations
from math import gcd
from sympy import isprime, primerange


def cubic(x, t):
    return x**3 - t*x*x - (t+3)*x - 1


def discriminant_root(t):
    return t*t + 3*t + 9


def split_seed(n, x):
    """Return a proper divisor if exposed, otherwise a good centered t or None."""
    if n <= 3 or n % 2 == 0:
        raise ValueError("expected an odd input greater than three")
    x %= n
    for v in (x, x+1, x*x+x+1, 3):
        common = gcd(n, v)
        if common != 1:
            return (common if common < n else None), None
    denominator = x*(x+1) % n
    t = ((x**3-3*x-1)*pow(denominator, -1, n) + n//2) % n - n//2
    return None, t


def roots_from_seed(n, x):
    b = (-(x+1)*pow(x, -1, n)) % n
    return (x % n, b, (-(b+1)*pow(b, -1, n)) % n)


def check_local_primes():
    cases = 0
    for t in range(-38, 39):
        conductor = discriminant_root(t)
        if not isprime(conductor):
            continue
        for s in primerange(5, 180):
            if s == conductor:
                continue
            roots = [x for x in range(s) if cubic(x, t) % s == 0]
            expected = 3 if pow(s, (conductor-1)//3, conductor) == 1 else 0
            assert len(roots) == expected
            if roots:
                assert set(roots_from_seed(s, roots[0])) == set(roots)
            cases += 1
    return cases


def check_fibers():
    cases = 0
    for p, q in combinations(primerange(5, 35), 2):
        n = p*q
        good = Counter()
        for x in range(n):
            divisor, t = split_seed(n, x)
            if divisor:
                assert divisor in (p, q)
            if t is not None:
                deg_three = x*x+x+1
                denominator = x*(x+1)
                assert (discriminant_root(t)*denominator**2 - deg_three**3) % n == 0
                assert gcd(discriminant_root(t), n) == 1
                for z in roots_from_seed(n, x):
                    assert cubic(z, t) % n == 0
                good[t] += 1
        assert set(good.values()) == {9}
        allowed = [(s-4 if s%3 == 1 else s-2)//3 for s in (p, q)]
        assert len(good) == allowed[0]*allowed[1]
        cases += 1
    return cases


def check_input():
    n = 713
    for seed, expected, conductor in ((7,-7,37), (19,-14,163)):
        divisor, t = split_seed(n, seed)
        assert divisor is None and t == expected and discriminant_root(t) == conductor
        assert all(cubic(z,t) % n == 0 for z in roots_from_seed(n, seed))
        assert isprime(conductor)
        residues = [r for r in range(1, conductor) if pow(r, (conductor-1)//3, conductor) == 1]
        assert len(residues) == (conductor-1)//3
        assert all(p % conductor in residues for p in (23, 31))
        print(f"N={n}, seed={seed}, t={t}, prime conductor={conductor}, "
              f"{len(residues)} permitted unit residues")
    # The input and seeds are the only inputs to split_seed.
    assert 23*31 == n and pow(31, 12, 37) == pow(23, 12, 37) == 1
    assert {pow(a,3,23) for a in range(23)} == set(range(23))


if __name__ == "__main__":
    print("Local character/root checks", check_local_primes())
    print("Exact composite nine-seed fiber checks", check_fibers())
    check_input()
