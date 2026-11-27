"""Test one proposed sum residue on a balanced-semiprime sum interval."""

from math import gcd, isqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.nativemod.reductions import _rigorous_ok, small_root_divisor


def query(n, low, high, prime, residue):
    """Assume distinct odd prime factors and auxiliary prime ≡3 mod 4, coprime to n.

    Return a matching factor or None after complete coverage. Raise if too wide.
    Bounds must lie above the real discriminant threshold.
    """
    if low > high or low*low < 4*n or low < 0 or prime % 4 != 3 or gcd(prime,n) != 1:
        raise ValueError('Invalid sum bounds or auxiliary modulus')
    disc = (residue*residue-4*n) % prime
    delta = pow(disc, (prime+1)//4, prime)
    if delta*delta % prime != disc:
        return None
    classes = {(residue+sign)*pow(2,-1,prime) % prime for sign in (delta, -delta)}
    return query_classes(n, low, high, prime, classes, prime, residue)


def query_classes(n, low, high, stride, classes, sum_modulus, residue):
    """Search certified progressions for the larger factor, checking the sum."""
    def upper_factor(s):
        return (s + isqrt(s*s-4*n))//2
    first_q = max(isqrt(n)+1, upper_factor(low))
    last_q = upper_factor(high)
    inverse = pow(stride, -1, n)
    for cls in classes:
        first = first_q + (cls-first_q) % stride
        if first > last_q:
            continue
        span = (last_q-first)//stride
        center = first + stride*(span//2)
        radius = max(1, (span+1)//2)
        if not _rigorous_ok(n, isqrt(n)+1, radius, 9, 4):
            raise ValueError('Progression exceeds the certified lattice radius')
        factor = small_root_divisor(n, center*inverse % n, radius, isqrt(n)+1, 9, 4)
        if factor:
            total = factor + n//factor
            if low <= total <= high and total % sum_modulus == residue % sum_modulus:
                return factor
    return None


if __name__ == '__main__':
    from check_lehman_window_baseline import cohorts
    from sympy import nextprime
    prime = int(nextprime(2**34))
    while prime % 4 != 3:
        prime = int(nextprime(prime))
    p, q = next(cases[0] for bits, planted, cases in cohorts() if bits == 48 and not planted)
    n = p*q
    low, high = isqrt(4*n-1)+1, isqrt(9*n//2)
    residue = (p+q) % prime
    assert query(n,low,high,prime,residue) in (p,q)
    false = next(s for s in range(100) if s != residue and
                 pow((s*s-4*n)%prime,(prime-1)//2,prime)==1)
    assert query(n,low,high,prime,false) is None
    print(n, prime, residue, false, high-low+1)
