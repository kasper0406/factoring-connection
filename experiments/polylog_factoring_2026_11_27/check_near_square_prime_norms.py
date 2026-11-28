"""Check band norm neighbors and single-seed decoding of general prime-edge cycles."""

from collections import defaultdict
from itertools import combinations, combinations_with_replacement
from math import gcd, isqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.genusleak.forms import small_primes, sqrt_mod_prime


def band_neighbor(n, a, t):
    """Return the unique possible partner and coordinate for the supplied root."""
    if a % 2 == 0 or gcd(n, a) != 1 or not n <= a*a <= 4*n or (t*t+n) % a:
        raise ValueError("expected an odd coprime band prime power and a square root of -n")
    x = t % a
    if x % 2:
        x = a-x
    b = (x*x+n)//a
    return (b, x) if n <= b*b <= 4*n else None


def seed_sign(n, ell, t):
    """Given an odd split prime and a root, return its seed sign modulo ±1 or None."""
    if ell < 3 or ell % 2 == 0 or gcd(n, ell) != 1 or (t*t+n) % ell:
        raise ValueError("expected a coprime odd prime and a square root of -n")
    t %= ell
    lift = t + ell*(-((t*t+n)//ell)*pow(2*t, -1, ell) % ell)
    a, b = (ell*ell, 0), (lift, 1)

    def norm(v):
        return v[0]*v[0] + n*v[1]*v[1]

    # Rank-two Gauss reduction for the exact weighted norm.
    while True:
        if norm(b) < norm(a):
            a, b = b, a
        aa = norm(a)
        ab = a[0]*b[0] + n*a[1]*b[1]
        k = (2*ab+aa)//(2*aa)
        b = (b[0]-k*a[0], b[1]-k*a[1])
        if norm(b) >= aa:
            break
    if aa != ell*ell:
        return None
    sign = a[0]*pow(ell, -1, n) % n
    assert sign*sign % n == 1
    return min(sign, n-sign)


def check_band():
    pairs = [(p, q) for p, q in combinations(small_primes(199), 2) if 7 <= p and q <= 2*p]
    vertices = edges = loops = prime_edges = prime_loops = 0
    for p, q in pairs:
        n = p*q
        powers = {}
        for ell in small_primes(isqrt(4*n)):
            a = ell
            while a*a <= 4*n:
                if n < a*a and a % 2 and gcd(a, n) == 1:
                    powers[a] = ell
                a *= ell
        vertices += len(powers)
        partners = {}
        for a, b in combinations_with_replacement(sorted(powers), 2):
            x = isqrt(a*b-n)
            if x*x != a*b-n:
                continue
            assert a not in partners and b not in partners
            partners[a] = b, x
            partners[b] = a, x
            edges += 1
            loops += a == b
            prime_edges += a == powers[a] and b == powers[b]
            prime_loops += a == b == powers[a]
            s = isqrt(a*b)
            if s*s == a*b:
                assert s == (p+q)//2 and {s-x, s+x} == {p, q}
        for a in powers:
            roots = [t for t in range(a) if (t*t+n) % a == 0]
            if not roots:
                assert a not in partners
            for t in roots:
                neighbor = band_neighbor(n, a, t)
                if neighbor and neighbor[0] in powers:
                    assert neighbor == partners[a]
                else:
                    assert a not in partners
    print(f"{len(pairs)} semiprimes, {vertices} prime-power vertices, {edges} edges, "
          f"{loops} loops; prime-only {prime_edges} edges, {prime_loops} loops")
    # A square product can also occur on a nonloop between even-power vertices.
    n = 2823013
    assert band_neighbor(n, 41**2, 534) == (43**2, 534)
    assert (41*43-534, 41*43+534) == (1229, 2297)
    assert (41*43)**2 - 534**2 == n


def check_cycles():
    seeds = cycles = mixed = empty = global_signs = 0
    for n in (35, 55, 77, 91, 143, 187, 221, 247):
        primes = [v for v in small_primes(151) if v % 2 and gcd(v, n) == 1
                  and pow(-n, (v-1)//2, v) == 1]
        graph = defaultdict(list)
        signs = {}
        for ell in primes:
            t = sqrt_mod_prime(-n, ell)
            sign = seed_sign(n, ell, t)
            assert sign == seed_sign(n, ell, -t)
            brute = set()
            for y in range(1, isqrt(ell*ell//n)+1):
                x = isqrt(ell*ell-n*y*y)
                if x*x+n*y*y == ell*ell:
                    z = x*pow(ell, -1, n) % n
                    brute.add(min(z, n-z))
            assert brute == (set() if sign is None else {sign})
            signs[ell] = sign
            seeds += 1
            empty += sign is None
            global_signs += sign == 1
        for a, b in combinations_with_replacement(primes, 2):
            for y in range(1, isqrt(a*b//n)+1):
                x = isqrt(a*b-n*y*y)
                if x*x+n*y*y == a*b:
                    graph[a].append((b, x))
                    if a != b:
                        graph[b].append((a, x))
        potential = {}
        for seed in primes:
            if seed in potential:
                continue
            potential[seed] = 1
            stack = [seed]
            while stack:
                a = stack.pop()
                assert signs[a] == signs[seed]
                for b, x in graph[a]:
                    new = potential[a]*x*pow(a, -1, n) % n
                    if b not in potential:
                        potential[b] = new
                        stack.append(b)
                        continue
                    z = new*pow(potential[b], -1, n) % n
                    assert z*z % n == 1 and min(z, n-z) in (1, signs[seed])
                    cycles += 1
                    if min(z, n-z) != 1:
                        assert 1 < gcd(z-1, n) < n
                        mixed += 1
    assert empty and global_signs and mixed
    print(f"{seeds} unrestricted split-prime seeds, {empty} nonprincipal squares, "
          f"{global_signs} global seed signs, {cycles} closure checks, {mixed} mixed")


if __name__ == '__main__':
    check_band()
    check_cycles()
