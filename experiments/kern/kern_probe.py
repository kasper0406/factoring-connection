"""Reproduce the preregistered tests with the repository's virtualenv, using -B."""
import json
import math
import random
from collections import Counter
from itertools import combinations
from pathlib import Path

import gmpy2
import numpy as np
from cypari2 import Pari
from flint import fmpz_mod_poly_ctx
from coef_theta_check import norm_representations


def coords(poly):
    return tuple(int(poly[i]) for i in range(3))


def cubic_norm(values):
    a, b, c = values
    e = a + c
    f = b - c
    return a * (e * e - b * f) + c * (b * e - c * f) - b * (b * b - c * e)


def probe(n, coins):
    ring = fmpz_mod_poly_ctx(n)
    x = ring.gen()
    cubic = x ** 3 - x + 1
    readouts = np.zeros(8, dtype=int)
    hits = np.zeros(8, dtype=int)
    factors = set()
    for triple in coins:
        a = ring(triple)
        norm = cubic_norm(triple) % n
        d = math.gcd(norm, n)
        if d != 1:
            if d != n:
                factors.add(d)
            continue
        b = (pow(a, 3, cubic) * pow(norm, -1, n)) % cubic
        assert cubic_norm(coords(b)) % n == 1
        for offset, seed in ((0, a), (4, b)):
            s = seed
            for j in range(4):
                s = pow(s, n, cubic)
                diff = coords(s - seed)
                norm_diff = cubic_norm(diff) % n
                hits[offset + j] += int(all(t == 0 for t in diff))
                for d in (math.gcd(norm_diff, n), math.gcd(n, *diff)):
                    if 1 < d < n:
                        factors.add(d)
                readouts[offset + j] += int(gmpy2.jacobi(norm_diff, n))
    assert all(n % d == 0 and 1 < d < n for d in factors)
    return readouts, hits, bool(factors)


def prime(bits, rng, seen, pari):
    while True:
        x = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if x not in seen and gmpy2.is_prime(x) and pari.isprime(x):
            return x


def trace_check(pari):
    gram = pari('[3,0,2;0,2,-3;2,-3,2]')
    assert int(pari.matdet(gram)) == -23
    diagonal = [pari(3), pari(2), pari(-23) / 6]
    primes = [r for r in range(5, 500) if r != 23 and pari.isprime(r)]
    data = {}
    for r in primes:
        roots = sum((x ** 3 - x + 1) % r == 0 for x in range(r))
        hasse = math.prod(int(pari.hilbert(a, b, r)) for a, b in combinations(diagonal, 2))
        assert hasse == 1
        data[r] = (roots, int(gmpy2.jacobi(-23, r)))
    cells = Counter()
    for p, q in combinations(primes, 2):
        n = p * q
        r_p, s_p = data[p]
        r_q, s_q = data[q]
        bit = int(norm_representations(n) == 4)
        assert bit == int(r_p == 0 and r_q == 0)
        signs = [int(pari.hilbert(n, -23, r)) for r in (p, q)]
        assert signs == [s_p, s_q]
        assert math.prod(signs) == int(gmpy2.jacobi(-23, n))
        cells[','.join(map(str, sorted(signs) + [bit]))] += 1
    return dict(cells)


def correlation_test(rows):
    keys = np.asarray([24 * n + mod for n, mod, _, _ in rows])
    cells = [np.where(keys == key)[0] for key in np.unique(keys)]
    a = np.asarray([v for _, _, v, _ in rows], dtype=float)
    b = np.asarray([bit for _, _, _, bit in rows], dtype=float)
    for g in cells:
        a[g] -= a[g].mean(axis=0)
        b[g] -= b[g].mean()
    a /= np.linalg.norm(a, axis=0)
    b /= np.linalg.norm(b)
    observed = float(abs(a.T @ b).max())
    rng = np.random.default_rng(683131)
    exceed = 0
    for _ in range(1999):
        perm = b.copy()
        for g in cells:
            perm[g] = rng.permutation(b[g])
        exceed += abs(a.T @ perm).max() >= observed
    return dict(rows=len(rows), max_correlation=observed, permutation_p=(1 + exceed) / 2000)


def main():
    pari = Pari()
    curve = pari.ellinit([0, 0, 0, -1, 1])
    f = pari('x^3-x+1')
    prime_rng = random.Random(112739)
    read_rng = random.Random(971303)
    rows = []
    result = {}
    for width in (40, 60, 80):
        used = set()
        hits = np.zeros(8, dtype=int)
        splits = nontrivial = odd = 0
        while len(used) < 2400:
            p, q = (prime(width // 2, prime_rng, used, pari) for _ in range(2))
            n = p * q
            if p == q or n.bit_length() != width:
                continue
            used.update((p, q))
            coins = [[read_rng.randrange(n) for _ in range(3)] for _ in range(4)]
            vals, returns, split = probe(n, coins)
            splits += split
            hits += returns
            bit = int((p + 1 - int(pari.ellap(curve, p))) *
                      (q + 1 - int(pari.ellap(curve, q))) % 2)
            odd += bit
            degrees = [[int(pari.poldegree(g)) for g in pari.factormod(f, r)[0]]
                       for r in (p, q)]
            assert bit == int(degrees == [[3], [3]])
            kernel = math.prod(math.gcd(n, r ** k - 1) for r, deg in zip((p, q), degrees)
                               for k in deg)
            nontrivial += kernel > 1
            if not split and gmpy2.jacobi(-23, n) == 1:
                rows.append((width, n % 24, vals, bit))
        result[width] = dict(splits=splits, hits=hits.tolist(), kernel_gt1=nontrivial, odd=odd)
    output = dict(widths=result, signal=correlation_test(rows), trace=trace_check(pari))
    msg = json.dumps(output, indent=2, default=float)
    Path(__file__).with_name('kern_results.json').write_text(msg + '\n')
    print(msg)


if __name__ == '__main__':
    main()
