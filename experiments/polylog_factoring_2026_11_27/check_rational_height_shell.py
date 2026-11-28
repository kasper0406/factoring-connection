"""Exhaust the short rational representations of the matched random-square targets.

For a unit t modulo N consider nonzero integers u,v with u = tv mod N,
v a unit, v > 0 and |uv| <= N. Dividing by g=gcd(u,v) preserves this
congruence since g divides the unit v. It reduces the height and removes
g^2 from |uv|, preserving its squarefree kernel. If the original fraction
is smooth on both sides, its reduction is also smooth. Thus reduced
fractions suffice for both existence questions below.

Write x=t/N and a=(tv-u)/N. Reduction gives gcd(a,v)=1 and
|x-a/v| <= 1/v^2. Choose the last convergent p_i/q_i with q_i <= v,
so q_i <= v < q_(i+1). The terminal denominator N cannot occur since
v is a unit and v <= N. The determinant |a q_i - v p_i| is at most
q_i/v + v/q_(i+1) < 2. If 2q_i < v < q_(i+1)/2, this determinant
is below 1, forcing a/v=p_i/q_i and contradicting v>q_i.

If v <= 2q_i, write (a,v) in the unimodular basis of convergents i-1,i.
The preceding coefficient has absolute value at most 1 by the determinant
bound; the other has absolute value at most 3 by comparing denominators.
If v >= q_(i+1)/2, use the basis i,i+1 instead. The determinant with
convergent i+1 is less than 3 (at most 2 if that convergent is terminal),
so the i coefficient is at most 2 in absolute value; the other is at most 3.
Under (a,v) -> (tv-Na,v), these bases become consecutive extended-Euclid
rows. Their boundary rows are (+/-N,0) and (0,+/-N). Up to overall sign
it suffices to enumerate the coprime coefficients of magnitude at most 3.
There are O(log N) Euclidean steps, each taking a constant number of pairs.

Run directly with Python or as a module from the repository root.
"""

from functools import lru_cache
from math import gcd, isqrt, prod
from pathlib import Path
from random import Random
import sys

import gmpy2

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.instances import balanced_semiprime
from src.smoothbuild.core import primes_upto


COEFFICIENTS = [(1, 0)] + [
    (a, b) for a in range(-3, 4) for b in range(1, 4) if gcd(a, b) == 1
]


def representations(n, t):
    """Return all shell candidates and the ordinary convergents before unit screening."""
    previous_u, previous_v = n, 0
    u, v = t, 1
    candidates = set()
    convergents = set()
    while True:
        for a, b in COEFFICIENTS:
            num = a * previous_u + b * u
            den = a * previous_v + b * v
            if num and den and abs(num * den) <= n:
                candidates.add((num if den > 0 else -num, abs(den)))
        if u == 0:
            return candidates, convergents
        convergents.add((u if v > 0 else -u, abs(v)))
        quotient, remainder = divmod(previous_u, u)
        previous_u, u = u, remainder
        previous_v, v = v, previous_v - quotient * v


def check_small_moduli():
    checked = 0
    for n in range(3, 151):
        for t in range(1, n):
            if gcd(t, n) != 1:
                continue
            expected = {
                (num, den)
                for den in range(1, n + 1) if gcd(den, n) == 1
                for num in range(-(n // den), n // den + 1)
                if num and gcd(num, den) == 1 and (num - t * den) % n == 0
            }
            candidates, _ = representations(n, t)
            assert {(u, v) for u, v in candidates if gcd(v, n) == 1} == expected
            checked += 1
    return checked


def benchmark(bits):
    base_product = gmpy2.mpz(prod(primes_upto(bits * bits)))

    @lru_cache(maxsize=120_000)
    def classify(value):
        """Return 2 for smooth, 1 for smooth up to a square, and 0 otherwise."""
        common = int(gmpy2.gcd(base_product, value))
        while common > 1:
            value //= common
            common = gcd(common, value)
        if value == 1:
            return 2
        return int(isqrt(value) ** 2 == value)

    convergent_hits = shell_hits = kernel_hits = candidate_count = 0
    for seed in range(4):
        n = balanced_semiprime(bits, seed).N
        rng = Random(f"relation-baseline:{bits}:{seed}")
        for _ in range(1024):
            root = rng.randrange(1, n)
            while gcd(root, n) != 1:
                root = rng.randrange(1, n)
            t = root * root % n
            candidates, convergents = representations(n, t)
            candidate_count += len(candidates)
            best = best_convergent = 0
            for num, den in candidates:
                if gcd(den, n) != 1:
                    continue
                # Coprime lattice coefficients and a unit denominator imply
                # gcd(num, den)=1, so square cofactors can be tested separately.
                status = min(classify(abs(num)), classify(den))
                best = max(best, status)
                if (num, den) in convergents:
                    best_convergent = max(best_convergent, status)
            convergent_hits += best_convergent == 2
            shell_hits += best == 2
            kernel_hits += best >= 1
    return convergent_hits, shell_hits, kernel_hits, candidate_count


if __name__ == "__main__":
    print("Brute-force unit-target comparisons:", check_small_moduli())
    print("bits  convergent_hits  shell_hits  shell_hits_with_squares  candidates_before_gcd")
    for bits in (32, 48, 64, 80, 96):
        print(bits, *benchmark(bits), flush=True)
