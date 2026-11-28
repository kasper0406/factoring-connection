"""Check literal binomials, all prefix widths on small products, and larger oracle controls."""

from itertools import combinations
from math import comb, isqrt
from random import Random
from statistics import median

from sympy import isprime, primerange

if __package__:
    from .fixed_binomial_streams import factor_by_fixed_streams, fixed_menu
else:
    from fixed_binomial_streams import factor_by_fixed_streams, fixed_menu


SETTINGS = [(3, (1,)), (5, (0, 1, 0, 1, 1))]


def check_literal():
    cases = calls = 0
    for p, q in combinations(primerange(3, 137), 2):
        if q >= 2*p:
            continue
        n = p*q
        for h, schedule in SETTINGS:
            menu = fixed_menu(n, h, schedule)
            seen = []

            def literal(modulus, t):
                seen.append(t)
                return comb(2*t, t) % modulus == 0

            stats = {}
            assert factor_by_fixed_streams(n, literal, multiplier=h, schedule=schedule,
                                          stats=stats) == p
            assert seen == [t for _, _, t in menu]
            for j, width in enumerate(stats["widths"], 1):
                modulus = 1 << j
                masks = [(m, 1 << level) for level, m, _ in menu if level < j]
                direct = 0
                for x in range(1, modulus, 2):
                    y = n * pow(x, -1, modulus) % modulus
                    if all((m*x & m*y & bit) == (m*p & m*q & bit) for m, bit in masks):
                        direct += 1
                assert width == direct
            cases += 1
            calls += len(seen)
    return cases, calls


def divisible(t, prime):
    # Legendre's formula for the factorial valuations.
    a, b = 2*t, t
    while a:
        a //= prime
        b //= prime
        if a > 2*b:
            return True
    return False


def dyadic_prime(bits, rng):
    while True:
        candidate = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if isprime(candidate):
            return candidate


def check_larger():
    rng = Random(93157)
    for bits in (24, 48, 96, 128):
        samples = []
        while len(samples) < 16:
            p, q = sorted(dyadic_prime(bits, rng) for _ in range(2))
            if p != q:
                samples.append((p, q))
        for h, schedule in SETTINGS:
            peaks, trials, calls = [], [], []
            for p, q in samples:
                n = p*q
                stats = {}
                def valuation(modulus, t):
                    assert modulus == n
                    return divisible(t, p) and divisible(t, q)
                assert factor_by_fixed_streams(n, valuation, multiplier=h, schedule=schedule,
                                               stats=stats) == p
                lower = isqrt((n - 1)//2) + 1
                assert isqrt(n) - lower < 1 << (lower.bit_length()-1)
                peaks.append(max(stats["widths"]))
                trials.append(stats["trials"])
                calls.append(stats["queries"])
            print(bits, h, (min(calls), max(calls)), (median(peaks), max(peaks)), max(trials))


if __name__ == "__main__":
    print("literal input/menu pairs and calls", *check_literal())
    print("factor bits, h, calls (min,max), peak width (median,max), max final divisions")
    check_larger()
