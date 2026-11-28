"""Compare rational representations for a square t and the square z=t/(1+t).

Draw a unit s and let x=(s−1/s)/2, y=(s+1/s)/2 modulo odd N. Keep x,y
units, so y²=x²+1 and both t=x² and z=(x/y)² are units. Use only t,z in
the two representation searches. If t=u/v, then 1+t=(u+v)/v and
z=u/(u+v). Require both starting fractions |uv|, |v(u+v)|≤N and measure
the target height |u(u+v)|, which may exceed N. Negative numerators and
denominators are allowed. Count targets with a smooth squarefree kernel
and separately compare against every primitive short representation of z.

For each bit length generate four existing balanced moduli and use 1,024
accepted s per modulus. All tables count every candidate without early
stopping. 'extra' counts targets hit by the paired construction but missed
by the reference shell on z. Run directly or as a module, optionally
passing --targets or a list of bit lengths.
"""

import argparse
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, prod
from random import Random
from time import perf_counter

import gmpy2

if __package__:
    from . import check_rational_height_shell as shell
else:
    import check_rational_height_shell as shell


COLUMNS = (
    "targets", "rejected", "t_rows", "paired_low_rows", "paired_high_rows",
    "z_rows", "low_hits", "high_hits", "low_and_high", "paired_hits",
    "z_hits", "extra",
)


def make_kernel_test(primes):
    primorial = gmpy2.mpz(prod(primes))

    @lru_cache(maxsize=120_000)
    def accepted(value):
        common = int(gmpy2.gcd(primorial, value))
        while common > 1:
            value //= common
            common = gcd(common, value)
        return isqrt(value) ** 2 == value

    return accepted


def unit_fractions(n, target):
    candidates, _ = shell.representations(n, target)
    for numerator, denominator in candidates:
        if gcd(denominator, n) == 1:
            assert denominator > 0 and gcd(numerator, denominator) == 1
            assert (numerator - target * denominator) % n == 0
            yield numerator, denominator


def compare(n, t, z, accepted, counts):
    z_shell = set(unit_fractions(n, z))
    counts["z_rows"] += len(z_shell)
    reference = False
    for u, w in z_shell:
        reference |= accepted(abs(u)) and accepted(w)

    low = high = False
    for u, v in unit_fractions(n, t):
        counts["t_rows"] += 1
        w = u + v
        if abs(v * w) > n:
            continue
        assert w and gcd(w, n) == 1 and gcd(u, w) == 1
        assert (u - z * w) % n == 0
        if abs(u * w) <= n:
            counts["paired_low_rows"] += 1
            assert (u if w > 0 else -u, abs(w)) in z_shell
            low |= accepted(abs(u)) and accepted(abs(w))
        else:
            counts["paired_high_rows"] += 1
            high |= accepted(abs(u)) and accepted(abs(w))
    assert not low or reference
    counts["low_hits"] += low
    counts["high_hits"] += high
    counts["low_and_high"] += low and high
    counts["paired_hits"] += low or high
    counts["z_hits"] += reference
    counts["extra"] += high and not reference


def measure(bits, seed, pool, accepted):
    n = shell.balanced_semiprime(bits, seed).N
    rng = Random(f"consecutive-square-pool:{bits}:{seed}")
    inv2 = (n + 1) // 2
    counts = Counter()
    while counts["targets"] < pool:
        s = rng.randrange(1, n)
        if gcd(s, n) != 1:
            counts["rejected"] += 1
            continue
        inverse = pow(s, -1, n)
        x, y = (s - inverse) * inv2 % n, (s + inverse) * inv2 % n
        if gcd(x * y, n) != 1:
            counts["rejected"] += 1
            continue
        t = x * x % n
        z = t * pow(t + 1, -1, n) % n
        assert (y * y - t - 1) % n == 0 and pow(x * pow(y, -1, n), 2, n) == z
        counts["targets"] += 1
        compare(n, t, z, accepted, counts)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bits", nargs="*", type=int, default=[48, 64, 80, 96])
    parser.add_argument("--targets", type=int, default=1024)
    args = parser.parse_args()
    if min(args.bits) < 16 or args.targets < 1:
        parser.error("use at least 16 bits and a positive target count")
    print("bits", "seed", *COLUMNS, "seconds", sep=",", flush=True)
    for bits in args.bits:
        accepted = make_kernel_test(shell.primes_upto(bits ** 2))
        total = Counter()
        start = perf_counter()
        for seed in range(4):
            before = perf_counter()
            counts = measure(bits, seed, args.targets, accepted)
            total.update(counts)
            print(bits, seed, *(counts[k] for k in COLUMNS), f"{perf_counter()-before:.3f}",
                  sep=",", flush=True)
        print(bits, "sum", *(total[k] for k in COLUMNS), f"{perf_counter()-start:.3f}",
              sep=",", flush=True)


if __name__ == "__main__":
    main()
