"""Factor an odd integer with two distinct prime divisors via centered disk residues."""

import json
from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt


@dataclass(frozen=True)
class Request:
    form: tuple[int, int, int]
    cutoff: int
    root: int


def chart(n, m, r, baseline=False):
    """Send a primitive form and a unit discriminant root to the disk oracle."""
    B = (n - 1) // 2
    d = n * (B + 1) - 1 + (n * n if baseline else 0)
    a, b, h = m * m, m * (2 * r - n), r * r - n * r + d
    g = gcd(a, b, h)
    # m is a power of 2 and 2r-n is odd, so the content divides m.
    assert m % g == 0
    return Request((a // g, b // g, h // g), B // g, (2 * (m // g)) % n)


def hidden_parity(n, m, r, oracle):
    main = oracle(n, chart(n, m, r))
    axes = oracle(n, chart(n, m, r, True))
    public = int((n - 1) % m == r)
    value = (main - axes - 2 * public) % 4
    if value not in (0, 2):
        raise ArithmeticError("odd signed oracle count")
    return value // 2


def split(n, oracle):
    """The oracle takes (N, Request) and returns the signed nonzero count modulo 4."""
    if n < 3 or n % 2 == 0:
        raise ValueError(
            "expected an odd integer with at least two distinct prime divisors"
        )
    m, r, queries = 2, 0, 0
    while m < n:
        m *= 2
        queries += 2
        if not hidden_parity(n, m, r, oracle):
            r += m // 2
    divisor = gcd(r - 1, n)
    if not (1 < divisor < n and n % divisor == 0 and pow(r, 2, n) == 1):
        raise ArithmeticError("no verified split; check the promise or the oracle")
    return {
        "N": n,
        "root": r,
        "factors": sorted((divisor, n // divisor)),
        "queries": queries,
    }


def recover_prime_power(n):
    """For n promised to be a prime power, its largest exact exponent gives the prime."""
    bits = n.bit_length()
    for exponent in range(bits - 1, 1, -1):
        lo, hi = 2, 1 << ((bits + exponent - 1) // exponent)
        while lo <= hi:
            base = (lo + hi) // 2
            value = base**exponent
            if value == n:
                return base, exponent
            if value < n:
                lo = base + 1
            else:
                hi = base - 1
    return n, 1


def factor(n, oracle):
    """Fully factor any positive integer using the same centered count oracle."""
    if n < 1:
        raise ValueError("expected a positive integer")
    result = Counter()
    while n % 2 == 0:
        n //= 2
        result[2] += 1
    pending = [n] if n > 1 else []
    while pending:
        node = pending.pop()
        if hidden_parity(node, 2, 0, oracle):
            pending.extend(split(node, oracle)["factors"])
        else:
            p, exponent = recover_prime_power(node)
            result[p] += exponent
    return dict(sorted(result.items()))


def rows(n, req):
    """Exact row intervals of the supplied form by completing its square."""
    a, b, d = req.form
    height = n * req.cutoff
    delta = 4 * a * d - b * b
    vmax = isqrt(4 * a * height // delta)
    for v in range(-vmax, vmax + 1):
        w = isqrt(4 * a * height - delta * v * v)
        yield v, -((b * v + w) // (2 * a)), (w - b * v) // (2 * a)


def enumerated_oracle(n, req, work=None):
    """Reference evaluator from public data. The first queries can scan linearly many x."""
    a, b, d = req.form
    total = 0
    for v, lo, hi in rows(n, req):
        if work is not None:
            work["positions"] += max(0, hi - lo + 1)
        total += sum(
            (u or v) != 0 and (a * u * u + b * u * v + d * v * v) % n == 0
            for u in range(lo, hi + 1)
        )
    return total % 4


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument(
        "--complete",
        action="store_true",
        help="also handle powers and recursively split",
    )
    args = parser.parse_args()
    work = Counter()

    def oracle(n, q):
        work["queries"] += 1
        return enumerated_oracle(n, q, work)

    result = (
        {"N": args.N, "prime_factors": factor(args.N, oracle)}
        if args.complete
        else split(args.N, oracle)
    )
    print(json.dumps({**result, "work": dict(work)}, indent=2))
