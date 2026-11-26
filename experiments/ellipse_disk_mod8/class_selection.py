"""Propose reduced prime-class minima from short rational approximations on ell=1."""

from collections import Counter
from math import gcd, isqrt

from class_prefix import branch_from_n, reduced_forms_at, represented_factors


def denominator_budget(n):
    """Integer floor of the sixteenth root."""
    for _ in range(4):
        n = isqrt(n)
    return n


def propose(br, limit, gap, work):
    """Return candidate A values and factors certified by the enumerated weights."""
    if br.ell != 1 or limit < 1 or gap < 1:
        raise ValueError("requires ell=1 and positive denominator and gap budgets")
    minima, direct = set(), set()
    for u in range(1, limit + 1):
        for v in range(u, isqrt(2 * u * u) + 2):
            work["weight_pairs"] += 1
            if gcd(u, v) != 1:
                continue
            work["primitive_weights"] += 1
            m = u * v
            radicand = 4 * br.n * m * m
            first = isqrt(radicand) + 1  # N is nonsquare
            offset = br.c * m
            for h in range(first, min(first + gap, offset + br.cap + 1)):
                work["gap_positions"] += 1
                minima.add(h - offset)
                disc = h * h - radicand
                root = isqrt(disc)
                if root * root != disc:
                    continue
                for num in (h - root, h + root):
                    p, rem = divmod(num, 2 * v * v)
                    if rem == 0 and 1 < p < br.n and br.n % p == 0:
                        q = br.n // p
                        if p <= q <= 2 * p:
                            direct.add(p)
    return minima, direct


def classes_at(delta, minima, work):
    """Enumerate primitive reduced forms up to inversion at the supplied minima."""
    result = set()
    for a in sorted(minima):
        work["minima_probed"] += 1
        result.update(reduced_forms_at(delta, a, work))
    work["classes"] += len(result)
    return sorted(result)


def search_minima(br, minima, build, search):
    classes = classes_at(br.delta, minima, build)
    return represented_factors(br, classes, search)


def factor(n, gap=4):
    br, shortcut = branch_from_n(n)
    if shortcut:
        return shortcut
    minima, direct = propose(br, denominator_budget(n), gap, Counter())
    if direct:
        return min(direct)
    hits = search_minima(br, minima, Counter(), Counter())
    return hits[0] if hits else None


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("--gap", type=int, default=4)
    args = parser.parse_args()
    print(factor(args.N, args.gap))
