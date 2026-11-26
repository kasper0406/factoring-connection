"""Simultaneous mixed-slope search on the norm-N shell of a public principal form."""

import json
import sys
from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt
from pathlib import Path

from flint import fmpz_mod_poly_ctx

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ellipse_aggregate"))
from factor import branch_from_n
from row_batches import product_tree


class WindowProducts:
    """P_d(t) = product(t-i, i=1..d), with a monic product tree over Z/nZ."""

    def __init__(self, n, work):
        self.ring = fmpz_mod_poly_ctx(n)
        self.work = work
        self.polynomials = {}

    def evaluate(self, width, points):
        if width not in self.polynomials:
            self.polynomials[width] = product_tree(
                [self.ring([-i, 1]) for i in range(1, width + 1)]
            )
            self.work["linear_factors_built"] += width
            self.work["polynomials_built"] += 1
        self.work["multipoint_values"] += len(points)
        return self.polynomials[width].multipoint_evaluate(points)


def small_factor(n, work):
    """Deterministic fourth-root block search. Returns n itself for a prime."""
    if n % 2 == 0:
        return 2
    limit = isqrt(n)
    if limit < 3:
        return n
    width = isqrt(limit - 1) or 1
    # The windows contain the integers 2 through 1+num_windows*width.
    num_windows = (limit - 1 + width - 1) // width
    starts = [1 + b * width for b in range(num_windows)]
    values = WindowProducts(n, work).evaluate(width, [-s for s in starts])
    for start, value in zip(starts, values):
        work["block_gcds"] += 1
        if gcd(int(value), n) > 1:
            for d in range(start + 1, min(start + width, limit) + 1):
                work["trial_remainders"] += 1
                if n % d == 0:
                    return d
    return n


def discriminant_divisors(delta, work):
    factors = Counter()
    pending = [delta]
    while pending:
        n = pending.pop()
        if n == 1:
            continue
        work["factor_nodes"] += 1
        d = small_factor(n, work)
        if d == n:
            factors[d] += 1
        else:
            pending.extend((d, n // d))
    divisors = [1]
    for p, exponent in sorted(factors.items()):
        old = divisors[:]
        for e in range(1, exponent + 1):
            power = p**e
            divisors.extend(d * power for d in old)
    return sorted(divisors)


def linear_hits(n, long_bound, short_bound, ratio, windows, work):
    """Find all (i,j) in the rectangle for which i-ratio*j is a nonunit."""
    size = isqrt(long_bound * short_bound)
    width = 1 << (size.bit_length() - 1)
    pending = [
        (j, base)
        for j in range(1, short_bound + 1)
        for base in range(0, long_bound, width)
    ]
    while pending:
        values = windows.evaluate(
            width, [(ratio * j - base) % n for j, base in pending]
        )
        children = []
        for (j, base), value in zip(pending, values):
            work["block_gcds"] += 1
            if gcd(int(value), n) == 1:
                continue
            work["nonunit_blocks"] += 1
            if width == 1:
                yield base + 1, j
            else:
                half = width // 2
                children.append((j, base))
                if base + half < long_bound:
                    children.append((j, base + half))
        pending = children
        width //= 2


def mixed_points(br, prep, search):
    """Return all positive mixed (w,y) using only a public eligible branch."""
    n, c, ell, delta = br.n, br.c, br.ell, br.delta
    if (
        delta <= 4
        or delta % 4 not in (0, 3)
        or min(ell, c) <= 0
        or (c - delta * ell) % 2
        or c * c + delta * ell * ell != 4 * n
        or gcd(c * ell * delta, n) != 1
    ):
        raise ValueError("expected an eligible primitive public norm-N branch")
    bound = (4 if delta % 8 == 3 else 8) * isqrt(2 * n)
    windows = WindowProducts(n, search)
    result = set()
    for u in discriminant_divisors(delta, prep):
        v = delta // u
        xmax, ymax = isqrt(bound // u), isqrt(bound // v)
        search["divisors"] += 1
        if not xmax or not ymax:
            continue
        search["nonempty_rectangles"] += 1
        search["rectangle_area"] += xmax * ymax
        swap = xmax < ymax
        long_coeff, short_coeff = (c, ell * u) if swap else (ell * u, c)
        ratio = short_coeff * pow(long_coeff, -1, n) % n
        for i, j in linear_hits(
            n, max(xmax, ymax), min(xmax, ymax), ratio, windows, search
        ):
            x, y = (j, i) if swap else (i, j)
            search["leaf_hits"] += 1
            h, k = u * x * x, v * y * y
            denominator = h + k
            if denominator > bound or gcd(u * x, y) != 1 or gcd(x, v) != 1:
                continue
            w, rw = divmod(c * (h - k) + 2 * delta * ell * x * y, denominator)
            t, rt = divmod(2 * c * x * y - ell * (h - k), denominator)
            if rw or rt or min(w, t) <= 0 or (w - delta * t) % 2:
                continue
            if w * w + delta * t * t != 4 * n:
                raise ArithmeticError("invalid reflected norm")
            if (ell * w - c * t) % n and (ell * w + c * t) % n:
                result.add((w, t))
    return sorted(result)


def slope_residues(br, points, stride):
    """Sign pairs with positive odd y, counted in the x/y chart modulo a 2-power."""
    if stride < 2 or stride & (stride - 1):
        raise ValueError("stride must be a power of two at least two")
    counts = Counter()
    for w, y in points:
        if y % 2 == 0:
            continue
        for signed_w in (w, -w):
            x = (signed_w - (br.delta % 2) * y) // 2
            counts[x * pow(y, -1, stride) % stride] += 1
    return counts


def block_summary(br, hits, stride, residue, prefix_modulus):
    """Count the entire selected slope coset and its two binary children."""
    if (
        prefix_modulus < 2
        or prefix_modulus & (prefix_modulus - 1)
        or stride < 2 * prefix_modulus
    ):
        raise ValueError(
            "need a dyadic prefix and a stride large enough for its two lifts"
        )
    residue %= prefix_modulus
    mixed = slope_residues(br, hits, stride)
    public = slope_residues(br, [(br.c, br.ell)], stride)

    def count(table, child):
        return sum(k for r, k in table.items() if r % (2 * prefix_modulus) == child)

    children = []
    for child in (residue, residue + prefix_modulus):
        known, extra = count(public, child), count(mixed, child)
        children.append(
            {
                "residue": child,
                "public_pairs": known,
                "mixed_pairs": extra,
                "signed_mod16": 2 * (known + extra) % 16,
            }
        )
    return {
        "residue": residue,
        "prefix_modulus": prefix_modulus,
        "fine_slopes_covered": stride // prefix_modulus,
        "mixed_pairs": sum(row["mixed_pairs"] for row in children),
        "signed_mod16": sum(row["signed_mod16"] for row in children) % 16,
        "children": children,
        "mixed_residues": sorted(r for r in mixed if r % prefix_modulus == residue),
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("--stride", type=int, default=65536)
    parser.add_argument("--prefix-modulus", type=int, default=8)
    parser.add_argument(
        "--residue", type=int, help="default: slope of the positive public point"
    )
    args = parser.parse_args()
    branch, shortcut = branch_from_n(args.N)
    if shortcut:
        print(json.dumps({"public_branch_gcd": shortcut}))
    else:
        prep, search = Counter(), Counter()
        hits = mixed_points(branch, prep, search)
        residue = (
            (branch.c - branch.delta % 2) // 2 if args.residue is None else args.residue
        )
        print(
            json.dumps(
                {
                    "N": branch.n,
                    "delta": branch.delta,
                    "mixed_points": hits,
                    "mixed_slope_multiplicities": slope_residues(
                        branch, hits, args.stride
                    ),
                    "raw_signed_mod16": 4 * (1 + len(hits)) % 16,
                    "group": block_summary(
                        branch, hits, args.stride, residue, args.prefix_modulus
                    ),
                    "discriminant_work": prep,
                    "search_work": search,
                },
                indent=2,
            )
        )
