"""Factor-or-empty decision for a prefix of mixed norms, via public classes."""

import json
import os
import sys
from argparse import ArgumentParser
from collections import Counter, defaultdict
from math import gcd, isqrt
from pathlib import Path

from bisector import WindowProducts, branch_from_n, public_factorization
from sympy.ntheory.residue_ntheory import sqrt_mod

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ellipse_global.benchmark import all_lifts, pari_points


def canonical(f):
    a, b, d = f
    return a, abs(b), d


def locally_square(f, factors):
    """Necessary for a square class: a represented square coprime to the discriminant."""
    a, b, d = f
    for p in factors:
        if p == 2:
            if not any(
                (a * x * x + b * x * y + d * y * y) % 8 == 1
                for x in range(8)
                for y in range(8)
            ):
                return False
        elif pow(a % p or d % p, (p - 1) // 2, p) != 1:
            return False
    return True


def reduced_forms_at(delta, a, work):
    """Generate primitive reduced forms at this leading coefficient up to inversion."""
    for root in sqrt_mod(-delta, 4 * a, all_roots=True):
        work["modular_roots"] += 1
        b = root % (2 * a)
        b = min(b, 2 * a - b)
        d = (b * b + delta) // (4 * a)
        if d >= a and gcd(a, b, d) == 1:
            yield a, b, d


def candidate_lifts(delta, lower, upper, gp, work):
    """Classes whose reduced square or inverse has minimum in [lower,upper]."""
    if lower > upper:
        return []
    factors = public_factorization(delta, work)
    candidates = set()
    for a in range(lower, upper + 1):
        work["indices"] += 1
        candidates.update(reduced_forms_at(delta, a, work))
    work["reduced_squares_to_test"] += len(candidates)
    problems = [(f, factors) for f in sorted(candidates) if locally_square(f, factors)]
    work["conics"] += len(problems)
    points, _ = pari_points(problems, gp)
    result = set()
    for (f, _), point in zip(problems, points, strict=True):
        if len(point) == 1:  # PARI's local obstruction to rational solubility
            work["rational_obstructions"] += 1
            continue
        roots, _ = all_lifts(f, factors, point)
        result.update(canonical(r) for r in roots)
    work["classes"] += len(result)
    return sorted(result)


def represented_factors(br, classes, work):
    """Batched verification; only exact prime values on the list count as hits."""
    n = br.n
    bound = isqrt(2 * n)
    root = br.c * pow(br.ell, -1, n) % n
    routes = []
    hits = set()
    for form in classes:
        a, b, d = form
        for value in (a, d):
            if 1 < value < n and n % value == 0:
                hits.add(value)
        xmax, ymax = isqrt(2 * bound // a), isqrt(2 * bound // d)
        if not ymax:
            continue
        inv = pow(2 * a, -1, n)  # a <= the public cap < both primes
        ratios = {sign * (b + u) * inv % n for sign in (-1, 1) for u in (root, -root)}
        routes.extend((form, ratio, xmax, ymax) for ratio in ratios)
    if not routes:
        return sorted(hits)
    area = sum(x * y for _, _, x, y in routes)
    side = sum(y for _, _, _, y in routes)
    work["signed_rectangle_area"] += area
    work["short_side_sum"] += side
    work["shape_budget"] += isqrt(area) + side
    pending = defaultdict(list)
    for f, rho, xmax, ymax in routes:
        width = 1 << (min(xmax, isqrt(area)).bit_length() - 1)
        for y in range(1, ymax + 1):
            target = rho * y % n
            pending[width].extend(
                (f, xmax, y, target, base) for base in range(0, xmax, width)
            )
    windows = WindowProducts(n, work)
    while pending:
        width = max(pending)
        items = pending.pop(width)
        vals = windows.evaluate(
            width, [(target - base) % n for _, _, _, target, base in items]
        )
        for (f, xmax, y, target, base), value in zip(items, vals, strict=True):
            if gcd(n, int(value)) == 1:
                continue
            if width == 1:
                x = base + 1
                a, b, d = f
                for z in (
                    a * x * x + b * x * y + d * y * y,
                    a * x * x - b * x * y + d * y * y,
                ):
                    if 1 < z < n and n % z == 0:
                        hits.add(z)
                work["leaf_hits"] += 1
            else:
                half = width // 2
                pending[half].append((f, xmax, y, target, base))
                if base + half < xmax:
                    pending[half].append((f, xmax, y, target, base + half))
    return sorted(hits)


def decision(br, low, high, gp):
    if not 1 <= low <= br.cap + 1 or not 0 <= high <= br.cap:
        raise ValueError("outside the eligible mixed-minimum cap")
    work = Counter()
    forms = candidate_lifts(br.delta, low, high, gp, work)
    hits = represented_factors(br, forms, work)
    return {"factor": hits[0] if hits else None, "work": dict(work)}


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--cutoff", type=int, help="default: full public cap")
    group.add_argument(
        "--only-minimum", type=int, help="test a single externally proposed minimum"
    )
    parser.add_argument("--gp", default=os.environ.get("GP", "gp"))
    opt = parser.parse_args()
    br, shortcut = branch_from_n(opt.N)
    if shortcut:
        answer = {"factor": shortcut}
    elif opt.only_minimum is not None:
        answer = decision(br, opt.only_minimum, opt.only_minimum, opt.gp)
    else:
        answer = decision(br, 1, br.cap if opt.cutoff is None else opt.cutoff, opt.gp)
    print(json.dumps(answer, indent=2))
