"""Parity in the primitive nonsquare upper band, via public form classes."""

import json
import os
from argparse import ArgumentParser
from collections import Counter, defaultdict
from math import gcd, isqrt, prod
from time import perf_counter

from class_prefix import (
    branch_from_n,
    canonical,
    lift_selected_squares,
    public_factorization,
    reduced_forms_at,
    represented_factors,
)
from rectangle_rays import Rectangles

from genusleak.forms import compose, reduce_form


def catalogue(delta, work):
    forms = set()
    for a in range(1, isqrt(delta // 3) + 1):
        work["parent_minimum_probes"] += 1
        forms.update(reduced_forms_at(delta, a, work))
    squares = defaultdict(set)
    for f in forms:
        squares[canonical(compose(f, f))].add(f)
    work["parent_classes_up_to_inverse"] += len(forms)
    work["parent_squares_up_to_inverse"] += len(squares)
    return squares


def restrict_form(form, line, prime):
    a, b, c = form
    if line == prime:  # y = 0 mod prime
        restricted = (a, b * prime, c * prime**2)
    else:  # x = line*y mod prime
        restricted = (
            a * prime**2,
            prime * (2 * a * line + b),
            a * line**2 + b * line + c,
        )
    if gcd(*restricted) != 1:
        return None
    return canonical(reduce_form(restricted))


def band_vectors(form, parent_delta, cutoff, bad_primes, work):
    """Parent vectors in the nonsquare band, up to sign and primitive away from E."""
    a, b, c = form
    for y in range(isqrt(4 * a * cutoff // parent_delta) + 1):
        s = isqrt(4 * a * cutoff - parent_delta * y * y)
        lo, hi = -((s + b * y) // (2 * a)), (s - b * y) // (2 * a)
        if y == 0:
            lo = max(1, lo)
        work["parent_rows"] += 1
        for x in range(lo, hi + 1):
            work["parent_vectors"] += 1
            t = a * x * x + b * x * y + c * y * y
            if 4 * t <= cutoff or isqrt(t) ** 2 == t:
                continue
            content = gcd(x, y)
            g = content
            while (d := gcd(g, bad_primes)) > 1:
                g //= d
            if g != 1:
                continue
            work["retained_vectors"] += 1
            yield t, x, y, content


def children_with_odd_band(form, parent_delta, prime, cutoff, bad_primes, work):
    """Route all short parent vectors to their incident child lines, up to sign."""
    odd_lines = set()
    common = 0
    directions = {}
    for _, x, y, content in band_vectors(form, parent_delta, cutoff, bad_primes, work):
        if prime == 1:
            line = 0
        elif y % prime:
            line = x * pow(y, -1, prime) % prime
        elif x % prime:
            line = prime
        else:
            common ^= 1
            work["vectors_on_all_lines"] += 1
            continue
        direction = (x // content, y // content)
        if line in directions and directions[line] != direction:
            work["nonparallel_revisits"] += 1
        directions[line] = direction
        odd_lines.symmetric_difference_update((line,))
    if prime == 1:
        return {form} if odd_lines else set()
    if common:
        work["complement_lines_materialized"] += prime + 1
        odd_lines = set(range(prime + 1)) - odd_lines
    result = set()
    work["odd_lines_before_primitivity"] += len(odd_lines)
    for line in odd_lines:
        child = restrict_form(form, line, prime)
        if child is not None:
            result.add(child)
    return result


def parent_catalogue(br, work):
    factors = public_factorization(br.delta, work)
    # E uses the physical coefficient, Delta/4 on the even branch.
    bad_primes = prod(p for p, e in factors.items() if e - (2 if p == 2 else 0) >= 2)
    prime = max((p for p, e in factors.items() if p != 2 and e >= 2), default=1)
    parent_delta = br.delta // prime**2
    squares = catalogue(parent_delta, work)
    return factors, prime, parent_delta, bad_primes, squares


def candidates(br, cutoff, gp, work):
    """Public class selector; no tests modulo N at this stage."""
    if not 0 <= cutoff <= br.cap:
        raise ValueError("cutoff outside first-minimum range")
    if not cutoff:
        return [], 1
    factors, prime, parent_delta, bad_primes, squares = parent_catalogue(br, work)
    selected = set()
    for parent in squares:
        selected.update(
            children_with_odd_band(
                parent, parent_delta, prime, cutoff, bad_primes, work
            )
        )
    work["selected_square_candidates"] = len(selected)
    if prime == 1:
        result = sorted({root for f in selected for root in squares[f]})
        work["classes"] = len(result)
    else:
        result = lift_selected_squares(selected, factors, gp, work)
    return result, prime


def verify_levels(br, cutoff):
    if not 0 <= cutoff <= br.cap:
        raise ValueError("cutoff outside first-minimum range")
    work = Counter()
    if cutoff < 2:
        return {
            "parity": 0,
            "exact_count": 0,
            "conductor_step": None,
            "build": {},
            "verification": {},
            "build_seconds": 0.0,
            "verification_seconds": 0.0,
        }
    start = perf_counter()
    factors, prime, parent_delta, bad_primes, squares = parent_catalogue(br, work)
    levels = set()
    for f in squares:
        for t, x, y, _ in band_vectors(f, parent_delta, cutoff, bad_primes, work):
            # If the line is isotropic, its restriction is imprimitive. A vector
            # zero modulo the conductor prime is incident to every line.
            if prime == 1 or t % prime or (x % prime == 0 and y % prime == 0):
                levels.add(t)
    work["proposed_levels"] = len(levels)
    d_factors = factors.copy()
    if br.delta % 4 == 0:
        d_factors[2] -= 2
        if d_factors[2] == 0:
            del d_factors[2]
    middle = perf_counter()
    rectangles = Rectangles(br, d_factors=d_factors)
    exact = rectangles.count(
        cutoff // 4 + 1,
        cutoff,
        modulus=None,
        base_nonsquares=True,
        proposed_levels=levels,
    )
    end = perf_counter()
    return {
        "parity": exact % 2,
        "exact_count": exact,
        "conductor_step": prime,
        "build": dict(work),
        "verification": rectangles.metrics(),
        "build_seconds": middle - start,
        "verification_seconds": end - middle,
    }


def evaluate(br, cutoff, gp):
    build, search = Counter(), Counter()
    start = perf_counter()
    classes, prime = candidates(br, cutoff, gp, build)
    middle = perf_counter()
    hits = represented_factors(br, classes, search)
    end = perf_counter()
    return {
        "parity": int(bool(hits)),
        "conductor_step": prime,
        "factor_witness": hits[0] if hits else None,
        "build": dict(build),
        "membership": dict(search),
        "build_seconds": middle - start,
        "membership_seconds": end - middle,
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("cutoff", type=int)
    parser.add_argument("--class-membership", action="store_true")
    parser.add_argument("--gp", default=os.environ.get("GP", "gp"))
    args = parser.parse_args()
    branch, shortcut = branch_from_n(args.N)
    if branch is None:
        print(json.dumps({"branch_factor": shortcut}))
    else:
        result = (
            evaluate(branch, args.cutoff, args.gp)
            if args.class_membership
            else verify_levels(branch, args.cutoff)
        )
        print(json.dumps(result, indent=2))
