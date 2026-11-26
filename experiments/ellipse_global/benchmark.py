"""Public ellipse counts, inverse Gauss duplication, and independent toy controls.

Run at repo root with Python 3.10+, NumPy, SymPy and PARI/GP 2.19.1.
Set GP to the executable path if it is outside PATH. See the accompanying note.
"""

import argparse
import ast
import json
import os
import random
import subprocess
import sys
from bisect import bisect_right
from collections import Counter
from itertools import combinations
from math import gcd, isqrt, prod
from pathlib import Path
from time import perf_counter

import numpy as np
from sympy import cyclotomic_poly, factorint, nextprime, primerange, symbols
from sympy.polys.polytools import Poly

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments/nov26_consolidation")]
from genusleak.forms import (
    compose,
    disc,
    identity,
    inverse,
    is_reduced,
    power,
    prime_form,
    reduce_form,
    sqrt_mod_prime,
)
from check_identities import count, extract, first, forms


def public_shells(n, delta, limit):
    """Return raw signed counts for each j. No roots modulo n or divisors of n."""
    parity = delta % 2
    shells = [0] * (limit + 1)
    tests = 0
    for j in range(1, limit + 1):
        radicand = 4 * n * j
        for y in range(isqrt(radicand // delta) + 1):
            v = radicand - delta * y * y
            s = isqrt(v)
            tests += 1
            if s * s == v and (s - parity * y) % 2 == 0:
                shells[j] += (1 if y == 0 else 2) * (1 if s == 0 else 2)
    return shells, tests


def direct_shells(n, delta, limit):
    """Independent point enumeration on toy inputs; signed vectors including both y signs."""
    b = delta % 2
    d = (b * b + delta) // 4
    result = np.zeros(limit + 1, dtype=np.int64)
    bound = isqrt(4 * n * limit // delta)
    for y in range(-bound, bound + 1):
        s = isqrt(4 * n * limit - delta * y * y)
        x = np.arange(-((s + b * y) // 2), (s - b * y) // 2 + 1, dtype=np.int64)
        values = x * x + b * y * x + d * y * y
        hits = values[(values > 0) & (values % n == 0)] // n
        result += np.bincount(hits, minlength=limit + 1)
    return result.tolist()


def principal_pairs(delta, limit):
    """Exact count up to sign, using intervals, including ties."""
    b = delta % 2
    total = isqrt(limit)
    for y in range(1, isqrt(4 * limit // delta) + 1):
        s = isqrt(4 * limit - delta * y * y)
        total += (s - b * y) // 2 + (s + b * y) // 2 + 1
    return total


def decode(delta, raw):
    def residue(b):
        v = (raw[b] - 4 * principal_pairs(delta, b)) % 16
        assert v % 4 == 0
        return v // 4

    hi = 1
    while residue(hi) == 0:
        hi *= 2
    a = first(0, hi, lambda b: residue(b) != 0)
    if a == 1:
        return identity(-delta)
    hi = (delta + a * a) // (4 * a)
    d = first(a - 1, hi, lambda b: (residue(b) - isqrt(b // a)) % 4 != 0)
    b = isqrt(4 * a * d - delta)
    assert b * b == 4 * a * d - delta
    return (a, b, d)


def clipped_control(n, ell, c, delta, raw, square_form):
    """Public enumeration simulates a clipped oracle; the reduction receives residues only."""
    m = (c - (delta % 2) * ell) // 2
    prefix_cache = {}
    radial_calls = 0
    clip_calls = 0

    def oracle(B, cutoff):
        nonlocal clip_calls
        clip_calls += 1
        if B not in prefix_cache:
            vals = []
            for j in range(1, B + 1):
                for y in range(
                    -isqrt(4 * n * j // delta), isqrt(4 * n * j // delta) + 1
                ):
                    v = 4 * n * j - delta * y * y
                    s = isqrt(v)
                    if s * s != v or (s - (delta % 2) * y) % 2:
                        continue
                    for k in (s, -s) if s else (0,):
                        x = (k - (delta % 2) * y) // 2
                        L = m * y - ell * x
                        if L > 0:
                            vals.append(L)
            prefix_cache[B] = sorted(vals)
        return bisect_right(prefix_cache[B], cutoff) % 8

    def end(B):
        return isqrt(4 * n * n * B // delta)

    def has_mixed(B):
        nonlocal radial_calls
        radial_calls += 1
        R = (2 * oracle(B, end(B)) + 2 * isqrt(B)) % 16
        assert R == raw[B] % 16
        return (R - 4 * principal_pairs(delta, B)) % 16 != 0

    hi = 1
    while not has_mixed(hi):
        hi *= 2
    a = first(0, hi, has_mixed)
    assert a == square_form[0] and 4 * a < delta

    def has_positive(L):
        principal = min(isqrt(a), L // (c * ell))
        return (oracle(a, L) - principal) % 8 != 0

    minimum = first(0, end(a), has_positive)
    factor = gcd(n, minimum)
    assert 1 < factor < n and n % factor == 0
    return dict(
        minimum_positive_coordinate=minimum,
        divisor=factor,
        total_calls=clip_calls,
        radial_calls=radial_calls,
    )


def bezout(a, b):
    x, y, u, v = 1, 0, 0, 1
    while b:
        q, r = divmod(a, b)
        a, b, x, u, y, v = b, r, u, x - q * u, v, y - q * v
    return (x, y) if a == 1 else (-x, -y)


def public_product_rows(n, ell, c, delta, limit):
    """Rejected scan in the public norm-n basis; row gcds may themselves reveal a divisor."""
    b0 = delta % 2
    m = (c - b0 * ell) // 2
    t, s = bezout(m, ell)
    assert m * t + ell * s == 1
    r = -s
    d0 = (delta + b0 * b0) // 4
    h = 2 * m * r + b0 * (m * t + ell * r) + 2 * d0 * ell * t
    k = r * r + b0 * r * t + d0 * t * t
    assert 4 * n * k - h * h == delta and gcd(h, n) == 1
    hinv = pow(h, -1, n)
    total = 2 * isqrt(limit)
    rows = isqrt(4 * n * n * limit // delta)
    for y in range(1, rows + 1):
        mod = n // gcd(y, n)
        residue = (-k * y * hinv) % mod
        rad = isqrt(4 * n * n * limit - delta * y * y)
        lo = -((rad + h * y) // (2 * n))
        hi = (rad - h * y) // (2 * n)
        total += 2 * ((hi - residue) // mod - (lo - 1 - residue) // mod)
    return total, rows


def branches(n):
    top = isqrt(4 * n)
    for ell in (3, 5, 7, 11, 13, 17, 19):
        if n % ell == 0:
            continue
        step = ell * ell
        for r in range(1, step):
            if (r * r - 4 * n) % step:
                continue
            c = top - (top - r) % step
            if c <= 0:
                continue
            delta = (4 * n - c * c) // step
            if delta > 4 and gcd(n, delta) == 1:
                yield (ell, c, delta)


def cohort():
    rng = random.Random(20261126)
    pairs = [
        (int(p), int(q))
        for p, q in combinations(primerange(41, 242), 2)
        if q <= 2 * p and p * q not in (11009, 28829)
    ]
    rng.shuffle(pairs)
    cases = [(101, 109, 11, 204, 20), (127, 227, 11, 336, 20)]
    for p, q in pairs:
        options = list(branches(p * q))
        if options:
            cases.append((p, q, *rng.choice(options)))
            if len(cases) == 194:
                break
    assert len(cases) == 194
    return cases


def ambiguous(delta, factors):
    """Candidates from the divisors of supplied factored delta. Includes nonmaximal orders."""
    ds = [1]
    for p, e in factors.items():
        exponents = range(e + 1) if p == 2 else (0, e)
        ds = [d * p**k for d in ds for k in exponents]
    out = set()
    for u in ds:
        v = delta // u
        candidates = []
        if v % 4 == 0:
            candidates.append((u, 0, v // 4))
        if (u + v) % 4 == 0:
            candidates.append((u, u, (u + v) // 4))
            if v >= u and (v - u) % 2 == 0:
                a = (u + v) // 4
                candidates.append((a, (v - u) // 2, a))
        for f in candidates:
            if is_reduced(f) and gcd(*f) == 1:
                out.add(f)
    return sorted(out)


def pari_points(problems, gp):
    """Ask PARI for rational conic points using the supplied factorization of 2*Delta."""
    if not problems:
        return [], 0
    src = []
    for (a, b, c), factors in problems:
        detfac = factors.copy()
        detfac[2] = 1 + detfac.get(2, 0)
        f = "Mat([" + ";".join(f"{p},{e}" for p, e in sorted(detfac.items())) + "])"
        g = f"[{2 * a},{b},0;{b},{2 * c},0;0,0,-2]"
        src.append(f"print(Vec(qfsolve([{g},{f}])))")
    t = perf_counter()
    proc = subprocess.run(
        [gp, "-fq"],
        input="\n".join(src) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    elapsed = perf_counter() - t
    if proc.stderr.strip():
        raise RuntimeError(proc.stderr)
    points = [ast.literal_eval(line) for line in proc.stdout.splitlines()]
    assert len(points) == len(problems)
    return points, elapsed


def root_from_point(f, point):
    """Inverse duplication at a primitive square representation."""
    a, b, c = f
    x, y, z = point
    g = gcd(x, y)
    assert z % g == 0
    x, y, z = x // g, y // g, abs(z) // g
    assert a * x * x + b * x * y + c * y * y == z * z
    u, v = bezout(x, y)
    r, s = -v, u
    assert x * s - y * r == 1
    middle = 2 * a * x * r + b * (x * s + y * r) + 2 * c * y * s
    last = a * r * r + b * r * s + c * s * s
    if gcd(z, middle) != 1:
        return None
    candidate = reduce_form((z, middle, z * last))
    assert compose(candidate, candidate) == f
    return candidate


def unit_sqrt(value, p, exponent):
    """Hensel lift a square unit at one supplied prime. At 2, require 1 mod 8."""
    if p == 2:
        assert value % 8 == 1
        x = 1
        for k in range(3, exponent):
            if (x * x - value) % (1 << (k + 1)):
                x += 1 << (k - 1)
    else:
        x = sqrt_mod_prime(value, p)
        inv = pow(2 * x, -1, p)
        mod = p
        for _ in range(1, exponent):
            x += mod * (((value - x * x) // mod * inv) % p)
            mod *= p
    assert (x * x - value) % (p**exponent) == 0
    return x


def first_root(f, point, factors):
    root = root_from_point(f, point)
    if root:
        return root, 0
    # Put the primitive conic point at (1,0,z); the binary form becomes [z^2,b,d].
    a, b0, c = f
    x, y, z = point
    g = gcd(x, y)
    x, y, z = x // g, y // g, abs(z) // g
    u, v = bezout(x, y)
    r, t = -v, u
    b = 2 * a * x * r + b0 * (x * t + y * r) + 2 * c * y * t
    d = a * r * r + b0 * r * t + c * t * t
    # CRT a secant slope S/T on the conic to obtain a primitive square norm
    # coprime to Delta, even when the original norm meets the conductor.
    S, T, mod = 0, 0, 1
    bad_count = 0
    for p, e in factors.items():
        if z % p:
            M, local_s, local_t = p, 1, 0
        else:
            bad_count += 1
            M = p ** (e + 1)
            w = 0
            if p == 2:
                w = next(w for w in range(8) if (z * z * w * w + b * w + d) % 8 == 1)
            s = unit_sqrt(z * z * w * w + b * w + d, p, e + 1)
            local_s, local_t = s - z * w, 1
        inv = pow(mod, -1, M)
        S += mod * (((local_s - S) * inv) % M)
        T += mod * (((local_t - T) * inv) % M)
        mod *= M
    X = d * T * T - S * S
    Y = 2 * z * S * T - b * T * T
    Z = z * (d * T * T + S * S) - b * S * T
    g = gcd(X, Y)
    assert Z % g == 0 and gcd(Z // g, -disc(f)) == 1
    root = root_from_point(f, (x * X + r * Y, y * X + t * Y, Z))
    assert root is not None
    return root, bad_count


def all_lifts(f, factors, point):
    seed, repairs = first_root(f, point, factors)
    torsion = ambiguous(-disc(f), factors)
    roots = sorted({compose(seed, t) for t in torsion})
    assert len(roots) == len(torsion) and all(compose(r, r) == f for r in roots)
    return roots, repairs


def fourier_check(n):
    """Exact cyclotomic DFT of the public multiplicity minus the two public lines."""
    delta = 4 * n - 1  # test root 1, split and unramified at both primes
    b = 1
    d = (delta + 1) // 4
    pts = []
    weight = {}
    for x in range(n):
        for y in range(n):
            if (x * x + b * x * y + d * y * y) % n:
                continue
            g = gcd(x, y, n)
            w = (
                (4 if g == n else 1 if g == 1 else 2)
                - int(x % n == 0)
                - int((x + y) % n == 0)
            )
            if w:
                pts.append((x, y, w))
                weight[x, y] = w
    t = symbols("t")
    coeff = list(map(int, reversed(Poly(cyclotomic_poly(n, t), t).all_coeffs())))
    degree = len(coeff) - 1
    active = [(i, x) for i, x in enumerate(coeff) if x]
    for u in range(n):
        for v in range(n):
            hist = [0] * n
            for x, y, w in pts:
                hist[(u * y - v * x) % n] += w
            hist[0] -= n * weight.get((u, v), 0)
            for j in range(n - 1, degree - 1, -1):
                c = hist[j]
                if c:
                    for i, k in active:
                        hist[j - degree + i] -= k * c
            assert not any(hist)
    return {
        "N": n,
        "frequencies_verified": n * n,
        "unit_mixed_coefficients_on_primitive_frequencies": sum(
            gcd(x, y, n) == 1 and w == 1 for x, y, w in pts
        ),
    }


def lifting_stress(gp):
    # Small ground truth enumerates every class, independent of the conic solver.
    problems, expected = [], []
    for delta in range(5, 1201):
        fibers = {}
        for r in forms(delta):
            fibers.setdefault(compose(r, r), set()).add(r)
        if fibers:
            fac = dict(factorint(delta))
            for f, preimages in fibers.items():
                problems.append((f, fac))
                expected.append(preimages)
    points, _ = pari_points(problems, gp)
    repaired, size = 0, 0
    for (f, fac), want, point in zip(problems, expected, points):
        roots, corrections = all_lifts(f, fac, point)
        assert set(roots) == want
        repaired += int(corrections > 0)
        size += len(roots)

    # Larger orders have supplied, known factorizations. Choose the seed class
    # before computing its square, and keep that seed out of the lifting input.
    rng = random.Random(891463)
    big, control = [], []
    for bits in (32, 64, 96, 128):
        for i in range(8):
            num = (2, 3, 4, 5)[i % 4]
            exps = [1] * num
            if i % 3 != 0:
                exps[0] = 2
            width = max(5, (bits - 3) // sum(exps))
            fac = {}
            for e in exps:
                while True:
                    p = int(nextprime(rng.randrange(1 << (width - 1), 1 << width)))
                    if p not in fac:
                        break
                fac[p] = e
            delta = prod(p**e for p, e in fac.items())
            if i % 2 or delta % 4 != 3:
                fac[2] = 2 + i % 4
                delta *= 2 ** fac[2]
            root = identity(-delta)
            split = [
                s for s in primerange(3, 400) if pow(-delta % s, (s - 1) // 2, s) == 1
            ]
            for s in split[:4]:
                root = compose(
                    root, power(prime_form(-delta, s), rng.randrange(17, 2000))
                )
            big.append((compose(root, root), fac))
            control.append(root)
    # A separate higher-2-rank order measures output-size cost.
    fac = {p: 1 for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)}
    fac[2] = 3
    delta = prod(p**e for p, e in fac.items())
    root = prime_form(
        -delta,
        next(s for s in primerange(41, 400) if pow(-delta % s, (s - 1) // 2, s) == 1),
    )
    big.append((compose(root, root), fac))
    control.append(root)
    points, elapsed = pari_points(big, gp)
    records = []
    for (f, fac), want, point in zip(big, control, points):
        start = perf_counter()
        roots, corrections = all_lifts(f, fac, point)
        seconds = perf_counter() - start
        assert want in roots
        records.append(
            dict(
                delta=-disc(f),
                bits=(-disc(f)).bit_length(),
                squared_form=f,
                roots=len(roots),
                repaired_primes=corrections,
                seconds_after_conic=seconds,
            )
        )
    return dict(
        small_discriminant_limit=1200,
        small_square_classes=len(problems),
        small_preimages=size,
        small_cases_repaired=repaired,
        larger_conics_seconds=elapsed,
        larger=records,
    )


def benchmark(gp):
    cases = cohort()
    weights = [(u, v) for u in range(1, 25) for v in range(1, 25) if gcd(u, v) == 1]
    rows = []
    problems = []
    factor_time = 0
    for idx, (p, q, ell, c, delta) in enumerate(cases):
        n = p * q
        limit = min(
            isqrt(n // 2),
            max(20 if idx < 2 else 1, isqrt(4 * delta // 3) + 1, (delta + 4) // 8),
        )
        shells, steps = public_shells(n, delta, limit)
        assert shells == direct_shells(n, delta, limit)
        raw = list(map(int, np.cumsum(shells)))
        f = decode(delta, raw)
        # Holdout truth uses a prime ideal over p only; no mixed root is needed anywhere.
        truth = prime_form(-delta, p)
        square = compose(truth, truth)
        assert square in (f, inverse(f))
        assert all(
            raw[j] == 4 * (principal_pairs(delta, j) + count(f, j))
            for j in range(1, limit + 1)
        )
        alt, ys = public_product_rows(n, ell, c, delta, limit)
        assert alt == raw[-1]
        start = perf_counter()
        factors = dict(factorint(delta))
        factor_time += perf_counter() - start
        problems.append((f, factors))
        rows.append(
            dict(
                N=n,
                p=p,
                q=q,
                ell=ell,
                c=c,
                delta=delta,
                limit=limit,
                signed_counts=raw[1:],
                square_form=f,
                true_min=truth[0],
                quotient_tests=steps,
                product_rows=ys,
                clipped=clipped_control(n, ell, c, delta, raw, f),
            )
        )
    points, pari_time = pari_points(problems, gp)
    lift_time = verify_time = 0
    for row, (f, factors), point in zip(rows, problems, points):
        begin = perf_counter()
        roots, repairs = all_lifts(f, factors, point)
        lift_time += perf_counter() - begin
        assert set(roots) == {r for r in forms(row["delta"]) if compose(r, r) == f}
        mins = sorted({r[0] for r in roots})
        assert row["true_min"] in mins
        begin = perf_counter()
        witnesses = {}
        for j in mins:
            for u, v in weights:
                split = extract(row["N"], row["ell"], row["c"], j, u, v)
                if split is not None:
                    assert split in (row["p"], row["q"])
                    witnesses.setdefault(j, [u, v, split])
        verify_time += perf_counter() - begin
        row.update(
            lifts=roots,
            minima=mins,
            repaired_primes=repairs,
            witnesses=witnesses,
            integer_square_tests=len(mins) * len(weights),
        )

    def summary(rows):
        return dict(
            semiprimes=len(rows),
            cutoffs=sum(r["limit"] for r in rows),
            square_tests=sum(r["integer_square_tests"] for r in rows),
            hit_branches=sum(bool(r["witnesses"]) for r in rows),
            nonprincipal_squares=sum(r["square_form"][0] > 1 for r in rows),
            clipped_extractions=len(rows),
            maximum_clipped_calls=max(r["clipped"]["total_calls"] for r in rows),
            lift_histogram=dict(sorted(Counter(len(r["lifts"]) for r in rows).items())),
            quotient_tests=sum(r["quotient_tests"] for r in rows),
            product_rows=sum(r["product_rows"] for r in rows),
            branches_with_nonzero_mixed_residue=sum(
                any(
                    (x - 4 * principal_pairs(r["delta"], j)) % 16
                    for j, x in enumerate(r["signed_counts"], 1)
                )
                for r in rows
            ),
            branches_with_multiple_passing_minima=sum(
                len(r["witnesses"]) > 1 for r in rows
            ),
        )

    return dict(
        repo_revision="9e70736b4abaacff226c93e4fdddaf5e6fcad5f5",
        preserved=summary(rows[:2]),
        replacement=summary(rows[2:]),
        rows=rows,
        fourier=[fourier_check(n) for n in (15, 35, 55)],
        lifting_stress=lifting_stress(gp),
        timings_seconds=dict(
            factoring_toy_discriminants=factor_time,
            conics_including_one_gp_startup=pari_time,
            inverse_duplication_torsion_and_verification=lift_time,
            full_weight_grid=verify_time,
        ),
        weights=len(weights),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = benchmark(os.environ.get("GP", "gp"))
    args.output.write_text(json.dumps(data, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in data.items() if k not in ("rows", "lifting_stress")},
            indent=2,
        )
    )
