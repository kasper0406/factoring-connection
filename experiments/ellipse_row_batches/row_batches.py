"""Whole-block and quadratic-moment experiments. No prime factors in evaluator inputs."""

import argparse
import json
from bisect import bisect_left, bisect_right
from collections import Counter
from dataclasses import dataclass
from functools import cache
from math import gcd, isqrt
from pathlib import Path

from flint import fmpz_mod_poly_ctx
from arc_batches import arc_shell, arc_cost


@dataclass(frozen=True)
class Branch:
    n: int
    ell: int
    c: int
    delta: int

    @property
    def cap(self):
        return isqrt(self.delta // 3)

    def height(self, B):
        return isqrt(4 * self.n**2 * B // self.delta)


def ceil_div(a, b):
    return -((-a) // b)


def interval_union(intervals):
    merged = []
    for lo, hi in sorted((a, b) for a, b in intervals if a <= b):
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
        else:
            merged.append((lo, hi))
    return merged


def row_intervals(br, j, lo, hi):
    n, ell, c, D = br.n, br.ell, br.c, br.delta
    hi = min(hi, br.height(j))
    lo = max(lo, 1)
    if lo > hi:
        return []
    left = 4 * n * n * j - D * lo * lo
    right = 4 * n * n * j - D * hi * hi
    # Outward square-root roundings; monotonicity on L>=0.
    big = isqrt(left)
    big += big * big < left
    small = isqrt(right)
    Y = isqrt(4 * n * j // D)
    intervals = [
        (ceil_div(c * lo + ell * small, 2 * n), (c * hi + ell * big) // (2 * n)),
        (ceil_div(c * lo - ell * big, 2 * n), (c * hi - ell * small) // (2 * n)),
    ]
    return interval_union([(max(-Y, l), min(Y, h)) for l, h in intervals])


def projected_shell(br, j, lo, hi, work, intervals=None):
    """Exact positive determinant rows in [lo,hi], with multiplicity."""
    n, ell, c, D = br.n, br.ell, br.c, br.delta
    if intervals is None:
        intervals = row_intervals(br, j, lo, hi)
    work["shell_intervals"] += 1
    ans = []
    for low, high in intervals:
        for y in range(low, high + 1):
            work["square_tests"] += 1
            v = 4 * n * j - D * y * y
            w = isqrt(v)
            if w * w != v or (w - D * y) % 2:
                continue
            for sign_w in (w, -w) if w else (0,):
                L = (c * y - ell * sign_w) // 2
                if lo <= L <= hi:
                    ans.append(L)
    return ans


def hybrid_shell(br, j, lo, hi, work):
    intervals = row_intervals(br, j, lo, hi)
    if sum(h - l + 1 for l, h in intervals) <= arc_cost(br, j):
        work["projected_shells"] += 1
        return projected_shell(br, j, lo, hi, work, intervals)
    work["arc_shells"] += 1
    return arc_shell(br, j, lo, hi, work)


def disk_query(br, B, interval=None):
    lo, hi = interval if interval is not None else (1, br.height(B))
    work = Counter()
    count = moment = 0
    for j in range(1, B + 1):
        rows = hybrid_shell(br, j, lo, hi, work)
        count += len(rows)
        moment = (moment + sum(L * L for L in rows)) % br.n
    answer = dict(
        positive_count_mod8=count % 8,
        positive_second_moment_mod_N=moment,
        work=dict(work),
    )
    if interval is None:
        answer["signed_count_mod16"] = (2 * count + 2 * isqrt(B)) % 16
    return answer


def shell_reference(br, j):
    """Independent public enumeration using x rather than projection onto y."""
    n, D = br.n, br.delta
    b = D % 2
    d = (D + b * b) // 4
    X = isqrt(4 * n * j * d // D)
    ans = []
    for x in range(-X, X + 1):
        v = 4 * d * n * j - D * x * x
        s = isqrt(v)
        if s * s != v:
            continue
        for w in (s, -s) if s else (0,):
            y, rem = divmod(w - b * x, 2 * d)
            if rem:
                continue
            L = (br.c - b * br.ell) // 2 * y - br.ell * x
            if L > 0:
                ans.append(L)
    return sorted(ans)


def product_tree(nodes):
    while len(nodes) > 1:
        nxt = [nodes[i] * nodes[i + 1] for i in range(0, len(nodes) - 1, 2)]
        if len(nodes) % 2:
            nxt.append(nodes[-1])
        nodes = nxt
    return nodes[0]


def screened_shell(br, j, lo, hi, size, work):
    """Boolean block-unit filter, followed by shell arcs on the surviving intervals."""
    assert 1 <= j <= br.cap
    hi = min(hi, br.height(j))
    lo = max(1, lo)
    if hi < lo:
        return []
    unit_step = br.c * br.ell
    k = isqrt(j)
    principal = k * unit_step if k * k == j else None
    result = [principal] if principal is not None and lo <= principal <= hi else []
    full, tail = divmod(hi - lo + 1, size)
    blocks = []
    if full:
        ring = fmpz_mod_poly_ctx(br.n)
        F = product_tree([ring([i, 1]) for i in range(1, size + 1)])
        starts = list(range(lo - 1, lo - 1 + full * size, size))
        vals = F.multipoint_evaluate(starts)
        blocks.extend((u + 1, u + size, int(v)) for u, v in zip(starts, vals))
        work["polynomial_factors"] += size
        work["multipoint_values"] += full
    if tail:
        start = hi - tail + 1
        v = 1
        for u in range(start, hi + 1):
            v = v * u % br.n
        blocks.append((start, hi, v))
        work["tail_modular_products"] += tail
    for l, h, v in blocks:
        work["blocks"] += 1
        if gcd(v, br.n) == 1:
            work["unit_blocks"] += 1
            continue
        work["dirty_blocks"] += 1
        vals = projected_shell(br, j, l, h, work)
        # A surviving block may also contain the public conjugate vector.
        if principal is not None and l <= principal <= h:
            vals.remove(principal)
        if not vals:
            work["dirty_empty"] += 1
        result.extend(vals)
    return sorted(result)


def first_visible(br, scalar):
    """Uses only the normalized mixed second moment, at capped public cutoffs."""
    calls = []

    @cache
    def value(B):
        calls.append(B)
        return scalar(B)

    for B in range(1, min(7, br.cap) + 1):
        if value(B):
            return B, value(B), calls
    lo, hi = 7, min(8, br.cap)
    while hi > lo:
        if value(hi):
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if value(mid):
                    hi = mid
                else:
                    lo = mid
            return hi, value(hi), calls
        lo, hi = hi, min(2 * hi, br.cap)
    raise ArithmeticError("first shell is outside the promised search range")


def shell_ys(br, a, S):
    """Possible positive y coordinates. Exact tests cover both shell multiplicities."""
    t = isqrt(S)
    if t * t == S and t:
        yield t
    b2 = 4 * a * a - br.delta
    if b2 <= 0:
        return
    b = isqrt(b2)
    if b * b != b2 or b >= a:
        return
    P, rem = divmod(a * S - br.n, b)
    gap = S * S - 4 * P * P
    if rem or gap < 0:
        return
    w = isqrt(gap)
    if w * w != gap or (S - w) % 2:
        return
    for yy in ((S - w) // 2, (S + w) // 2):
        y = isqrt(yy)
        if yy and y * y == yy:
            yield y


def split(br, a, S):
    for y in shell_ys(br, a, S):
        v = 4 * br.n * a - br.delta * y * y
        if v < 0:
            continue
        w = isqrt(v)
        if w * w != v or (w - br.delta * y) % 2:
            continue
        L = (br.c * y - br.ell * w) // 2
        g = gcd(L, br.n)
        if 1 < g < br.n:
            return g
    raise ArithmeticError("moment did not give a verified divisor")


def scalar_from_raw(br, B, raw_second):
    k = isqrt(B)
    principal = (br.c * br.ell) ** 2 * (k * (k + 1) * (2 * k + 1) // 6)
    return (raw_second - principal) * pow(br.c, -2, br.n) % br.n


def primes():
    p = 2
    while True:
        if all(p % d for d in range(2, isqrt(p) + 1)):
            yield p
        p += 1


def periodic_decode(br, a, oracle):
    """Recover moment from full signed residue buckets mod 16. No row endpoints."""
    T = br.height(a)
    step = br.c * br.ell
    root = isqrt(a)
    acc, mod, calls = 0, 1, 0
    for p in primes():
        histogram = [oracle(p, k) for k in range(p)]
        calls += p
        # On the exact shell only the +/- conjugates of k*alpha occur in
        # the two principal directions, and only when a=k^2.
        if a == root * root:
            for r in (0, 0, root * step, -root * step):
                histogram[r % p] -= 1
        histogram = [h % 16 for h in histogram]
        assert sum(histogram) in (4, 8)
        moment = sum(i * i * n for i, n in enumerate(histogram)) % p
        acc += mod * ((moment - acc) * pow(mod, -1, p) % p)
        mod *= p
        if mod > 8 * T * T:
            assert acc % 2 == 0
            S = acc // 2 * pow(br.c, -2, br.n) % br.n
            return S, calls, p


def native_moment_factor(br):
    """Streaming public computation, storing cumulative moments but no points across shells."""
    work = Counter()
    sums = [0]

    def value(B):
        for j in range(len(sums), B + 1):
            rows = hybrid_shell(br, j, 1, br.height(j), work)
            sums.append((sums[-1] + sum(L * L for L in rows)) % br.n)
        return scalar_from_raw(br, B, sums[B])

    a, S, calls = first_visible(br, value)
    return dict(
        a=a,
        S=S,
        factor=split(br, a, S),
        moment_cutoffs=calls,
        shells_computed=len(sums) - 1,
        work=dict(work),
    )


def benchmark(br):
    assert br.c**2 + br.delta * br.ell**2 == 4 * br.n and gcd(br.c * br.ell, br.n) == 1
    # cap < p holds in the balanced large-input regime and on these toy cases.
    direct = [[]] + [shell_reference(br, j) for j in range(1, br.cap + 1)]
    proj_work = Counter()
    native = [[]] + [
        sorted(projected_shell(br, j, 1, br.height(j), proj_work))
        for j in range(1, br.cap + 1)
    ]
    assert native == direct
    arc_work = Counter()
    for j in range(1, br.cap + 1):
        assert arc_shell(br, j, 1, br.height(j), arc_work) == direct[j]
    raw_moment = [0]
    for vals in native[1:]:
        raw_moment.append((raw_moment[-1] + sum(L * L for L in vals)) % br.n)
    a, S, calls = first_visible(br, lambda B: scalar_from_raw(br, B, raw_moment[B]))
    g = split(br, a, S)
    streaming = native_moment_factor(br)
    assert (streaming["a"], streaming["S"], streaming["factor"]) == (a, S, g)
    expected_a = next(
        j for j in range(1, br.cap + 1) if len(direct[j]) > int(isqrt(j) ** 2 == j)
    )
    assert a == expected_a
    h = (len(direct[a]) - int(isqrt(a) ** 2 == a)) // 2
    assert h in (1, 2)
    # Reference oracle gives signed shell residues. L=0 came from the omitted axis.
    signed = (
        direct[a] + [-L for L in direct[a]] + ([0, 0] if isqrt(a) ** 2 == a else [])
    )

    @cache
    def residues(p):
        hist = Counter(v % p for v in signed)
        return tuple(hist[i] % 16 for i in range(p))

    S2, buckets, last_prime = periodic_decode(br, a, lambda p, k: residues(p)[k])
    assert S2 == S
    T = br.height(a)
    configs = {}
    for name, size in [
        ("fourth_root", max(2, isqrt(isqrt(br.n)))),
        ("sqrt_N", isqrt(br.n)),
        ("one_block", T),
    ]:
        w = Counter()
        if name == "one_block":
            # One interval with no coprimality test (projection control).
            vals = projected_shell(br, a, 1, T, w)
        else:
            vals = screened_shell(br, a, 1, T, size, w)
        assert sorted(vals) == direct[a]
        configs[name] = dict(w)
    # Queries crossing block boundaries plus +/-1 at actual hits, interior
    # windows and the full endpoint. Verify count and moment, including collisions.
    size = max(2, isqrt(isqrt(br.n)))
    endpoints = {0, 1, T // 4, T // 2, 3 * T // 4, T, max(0, size - 1), size, size + 1}
    endpoints.update(min(T, max(0, L + e)) for L in direct[a] for e in (-1, 0, 1))
    intervals = set((1, t) for t in endpoints if t <= T)
    cuts = sorted(t for t in endpoints if t <= T)
    intervals.update((l + 1, t) for l, t in zip(cuts, cuts[1:]) if l < t)
    for l, t in sorted(intervals):
        expect = direct[a][bisect_left(direct[a], l) : bisect_right(direct[a], t)]
        assert sorted(projected_shell(br, a, l, t, Counter())) == expect
    # One genuinely truncated, unaligned sieve request per branch.
    left = max(1, T // 7)
    right = max(left, T * 5 // 7)
    expect = direct[a][bisect_left(direct[a], left) : bisect_right(direct[a], right)]
    assert screened_shell(br, a, left, right, size, Counter()) == expect
    assert arc_shell(br, a, left, right, Counter()) == expect

    for B, window in ((a, (left, right)), (br.cap, None)):
        reference = [
            v
            for j in range(1, B + 1)
            for v in direct[j]
            if window is None or window[0] <= v <= window[1]
        ]
        got = disk_query(br, B, window)
        assert got["positive_count_mod8"] == len(reference) % 8
        assert (
            got["positive_second_moment_mod_N"] == sum(v * v for v in reference) % br.n
        )
        if window is None:
            assert got["signed_count_mod16"] == (2 * len(reference) + 2 * isqrt(B)) % 16

    return dict(
        N=br.n,
        ell=br.ell,
        c=br.c,
        delta=br.delta,
        a=a,
        orbits=h,
        S=S,
        factor=g,
        raw_positive_moment_mod_N=raw_moment[a],
        calls=calls,
        shell_positive_rows=direct[a],
        height=T,
        intervals=len(intervals),
        periodic_buckets=buckets,
        largest_auxiliary_prime=last_prime,
        full_shell_batches=configs,
        all_shell_projection=dict(proj_work),
        all_shell_arcs=dict(arc_work),
    )


def run(path):
    cases = [
        Branch(*r)
        for r in json.loads(Path(__file__).with_name("public_cases.json").read_text())
    ]
    records = [benchmark(br) for br in cases]
    # Deliberately selected tied minima. Public fixture data enter the decoder.
    tie_cases = [
        (3901, 7, 85, 171),
        (69689, 13, 474, 320),
        (55973, 11, 432, 308),
        (12827, 11, 123, 299),
        (7571, 5, 153, 275),
        (13589, 7, 204, 260),
        (22331, 7, 257, 475),
        (36031, 5, 357, 667),
        (27661, 7, 284, 612),
        (56153, 11, 369, 731),
    ]
    ties = [benchmark(Branch(*params)) for params in tie_cases]
    assert all(item["orbits"] == 2 for item in ties)
    sums = {}
    for name in records[0]["full_shell_batches"]:
        w = Counter()
        for r in records:
            w.update(r["full_shell_batches"][name])
        sums[name] = dict(w)
    summary = dict(
        branches=len(records),
        orbits=dict(Counter(r["orbits"] for r in records)),
        projected_intervals=sum(r["intervals"] for r in records),
        whole_disk_rows=sum(r["height"] for r in records),
        independent_shells=sum(br.cap for br in cases),
        max_moment_calls=max(len(r["calls"]) for r in records),
        largest_auxiliary_prime=max(r["largest_auxiliary_prime"] for r in records),
        full_shell_batches=sums,
    )
    path.write_text(
        json.dumps(
            dict(summary=summary, records=records, additional_tie_controls=ties),
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--output", type=Path, help="Run the benchmark and save its measurements."
    )
    group.add_argument(
        "--branch",
        nargs=4,
        type=int,
        metavar=("N", "ELL", "C", "DELTA"),
        help="Factor an eligible branch with the native arc and moment algorithm.",
    )
    parser.add_argument(
        "--cutoff",
        type=int,
        help="Evaluate the disk through this nonnegative norm quotient.",
    )
    parser.add_argument(
        "--rows",
        nargs=2,
        type=int,
        metavar=("LO", "HI"),
        help="Restrict to an inclusive interval of positive determinant rows.",
    )
    args = parser.parse_args()
    if args.cutoff is not None and (args.branch is None or args.cutoff < 0):
        parser.error("--cutoff requires --branch and a nonnegative integer.")
    if args.rows and (args.cutoff is None or not 1 <= args.rows[0] <= args.rows[1]):
        parser.error("--rows requires --cutoff and 1 <= LO <= HI.")
    if args.branch:
        br = Branch(*args.branch)
        if (
            br.n <= 1
            or br.n % 2 == 0
            or br.c <= 0
            or br.ell <= 0
            or br.delta <= 4
            or br.delta % 4 not in (0, 3)
            or br.c**2 + br.delta * br.ell**2 != 4 * br.n
            or gcd(br.n, br.ell * br.delta) != 1
            or (br.c - (br.delta % 2) * br.ell) % 2
        ):
            parser.error("The public data fail the branch eligibility checks.")
        answer = (
            native_moment_factor(br)
            if args.cutoff is None
            else disk_query(br, args.cutoff, args.rows)
        )
        print(json.dumps(answer, indent=2))
    else:
        run(args.output)
