"""End-to-end ellipse factorer for a promised balanced product of two distinct odd primes."""

from argparse import ArgumentParser
from bisect import bisect_right
from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ellipse_row_batches"))
from row_batches import Branch, first_visible, native_moment_factor, split


def branch_from_n(n):
    """Use ell=1 for the weighted moment. A nonunit c gives a proper factor directly."""
    if n < 15 or n % 2 == 0 or isqrt(n) ** 2 == n:
        raise ValueError("expected a balanced product of two distinct odd primes")
    c = isqrt(4 * n)
    delta = 4 * n - c * c
    if delta <= 4:
        c -= 1
        delta = 4 * n - c * c
    g = gcd(c, n)
    if g > 1:
        return None, g
    return Branch(n, 1, c, delta), None


@dataclass
class Disk:
    indices: list
    counts: list
    moments: list

    def prefix(self, B):
        i = bisect_right(self.indices, B) - 1
        return self.counts[i], self.moments[i]


def positive_matches(br, B, work):
    """Yield positive (w,y) in one disk. B<=cap<min(p,q) is the promised range.

    Coordinates w=|2x+b*y|, y=|y| are strictly positive on nonzero hits.
    Each positive pair represents one four-element conjugation/sign orbit.
    """
    if not 0 <= B <= br.cap:
        raise ValueError("cutoff outside the first-minimum search range")
    if B >= isqrt(br.n // 2) + 1:
        raise ValueError("cutoff exceeds the public balanced-semiprime bound")
    mod = 4 * br.n
    radius = mod * B
    X = isqrt(radius)
    Y = isqrt(radius // br.delta)
    table = {}
    for y in range(1, Y + 1):
        dy2 = br.delta * y * y
        table.setdefault((-dy2) % mod, []).append((dy2, y))
    work["y_insertions"] += Y
    work["w_probes"] += X
    work["builds"] += 1
    work["peak_y_entries"] = max(work["peak_y_entries"], Y)
    # Four square roots mod N at most for any target; each w occurs at most
    # ceil(X/N) times per residue. No per-probe scan over the whole short list.
    for w in range(1, X + 1):
        ww = w * w
        bucket = table.get(ww % mod)
        if not bucket:
            continue
        while bucket and bucket[-1][0] + ww > radius:
            bucket.pop()
        for _, y in bucket:
            work["hits"] += 1
            yield w, y


def residue_disk(br, B, work):
    """Build all prefix counts and second moments through B from one join."""
    events = {}
    for w, y in positive_matches(br, B, work):
        j = (w * w + br.delta * y * y) // (4 * br.n)
        count, weight = events.get(j, (0, 0))
        events[j] = (count + 1, weight + y * y)
    indices = [0]
    counts = [0]
    moments = [0]
    for j, (count, weight) in sorted(events.items()):
        indices.append(j)
        counts.append(counts[-1] + count)
        moments.append((moments[-1] + weight) % br.n)
    return Disk(indices, counts, moments)


def join_factor(br):
    capacity = -1
    disk = None
    work = Counter()

    def value(B):
        nonlocal capacity, disk
        if capacity < B:
            capacity = min(br.cap, max(B, 2 * capacity))
            disk = residue_disk(br, capacity, work)
        _, moment = disk.prefix(B)
        k = isqrt(B)
        return (moment - br.ell**2 * (k * (k + 1) * (2 * k + 1) // 6)) % br.n

    a, S, calls = first_visible(br, value)
    return dict(
        a=a,
        S=S,
        factor=split(br, a, S),
        moment_cutoffs=calls,
        maximum_cutoff_built=capacity,
        work=dict(work),
    )


def factor(n, method="shell"):
    if method not in ("shell", "join"):
        raise ValueError("unknown evaluator")
    br, direct = branch_from_n(n)
    if direct is not None:
        result = {"factor": direct, "route": "branch_gcd", "work": {}}
    else:
        result = native_moment_factor(br) if method == "shell" else join_factor(br)
        result.update(
            route=method, ell=br.ell, c=br.c, delta=br.delta, search_cap=br.cap
        )
    d = result.pop("factor")
    if not (1 < d < n and n % d == 0):
        raise ArithmeticError("factor failed exact verification")
    return dict(N=n, factors=sorted((d, n // d)), **result)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("--method", choices=("shell", "join"), default="shell")
    args = parser.parse_args()
    print(json.dumps(factor(args.N, args.method), indent=2))
