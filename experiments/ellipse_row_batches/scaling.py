"""Compare two series of complete exact first shells, with public cutoff supplied."""

from collections import Counter
from math import isqrt
from pathlib import Path
from time import perf_counter
import argparse, json
from row_batches import Branch, projected_shell, arc_shell, split

# Inputs (N,ell,c,Delta,j). The growing-discriminant cutoffs were selected
# using known fixture factors. Timings measure a single public shell only.
CASES = [
    (85119121, 3, 18452, 20, 1),
    (27380982829, 3, 330944, 20, 1),
    (7241728574161, 3, 5382092, 20, 1),
    (304115721499441, 3, 34877828, 20, 1),
    (94165439, 5, 19384, 36892, 73),
    (24057426803, 7, 310194, 191624, 219),
    (6157336418887, 3, 4962793, 3479411, 339),
    (394065987449677, 7, 39702156, 56301028, 3173),
]


def run():
    rows = []
    for args in CASES:
        br = Branch(*args[:4])
        j = args[4]
        T = br.height(j)
        runs = {}
        points = None
        for key, method in [("project", projected_shell), ("arc", arc_shell)]:
            work = Counter()
            start = perf_counter()
            result = sorted(method(br, j, 1, T, work))
            seconds = perf_counter() - start
            if points is not None:
                assert points == result
            points = result
            runs[key] = {"seconds": seconds, **work}
        k = isqrt(j)
        correction = (k * br.c * br.ell) ** 2 if k * k == j else 0
        S = (sum(L * L for L in points) - correction) * pow(br.c, -2, br.n) % br.n
        divisor = split(br, j, S)
        # Check a clipped endpoint immediately before a point and at the point.
        for t in (points[-1] - 1, points[-1]):
            assert arc_shell(br, j, 1, t, Counter()) == [L for L in points if L <= t]
        rows.append(
            dict(
                N=br.n,
                bits=br.n.bit_length(),
                ell=br.ell,
                c=br.c,
                delta=br.delta,
                j=j,
                positive_rows=points,
                divisor=divisor,
                methods=runs,
            )
        )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    out = run()
    p.parse_args().output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
