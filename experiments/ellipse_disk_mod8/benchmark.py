"""Whole-disk modulo-eight checks, order comparison, and GF(2) rectangle-rank witnesses."""

from argparse import ArgumentParser
from math import isqrt
from pathlib import Path
from time import perf_counter
import json
from decode import branch_from_n, decode
from disk import AdaptiveDisk, DiskCounter, FixedDisks
from row_batches import Branch


def independent(br):
    # Direct positive-quadrant area scan; never given to the symbolic evaluator.
    t = 4 * br.n * br.cap
    hits = []
    for y in range(1, isqrt(t // br.delta) + 1):
        for w in range(1, isqrt(t - br.delta * y * y) + 1):
            j, r = divmod(w * w + br.delta * y * y, 4 * br.n)
            if r == 0:
                hits.append((j, y))
    return hits


def check(br, order="lsb", all_cutoffs=True, early=True):
    hits = independent(br)
    start = perf_counter()
    counter = DiskCounter(br, order, early=early)
    built = perf_counter()
    result = decode(br, counter.count)
    used = counter.metrics()
    elapsed = perf_counter() - start

    parity_counter = DiskCounter(br, order, early=True, parity=True)
    parity_result = decode(br, parity_counter.count, parity=True)
    assert (
        parity_result["a"] == result["a"]
        and parity_result["factors"] == result["factors"]
    )
    parity_metrics = parity_counter.metrics()
    counts = 0
    for B in range(br.cap + 1) if all_cutoffs else [br.cap]:
        good = [y for j, y in hits if j <= B]
        for bit in (None, *range(len(counter.y))):
            exact = len(good) if bit is None else sum((y >> bit) & 1 for y in good)
            assert counter.count(B, bit) == exact % 8
            assert parity_counter.count(B, bit) == exact % 2
            counts += 1
    return dict(
        N=br.n,
        ell=br.ell,
        c=br.c,
        delta=br.delta,
        a=result["a"],
        mixed_orbits=result["mixed_orbits"],
        y_sum=result["y_sum"],
        y_xor=parity_result["y_xor"],
        factors=result["factors"],
        verified_queries=counts,
        work=used,
        parity_work=parity_metrics,
        seconds=elapsed,
        build_seconds=built - start,
    )


def rank_controls(n, d):
    # Positive rectangle, same congruence modulo N; no ellipse or cutoff constraint.
    groups = {}
    for y in range(1, n + 1):
        groups.setdefault((-d * y * y) % n, y)
    squares = {}
    for w in range(1, n + 1):
        squares.setdefault((w * w) % n, w)
    assert groups.keys() == squares.keys()
    # Pair one row and one column at each of the keys. Equality is an identity matrix.
    return len(groups)


def timed_factor(n, early, parity=False):
    start = perf_counter()
    br, d = branch_from_n(n)
    assert d is None
    counter = AdaptiveDisk(br, early=early, parity=parity)
    result = decode(br, counter.count, parity=parity)
    return dict(**result, work=counter.metrics(), seconds=perf_counter() - start)


def run():
    base = Path(__file__).resolve().parents[1]
    cases = json.loads((base / "ellipse_row_batches" / "public_cases.json").read_text())
    rows = [check(Branch(*case)) for case in cases]
    ties = [
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
    extra = [check(Branch(*case)) for case in ties]
    assert all(r["mixed_orbits"] == 2 for r in extra)
    # Area scans use no factors; the fixture controls only set the sampling order.
    scales = json.loads((base / "ellipse_aggregate" / "results.json").read_text())[
        "scaling"
    ]
    ns = (
        [11009, 28829]
        + [s["N"] for s in scales if s["bits"] <= 20]
        + [next(s["N"] for s in scales if s["bits"] == 23)]
    )
    scaling = []
    for n in ns:
        plain = timed_factor(n, False)
        early = timed_factor(n, True)
        assert plain["a"] == early["a"] and plain["y_sum"] == early["y_sum"]
        print("scaling", n, early["work"], round(early["seconds"], 3), flush=True)
        par = timed_factor(n, True, True)
        assert par["a"] == early["a"] and par["factors"] == early["factors"]
        scaling.append(dict(N=n, plain=plain, early=early, parity=par))
    orders = {}
    for n in [11009, 28829, 438751]:
        br, _ = branch_from_n(n)
        orders[n] = {}
        for key in ("lsb", "msb", "grouped"):
            start = perf_counter()
            counter = AdaptiveDisk(br, key)
            ans = decode(br, counter.count)
            orders[n][key] = dict(
                **counter.metrics(), a=ans["a"], seconds=perf_counter() - start
            )
    both = []
    for n in ns:
        start = perf_counter()
        br, d = branch_from_n(n)
        counter = FixedDisks(br)
        r = decode(br, counter.count, parity=True)
        assert r["factors"][0] * r["factors"][1] == n
        both.append(
            dict(
                N=n,
                a=r["a"],
                factors=r["factors"],
                work=counter.metrics(),
                seconds=perf_counter() - start,
            )
        )
    ranks = []
    for p, q in [(5, 7), (11, 13), (31, 43), (101, 109)]:
        n = p * q
        d = 4 * n - 1
        rank = rank_controls(n, d)
        assert rank == (p + 1) * (q + 1) // 4
        ranks.append(dict(N=n, rank_mod2=rank))

    # Exercise a real wrap in an enlarged disk below p, including nonaxis
    # principal points, outside the decoder cap. Only count equality is claimed.
    class LargerDisk:
        def __init__(self, br, cap):
            self.n, self.c, self.ell, self.delta, self.cap = (
                br.n,
                br.c,
                br.ell,
                br.delta,
                cap,
            )

    wrap = []
    for limit in (16, 32, 64):
        br = LargerDisk(Branch(*cases[0]), limit)
        pts = independent(br)
        counter = DiskCounter(br, early=True)
        got = counter.count(limit)
        assert got == len(pts) % 8 and len(pts) >= 8
        wrap.append(dict(N=br.n, B=limit, exact=len(pts), mod8=got))
    return dict(
        records=rows,
        tie_controls=extra,
        scaling=scaling,
        orders=orders,
        ranks=ranks,
        quotient_projection=both,
        wrapped=wrap,
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run()
    args.output.write_text(json.dumps(out, indent=2) + "\n")
