"""Check the signed ordinate bound independently on known-factor fixtures."""

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path

from benchmark import odd_primes, trial
from factor import Branch, branch_from_n, positive_matches
from weighted_recurrence_single import decode_trace


def fixtures():
    for p, q in combinations(odd_primes(503), 2):
        if q <= 2 * p:
            br, d = branch_from_n(p * q)
            if br is None:
                assert d in (p, q)
            else:
                yield br, p
    records = Path(__file__).resolve().parents[1] / "ellipse_row_batches/results.json"
    data = json.loads(records.read_text())
    for record in data["records"] + data["additional_tie_controls"]:
        n = record["N"]
        yield Branch(n, record["ell"], record["c"], record["delta"]), trial(n)[0]


def exact_shells(br):
    for j in range(1, br.cap + 1):
        height = 4 * br.n * j
        for y in range(1, isqrt(height // br.delta) + 1):
            w2 = height - br.delta * y * y
            w = isqrt(w2)
            if w > 0 and w * w == w2:
                yield w, y


def audit(br, p):
    # Only the fixture constructs the mixed root. The decoder gets the raw trace.
    n, q = br.n, br.n // p
    public = br.c * pow(br.ell, -1, n) % n
    mixed = (public + p * ((-2 * public * pow(p, -1, q)) % q)) % n
    points = sorted(exact_shells(br))
    assert points == sorted(positive_matches(br, br.cap, Counter()))
    raw = sum(w for w, _ in points) % n
    d, _ = decode_trace(br, raw)
    assert d in (p, q)
    signs = set()
    signed_y = trace = count = 0
    for w, y in points:
        if br.ell * w == br.c * y:
            continue
        assert gcd(y, n) == 1
        sign = 1 if (w - mixed * y) % n == 0 else -1
        assert (w - sign * mixed * y) % n == 0
        signs.add(sign)
        signed_y += sign * y
        trace += w
        count += 1
    normalized = trace * pow(public, -1, n) % n
    assert 0 < signed_y**2 < n
    assert normalized**2 % n == signed_y**2
    return Fraction(signed_y**2, n), count == 2 and len(signs) == 2


def main():
    maximum = Fraction(0)
    total = opposite = 0
    for br, p in fixtures():
        ratio, opposed_pair = audit(br, p)
        maximum = max(maximum, ratio)
        total += 1
        opposite += opposed_pair
    assert (total, opposite, maximum) == (2128, 100, Fraction(21609, 31313))
    print(f"{total} branches, {opposite} opposite pairs, maximum D²/N = {maximum}")


if __name__ == "__main__":
    main()
