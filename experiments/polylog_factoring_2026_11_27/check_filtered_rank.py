"""Check the support and weights of the CRT matrix after filtering."""

from itertools import combinations
from math import gcd, isqrt

from check import modular_rank


def mixed_points(n, c, delta, cap):
    for level in range(1, cap + 1):
        height = 4 * n * level
        for y in range(1, isqrt((height - 1) // delta) + 1):
            square = height - delta * y * y
            w = isqrt(square)
            if w * w == square and w != c * y:
                yield w, y


def check(p, q):
    n = p * q
    c = isqrt(4 * n)
    if 4 * n - c * c <= 4:
        c -= 1
    delta = 4 * n - c * c
    cap = isqrt(delta // 3)
    if cap < 3 or gcd(c, n) != 1:
        return None
    entries = list(mixed_points(n, c, delta, cap))
    rows = {(w % p, y % p) for w, y in entries}
    columns = {(w % q, y % q) for w, y in entries}
    assert len(entries) == len(rows) == len(columns)
    assert all(w % p and w % q and 0 < y < p for w, y in entries)
    size = len(entries)
    compressed = [[0] * size for _ in range(size)]
    row_index = {pair: i for i, pair in enumerate(sorted(rows))}
    column_index = {pair: i for i, pair in enumerate(sorted(columns))}
    for w, y in entries:
        compressed[column_index[w % q, y % q]][row_index[w % p, y % p]] = w % p
    assert modular_rank(compressed, p) == size
    d = (q - p) // 2
    m = (p + q) // 2
    if c == 2 * m - 1 and m > 2 * d * d:
        w, y = 2 * m - 4 * d * d, 2 * d
        assert {(k * w, k * y) for k in range(1, isqrt(cap) + 1)} <= set(entries)
        assert size == isqrt(cap) if 4 * cap < delta else size >= isqrt(cap)
    return cap, size


if __name__ == "__main__":
    primes = [
        k for k in range(3, 242, 2)
        if all(k % d for d in range(3, isqrt(k) + 1, 2))
    ]
    checked = 0
    for p, q in combinations(primes, 2):
        if q <= 2 * p and check(p, q) is not None:
            checked += 1
    print(f"{checked} small branches")
    for p, q in [(101, 109), (1009, 1013), (10007, 10009), (153247, 254119)]:
        cap, rank = check(p, q)
        print(f"N={p*q}, cap={cap}, exact rank={rank}")
