"""Check exact threshold positions for the shifted mixed-root window."""

from itertools import combinations
from math import isqrt

from fixed_window_split import mixed_root


def check(p, q):
    n = p * q
    half = (n - 1) // 2
    limit = half // 2
    root = 1 + p * (-2 * pow(p, -1, q) % q)
    u = min(root, n-root)
    threshold = min(u, half+1-u)
    assert mixed_root(n, u)
    assert u >= p-1
    assert 2*threshold >= p-1
    assert 4*(limit-threshold) >= p-9
    for start in {
        0, limit, max(0, threshold-1), min(threshold, limit)
    }:
        observed = sum(start < x <= start+half for x in (u, n-u)) % 2
        assert observed == int(start < threshold)


if __name__ == "__main__":
    primes = [p for p in range(3, 1000, 2) if all(p % d for d in range(3, isqrt(p)+1, 2))]
    count = 0
    for p, q in combinations(primes, 2):
        check(p, q)
        count += 1
    print(f"{count} semiprimes checked")
