"""Split odd inputs using a fixed window and a fixed signed disk cutoff."""

from math import gcd, isqrt


def query_forms(n, start):
    half = (n - 1) // 2
    radius = half // 2
    origin = start + radius + 1
    leading = (n + 1)**4
    cross = (2 * origin) % n
    cross += n * ((cross - 1) % (n + 1))
    padding = n * (n + 1)
    base = leading * (radius + 1)**2
    last = base + (origin**2 - 1 - base) % n
    cutoff = (base + leading * radius**2 + padding * radius + 2*n - 1) // n
    first = (leading, cross, last)
    second = (leading, cross, last + 2*n*cutoff)
    return first, second, cutoff


def mixed_root(n, value):
    value %= n
    return value not in (1, n - 1) and (value * value - 1) % n == 0


def window_parity(n, start, disk_count):
    half = (n - 1) // 2
    first, second, cutoff = query_forms(n, start)
    residue = (disk_count(n, first, cutoff) - disk_count(n, second, cutoff)) % 4
    assert residue in (0, 2)
    parity = residue // 2
    if half % 2 == 0:
        endpoint = (start + half + 1) % n
        parity ^= int((endpoint * endpoint - 1) % n == 0)
    for root in (1, n - 1):
        parity ^= int(1 <= (root - start) % n <= half)
    return parity


def split_from_windows(n, parity):
    """Return a split if n has at least two distinct prime divisors."""
    half = (n - 1) // 2
    low, high = 0, half // 2
    if half % 2 and mixed_root(n, high + 1):
        return gcd(high, n)
    if parity(low) == 0:
        return None
    while high - low > 1:
        mid = (high + low) // 2
        if parity(mid):
            low = mid
        else:
            high = mid
    for root in (high, high + half):
        if mixed_root(n, root):
            divisor = gcd(root - 1, n)
            assert 1 < divisor < n
            return divisor
    raise AssertionError("Window parity changed without a root")


def find_factor(n, disk_count):
    return split_from_windows(n, lambda start: window_parity(n, start, disk_count))

def reference_count(n, form, cutoff):
    a, b, c = form
    budget = n * cutoff
    discriminant = 4*a*c - b*b
    count = 0
    for y in range(-isqrt(4*a*budget // discriminant), isqrt(4*a*budget // discriminant)+1):
        radius = isqrt(4*a*budget - discriminant*y*y)
        low = -((radius + b*y) // (2*a))
        high = (radius - b*y) // (2*a)
        for x in range(low, high + 1):
            value = a*x*x + b*x*y + c*y*y
            if 0 < value <= budget and value % n == 0:
                count += 1
    return count


if __name__ == "__main__":
    splits = prime_powers = windows = 0
    for n in range(5, 321, 2):
        half = (n - 1) // 2
        cutoffs = set()
        for start in range(half // 2 + 1):
            first, second, cutoff = query_forms(n, start)
            cutoffs.add(cutoff)
            for a, b, c in (first, second):
                assert gcd(gcd(a, b), c) == 1
                assert 3*cutoff**2 < 4*a*c - b*b
            expected = sum(mixed_root(n, x) for x in range(start+1, start+half+1)) % 2
            assert window_parity(n, start, reference_count) == expected
            windows += 1
        assert len(cutoffs) == 1
        divisor = find_factor(n, reference_count)
        has_mixed = any(mixed_root(n, x) for x in range(n))
        if has_mixed:
            assert divisor and 1 < divisor < n and n % divisor == 0
            splits += 1
        else:
            assert divisor is None
            prime_powers += 1
    print(f"{windows} windows, {splits} splits, {prime_powers} primes or prime powers")
