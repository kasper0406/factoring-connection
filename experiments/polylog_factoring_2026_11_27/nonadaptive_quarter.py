"""Split inputs with two distinct odd prime divisors by nonadaptive disk queries."""

from math import gcd, isqrt

from fixed_window_split import reference_count


def queries(n):
    radius = n // 4
    leading = (n + 1)**3
    base = leading * (radius + 1)**2
    padding = n * (n + 1)
    cutoff = (base + leading*radius**2 + padding + n - 1) // n
    reference = (leading, 0, n*cutoff + (1 - n*cutoff) % (n + 1) + n + 1)
    forms = []
    scale = 1
    for _ in range(n.bit_length() - 1):
        last = base + (-scale*scale-base) % n
        last += n * ((last - 1) % (n + 1))
        forms.append((scale, (leading, 0, last)))
        scale = 2 * scale % n
    return cutoff, reference, forms


def find_factor(n, disk_count):
    cutoff, reference, forms = queries(n)
    baseline = disk_count(n, reference, cutoff)
    bit, prefix = 0, 0
    for scale, form in forms:
        residue = (disk_count(n, form, cutoff) - baseline) % 8
        assert residue in (0, 4)
        known = int(min(scale, n - scale) <= n // 4)
        mixed_inside = (residue // 4) ^ known
        bit ^= 1 ^ mixed_inside
        prefix = 2*prefix + bit
    denominator = 1 << n.bit_length()
    root = (n*prefix + denominator - 1) // denominator
    assert (root*root - 1) % n == 0
    divisor = gcd(root - 1, n)
    assert 1 < divisor < n
    return divisor


if __name__ == "__main__":
    primes = [p for p in range(3, 80, 2) if all(p%d for d in range(3, isqrt(p)+1, 2))]
    cases = [(p, q) for p in primes for q in primes if p < q]
    cases += [(9, 5), (9, 25), (27, 49), (25, 7), (9, 121)]
    count = 0
    for p, q in cases:
        n = p*q
        cutoff, reference, forms = queries(n)
        for a, _, c in [reference] + [form for _, form in forms]:
            assert gcd(a, c) == 1 and 3*cutoff**2 < 4*a*c
        assert find_factor(n, reference_count) in (p, q)
        count += len(forms) + 1
    print(f"{len(cases)} splits using {count} signed disk queries")
