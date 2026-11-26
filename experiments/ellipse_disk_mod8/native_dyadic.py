"""Compute the root-prefix parity by floor sums and local quadratic intersections."""

from argparse import ArgumentParser
from math import gcd, isqrt


def one_form_split(n, disk):
    """Split any odd composite with a mod-8 disk oracle, returning None on primes."""
    if n < 3 or n % 2 == 0:
        raise ValueError("expected an odd integer >= 3")
    if n == 3:
        return None
    middle = (n - 1) // 2
    a = n + 1
    d = a * n * ((middle**2 + n - 1) // n) - 1

    def probe(radius, row_offset):
        cutoff = (a * radius**2 + row_offset + n - 1) // n
        return disk(n, a, d, cutoff) % 8

    # The first nonempty axis interval contains exactly one positive multiple.
    lo, hi = 1, 1
    while hi < middle:
        lo, hi = hi, min(2 * hi, middle)
        if probe(hi, 0):
            while hi - lo > 1:
                mid = (hi + lo) // 2
                if probe(mid, 0):
                    hi = mid
                else:
                    lo = mid
            if n % hi:
                raise ArithmeticError("invalid axis divisor")
            return hi

    # Now N is squarefree. Its only possible axis point here is the origin.
    if probe(middle, d) == 4:
        return None
    lo, hi = 1, middle
    # The prefix through lo always has parity 1, and through hi parity 0.
    while hi - lo > 1:
        mid = (hi + lo) // 2
        mid_parity = probe(mid, d) // 4
        if mid_parity:
            lo = mid
        else:
            hi = mid
    divisor = gcd(hi - 1, n)
    if not (1 < divisor < n and hi * hi % n == 1):
        raise ArithmeticError("invalid root from the disk oracle")
    return divisor


def floor_sum(count, modulus, slope, offset):
    """Sum floor((slope*j+offset)/modulus) for 0 <= j < count."""
    total = 0
    while count:
        quotient, slope = divmod(slope, modulus)
        total += quotient * count * (count - 1) // 2
        quotient, offset = divmod(offset, modulus)
        total += quotient * count
        height = slope * count + offset
        if height < modulus:
            break
        count, offset = divmod(height, modulus)
        modulus, slope = slope, modulus
    return total


def block_count(n, start, length):
    """Count roots in [start, start+length) using one enclosing parallelogram."""
    slope = 2 * start
    offset = start * start - 1
    width = (length - 1) ** 2
    assert width < n and 2 * (length - 1) * width < n

    def candidates(size):
        upper = floor_sum(size, n, slope, offset + width)
        return upper - floor_sum(size, n, slope, offset - 1)

    def select(rank):
        lo, hi = 0, length
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if candidates(mid) >= rank:
                hi = mid
            else:
                lo = mid
        index = hi - 1
        return start + index, (slope * index + offset + width) // n

    count = candidates(length)
    if not count:
        return 0
    x1, y1 = select(1)
    if count == 1:
        return int((x1 * x1 - 1) == n * y1)
    x2, y2 = select(2)
    dx, dy = x2 - x1, y2 - y1
    # Substitute the unique candidate line into N*y=x*x-1.
    a, b, c = dx, -n * dy, n * (dy * x1 - dx * y1) - dx
    disc = b * b - 4 * a * c
    if disc < 0:
        return 0
    radical = isqrt(disc)
    if radical * radical != disc:
        return 0
    hits = 0
    for numerator in {-b - radical, -b + radical}:
        x, rem = divmod(numerator, 2 * a)
        if rem == 0 and start <= x < start + length and (x * x - 1) % n == 0:
            hits += 1
    return hits


def block_size(n):
    length = 1
    while 2 * (2 * length - 1) ** 3 < n:
        length *= 2
    return length


def prefix(n, cutoff):
    """Native parity of #{1 <= x <= cutoff: x*x = 1 mod n}."""
    if n < 3 or n % 2 == 0 or not 0 <= cutoff <= n:
        raise ValueError("expected odd N >= 3 and 0 <= cutoff <= N")
    length = block_size(n)
    return (
        sum(
            block_count(n, start, min(length, cutoff + 1 - start))
            for start in range(1, cutoff + 1, length)
        )
        % 2
    )


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("T", type=int)
    args = parser.parse_args()
    length = block_size(args.N)
    print(
        f"N={args.N} T={args.T} parity={prefix(args.N, args.T)} "
        f"blocks={(args.T + length - 1) // length} block_length={length}"
    )
