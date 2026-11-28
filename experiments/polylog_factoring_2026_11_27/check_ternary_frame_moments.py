"""Read three ternary Gram moments by a one-coordinate exact square scan."""

from math import gcd, isqrt

from sympy import jacobi_symbol

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
else:
    from check_mixed_ternary_embeddings import _all_frames


def row_roots(n, x):
    """Yield (a,b,c) with b,c positive and b²+(na+xb)²+c²=n²."""
    for b in range(1, n):
        residue = x*b % n
        for t in (residue, residue - n):
            remainder = n*n - b*b - t*t
            if remainder <= 0:
                continue
            c = isqrt(remainder)
            if c*c == remainder:
                yield (t - x*b) // n, b, c


def frame_moments(n, x):
    """Return the three normalized quartic traces on a promised mixed semiprime.

    Each nonzero coordinate of a representative noncoordinate frame contributes
    once. Only n and x enter this deterministic O(n) row scan.
    """
    if n < 3 or n % 2 == 0 or jacobi_symbol(-(x*x + 1), n) != -1:
        raise ValueError("expected an odd semiprime with a mixed Jacobi symbol")
    sums = [0, 0, 0]
    for a, b, _ in row_roots(n, x):
        sums[0] += b**4
        sums[1] += a*b**3
        sums[2] += a*a*b*b
    return tuple(sums)


def brute_traces(n, x):
    total = [0, 0, 0]
    for u, v in _all_frames(n, x):
        for a, b in zip(u, v):
            total[0] += b**4
            total[1] += a*b**3
            total[2] += a*a*b*b
    coordinate = n**4, -x*n**3, x*x*n*n
    residual = [value - 24*known for value, known in zip(total, coordinate)]
    assert all(value % 48 == 0 for value in residual)
    return tuple(value // 48 for value in residual)


def cm_coefficient_square(n):
    """Weight-three CM coefficient at n² via generator squares in Z[sqrt(-2)]."""
    total = 0
    limit = isqrt(n*n // 2)
    for z in range(-limit, limit + 1):
        rem = n*n - 2*z*z
        y = isqrt(rem)
        if y*y == rem:
            total += (2 if y else 1)*(y*y - 2*z*z)
    return total // 2  # two generators per ideal


def run():
    checked = compared = row_count = 0
    fourth_missed = first_two_missed = cm_checks = 0
    for p, q in ((3, 5), (3, 7), (3, 11), (3, 13), (5, 7), (5, 11), (7, 11), (11, 13)):
        n = p*q
        mixed = [x for x in range(n) if jacobi_symbol(-(x*x + 1), n) == -1]
        for i, x in enumerate(mixed):
            checked += 1
            r, d = (p, q) if jacobi_symbol(-(x*x + 1), p) == 1 else (q, p)
            rows = list(row_roots(n, x))
            row_count += len(rows)
            assert len(rows) == (2 if x % r == 0 else 3)
            assert all(gcd(b, n) == d for _, b, _ in rows)
            moments = frame_moments(n, x)
            assert gcd(n, *moments) == d
            fourth_missed += gcd(n, moments[0]) != d
            if gcd(n, *moments[:2]) != d:
                assert r == 3
                first_two_missed += 1
            if x == 1:
                coefficient = cm_coefficient_square(n)
                assert coefficient == 8*moments[2] - 5*n*n
                assert gcd(n, coefficient) == d
                cm_checks += 1
            if i in (0, len(mixed)//2, len(mixed)-1):
                assert moments == brute_traces(n, x)
                compared += 1
    print(f"{checked} mixed Grams, {row_count} positive rows, "
          f"{compared} independent full-trace comparisons")
    print(f"First quartic missed {fourth_missed}; first two missed {first_two_missed}")
    print("N=77, x=1:", list(row_roots(77, 1)), frame_moments(77, 1))
    print(f"{cm_checks} CM comparisons; coefficient at 77²:", cm_coefficient_square(77))


if __name__ == "__main__":
    run()
