"""Check the fixed-cutoff CRT separation bound with factor-supplied controls."""

from collections import Counter, defaultdict
from math import gcd, isqrt


def modular_rank(columns, prime):
    basis = {}
    for column in columns:
        vector = column[:]
        for i, value in enumerate(vector):
            if not value:
                continue
            if i in basis:
                vector = [
                    (a - value * b) % prime
                    for a, b in zip(vector, basis[i], strict=True)
                ]
            else:
                inverse = pow(value, -1, prime)
                basis[i] = [entry * inverse % prime for entry in vector]
                break
    return len(basis)


def column_classes(p, q, length):
    shifts = defaultdict(set)
    inverse = pow(q, -1, p)
    for j in range(q):
        size = (length + q - 1 - j) // q
        shifts[size].add(j * inverse % p)
    return {size: len(indices) for size, indices in shifts.items()}


def exhaustive():
    cases = 0
    for p, q in ((3, 5), (5, 7), (7, 11), (11, 13), (13, 19)):
        for length in range(q, (p - 1) * q + 1):
            columns = [[0] * p for _ in range(q)]
            for s in range(length):
                columns[s % q][s % p] = 1
            classes = column_classes(p, q, length)
            assert all(0 < size < p for size in classes)
            bound = max(classes.values())
            assert bound >= min(p, (q + 1) // 2)
            assert modular_rank(columns, p) >= bound
            for first in (1, 2):
                weighted = [
                    [entry * ((first + 2 * i) % p) for i, entry in enumerate(col)]
                    for col in columns
                ]
                assert modular_rank(weighted, p) >= bound - 1
            cases += 1
    print(f"{cases} intervals checked by modular elimination")


def check_cap(n, p, ell, c, delta):
    q = n // p
    p, q = sorted((p, q))
    assert p * q == n and c * c + delta * ell * ell == 4 * n
    cap = isqrt(delta // 3)
    maximum = isqrt(4 * n * cap - delta)
    first = 1 if delta % 2 else 2
    length = (maximum - first) // 2 + 1
    assert q <= length <= (p - 1) * q
    bound = max(column_classes(p, q, length).values()) - 1
    trace, matches = 0, Counter()
    for shell in range(1, cap + 1):
        height = 4 * n * shell
        for y in range(1, isqrt((height - 1) // delta) + 1):
            square = height - delta * y * y
            w = isqrt(square)
            if w * w == square:
                trace += w
                matches[y] += 1
    normalized = (trace - c * isqrt(cap) * (isqrt(cap) + 1) // 2)
    normalized = normalized * ell * pow(c, -1, n) % n
    square = normalized * normalized % n
    t = isqrt(square)
    assert t * t == square
    divisor = gcd(normalized, n)
    if divisor == 1:
        divisor = gcd(normalized - t, n)
    assert divisor in (p, q)
    print(n, cap, length, bound, matches[1], trace % n, divisor)


if __name__ == "__main__":
    exhaustive()
    print("N cap row_length weighted_rank_lower_bound row_matches trace divisor")
    check_cap(11009, 109, 11, 204, 20)
    check_cap(28829, 227, 11, 336, 20)
    check_cap(38942974393, 153247, 1, 394679, 384531)
