"""Check middle-cut degrees and ranks on squarefree odd moduli."""

from collections import defaultdict
from math import isqrt, prod


def roots(primes):
    n = prod(primes)
    values = [0]
    for prime in primes:
        cofactor = n // prime
        basis = cofactor * pow(cofactor, -1, prime)
        values = [(value + sign * basis) % n for value in values for sign in (-1, 1)]
    return n, values


def induced_matching(rows):
    rows = {row: set(columns) for row, columns in rows.items()}
    columns = defaultdict(set)
    for row, neighbors in rows.items():
        for column in neighbors:
            columns[column].add(row)
    size = 0
    while rows:
        row = next(iter(rows))
        column = next(iter(rows[row]))
        blocked_rows, blocked_columns = list(columns[column]), list(rows[row])
        size += 1
        for blocked in blocked_rows:
            for neighbor in rows.pop(blocked):
                columns[neighbor].remove(blocked)
        for blocked in blocked_columns:
            for neighbor in list(columns[blocked]):
                rows[neighbor].remove(blocked)
                if not rows[neighbor]:
                    del rows[neighbor]
            columns[blocked].clear()
    return size


def rank_binary(rows):
    pivots = {}
    for row in rows:
        while row:
            lead = row.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = row
                break
            row ^= pivots[lead]
    return len(pivots)


def check(primes):
    n, values = roots(primes)
    width = 1 << (isqrt(n).bit_length() - 1)
    midpoint = (n - 1) // 2
    limit = ((midpoint + 1) // width) * width
    graph = defaultdict(set)
    for value in values:
        if value < limit:
            quotient, remainder = divmod(value, width)
            graph[remainder].add(quotient)
    row_degree = max(map(len, graph.values()))
    reverse = defaultdict(int)
    for columns in graph.values():
        for col in columns:
            reverse[col] += 1
    column_degree = max(reverse.values())
    matching = induced_matching(graph)
    assert sum(reverse.values()) == sum(value < limit for value in values)
    assert matching * 2 * max(row_degree, column_degree)**2 >= sum(reverse.values())
    if n < 1_000_000:
        root_rows, prefix_rows = [0] * width, [0] * width
        parity = 0
        for x in range(limit):
            col, row = divmod(x, width)
            hit = x*x % n == 1
            parity ^= hit
            root_rows[row] |= int(hit) << col
            prefix_rows[row] |= parity << col
        assert rank_binary(root_rows) >= matching
        assert rank_binary(prefix_rows) >= rank_binary(root_rows) - 1
    print(n.bit_length(), len(primes), row_degree, column_degree, matching)


if __name__ == "__main__":
    primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    print("bits factors max_row_degree max_column_degree matching")
    for k in range(2, len(primes) + 1):
        check(primes[:k])
