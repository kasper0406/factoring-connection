"""Evaluate mixed-root window parity using bounded univariate lattices."""

from collections import Counter
from math import isqrt

from sympy import Matrix, Poly, Symbol

from fixed_window_split import mixed_root, split_from_windows


def multiplicity(n, limit):
    for m in range(1, n.bit_length()+1):
        dimension = 2*m+1
        if dimension**dimension * (1 << (m*dimension)) * limit**(2*m*dimension) < n**(2*m*m):
            return m
    raise ValueError("Block exceeds the certified bound")


def block_square_roots(n, start, limit, m, square):
    """Find the roots of x² = square (mod n) in one certified block."""
    variable = Symbol("x")
    linear = 2*start % n
    constant = (start*start-square) % n
    dimension = 2*m+1
    powers = [1]
    scales = [limit**i for i in range(dimension)]
    rows = []
    for i in range(m+1):
        multiplier = n**(m-i)
        row = [multiplier*value for value in powers]
        if i < m:
            rows.append(row + [0]*(dimension-len(row)))
            rows.append([0] + row + [0]*(dimension-len(row)-1))
            new_powers = [0]*(len(powers)+2)
            for j, value in enumerate(powers):
                new_powers[j] += constant*value
                new_powers[j+1] += linear*value
                new_powers[j+2] += value
            powers = new_powers
        else:
            rows.append(row)
    lattice = Matrix([[value*scale for value, scale in zip(row, scales)] for row in rows])
    try:
        reduced = lattice.lll()
    except (AssertionError, OverflowError) as error:
        from sympy.external.gmpy import GROUND_TYPES

        if GROUND_TYPES == 'flint':
            raise
        raise RuntimeError('SymPy LLL failed on this basis. Retry with python-flint installed.') from error
    short = list(reduced[0, :])
    assert sum(abs(value) for value in short) < n**m
    roots = Poly.from_list(
        [value//scale for value, scale in reversed(list(zip(short, scales)))],
        gens=variable,
    ).ground_roots()
    return [
        start+int(value) for value in roots
        if value.is_Integer and 0 <= value <= limit
        and ((start+int(value))**2-square) % n == 0
    ]


def block_roots(n, start, limit, m):
    return [root for root in block_square_roots(n, start, limit, m, 1)
            if mixed_root(n, root)]


def parity(n, low, high, work):
    """Return mixed-root parity on the inclusive interval."""
    limit = max(1, isqrt(n)//8)
    m = multiplicity(n, limit) if limit > 1 else None
    result = 0
    for start in range(low, high+1, limit+1):
        local_limit = min(limit, high-start)
        work["blocks"] += 1
        if local_limit <= 1:
            result ^= sum(mixed_root(n, x) for x in range(start, start+local_limit+1)) % 2
        else:
            work["lattices"] += 1
            result ^= len(block_roots(n, start, local_limit, m)) % 2
    return result


def find_factor(n, work):
    return split_from_windows(n, lambda start: parity(n, start+1, start+(n-1)//2, work))


if __name__ == "__main__":
    from argparse import ArgumentParser

    cli = ArgumentParser()
    cli.add_argument("n", type=int)
    cli.add_argument("--start", type=int, default=0)
    cli.add_argument("--split", action="store_true")
    args = cli.parse_args()
    work = Counter()
    if args.split:
        result = find_factor(args.n, work)
    else:
        result = parity(args.n, args.start+1, args.start+(args.n-1)//2, work)
    print(result, dict(work))
