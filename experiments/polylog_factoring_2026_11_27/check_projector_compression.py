"""Exact checks of the Schur compression of a uniform fixed-rank idempotent."""

from collections import Counter, deque
from fractions import Fraction
from math import gcd, isqrt, prod

if __package__:
    from .idempotent_sampler import gl_orders
else:
    from idempotent_sampler import gl_orders


def compress(matrix, modulus):
    """For a good last diagonal entry return P_minus and R as row-major tuples."""
    d = isqrt(len(matrix))
    e = matrix[-1]
    if gcd(e*(1-e), modulus) != 1:
        return None
    inverse = pow(e*(1-e), -1, modulus)
    direction = tuple(matrix[i*d+d-1]*matrix[(d-1)*d+j]*inverse % modulus
                      for i in range(d-1) for j in range(d-1))
    minus = tuple((matrix[i*d+j] - (1-e)*direction[i*(d-1)+j]) % modulus
                  for i in range(d-1) for j in range(d-1))
    return minus, direction


def multiply(left, right, n):
    d = isqrt(len(left))
    return tuple(sum(left[i*d+h]*right[h*d+j] for h in range(d)) % n
                 for i in range(d) for j in range(d))


def trace(matrix):
    d = isqrt(len(matrix))
    return sum(matrix[i*d+i] for i in range(d))


def orbit(d, rank, n):
    """Generate the conjugacy orbit independently by elementary transvections."""
    diagonal = tuple(int(i == j and i < rank) for i in range(d) for j in range(d))
    pending = deque([diagonal])
    visited = {diagonal}
    while pending:
        mat = pending.popleft()
        yield mat
        for i in range(d):
            for j in range(d):
                if i == j:
                    continue
                # Conjugation by I+e_ij adds row j to i, then subtracts col i from j.
                result = list(mat)
                for h in range(d):
                    result[i*d+h] = (result[i*d+h] + result[j*d+h]) % n
                for h in range(d):
                    result[h*d+j] = (result[h*d+j] - result[h*d+i]) % n
                neighbor = tuple(result)
                if neighbor not in visited:
                    visited.add(neighbor)
                    pending.append(neighbor)


def check_orbit(d, k, powers):
    n = prod(p**e for p,e in powers)
    full = pair_count = unit_count = pivot_count = 1
    probability = Fraction(1)
    for p, exp in powers:
        gl = gl_orders(p, exp, d)
        full *= gl[d] // (gl[k]*gl[d-k])
        pair_count *= gl[d-1] // (gl[k-1]*gl[1]*gl[d-k-1])
        unit_count *= p**(exp-1)*(p-1)
        pivot_count *= p**(exp-1)*(p-2)
        probability *= Fraction((p**k-1)*(p**(d-k)-1)*(p-2), (p**d-1)*p)
    total = 0
    fibers = Counter()
    for mat in orbit(d, k, n):
        total += 1
        compressed = compress(mat, n)
        if compressed is None:
            continue
        minus, direction = compressed
        assert multiply(minus, minus, n) == minus
        assert multiply(direction, direction, n) == direction and trace(direction) % n == 1
        assert not any(multiply(minus, direction, n)) and not any(multiply(direction, minus, n))
        fibers[compressed, mat[-1]] += 1
    good = sum(fibers.values())
    assert total == full and Fraction(good,total) == probability
    assert len(fibers) == pair_count*pivot_count
    assert all(count == unit_count for count in fibers.values())
    return good, total


def fixture(pivot, n):
    # A rank-two conjugate of diag(1,1,0,0) at every prime of n.
    e = pivot
    return tuple(x % n for x in (1,0,0,0, 0,1-e,0,1, 0,0,0,0, 0,e*(1-e),0,e))


def mixed_scalar():
    n = 15
    minus, direction = compress(fixture(2,n), n)
    ranks = {}
    for t in range(n):
        result = tuple((x+t*y) % n for x,y in zip(minus, direction))
        if multiply(result, result, n) == result:
            assert (trace(result)-trace(minus)) % n == t
            # The displayed fixture has result diag(1,t,0).
            ranks[t] = tuple(1 + int(t % p != 0) for p in (3,5))
        assert (t in ranks) == (t*t % n == t)
    assert ranks == {0:(1,1), 1:(2,2), 6:(1,2), 10:(2,1)}
    for e, expected in ((0,(15,1)), (1,(1,15)), (3,(3,1)), (6,(3,5))):
        mat = fixture(e, n)
        assert multiply(mat, mat, n) == mat and compress(mat, n) is None
        assert (gcd(e,n), gcd(1-e,n)) == expected
    return ranks


if __name__ == "__main__":
    for d, k, powers in ((2,1,((2,1),)), (3,1,((3,1),)), (3,2,((3,1),)),
                         (4,2,((3,1),)), (2,1,((3,2),)), (2,1,((3,1),(5,1)))):
        n = prod(p**e for p,e in powers)
        good, total = check_orbit(d, k, powers)
        print(f"d={d}, k={k}, N={n}: {good} good pivots / {total} projectors")
    print("N=15 interpolation, scalar -> ranks mod (3,5):", mixed_scalar())
