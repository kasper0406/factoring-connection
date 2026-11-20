"""Exact checks specified in prereg.txt. Uses only the Python standard library."""

from collections import deque
from fractions import Fraction
from math import gcd, prod


N = 2773
q = 31
base = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)


def rank(rows, binary=False):
    pivots = {}
    for row in rows:
        v = [a % 2 if binary else Fraction(a) for a in row]
        for j in sorted(pivots):
            c = v[j]
            if c:
                v = [a - c * b for a, b in zip(v, pivots[j])]
                if binary:
                    v = [a % 2 for a in v]
        j = next((i for i, a in enumerate(v) if a), None)
        if j is not None:
            pivots[j] = v if binary else [a / v[j] for a in v]
    return len(pivots)


def base_ranks(rows):
    return tuple(rank(rows, binary) - rank([r[-1:] for r in rows], binary)
                 for binary in (False, True))


def check_table():
    q_powers = [1]
    r = q % N
    while r != 1:
        q_powers.append(r)
        r = r * q % N
    d = len(q_powers)
    m = len(base)
    labels = {}
    representatives = {}
    queue = deque()

    def add_coset(rep, word):
        representatives[rep] = word
        queue.append(rep)
        for k, power in enumerate(q_powers):
            residue = rep * power % N
            assert residue not in labels
            labels[residue] = (rep, k)

    add_coset(1, [0] * m)
    rows = [[0] * m + [d]]
    transitions = {}
    while queue:
        rep = queue.popleft()
        for j, b in enumerate(base):
            word = representatives[rep].copy()
            word[j] += 1
            z = rep * q**3 * b % N
            if z not in labels:
                add_coset(z, word)
            dest, k = labels[z]
            v = [a - c for a, c in zip(word, representatives[dest])]
            row = v + [3 * sum(v) - k]
            assert prod(pow(p, a, N) for p, a in zip((*base, q), row)) % N == 1
            rows.append(row)
            transitions[rep, b] = dest

    z = transitions[1, 2]
    t = 1
    while z != 1:
        z = transitions[z, 2]
        t += 1
    order_multiple = t * d
    assert pow(2, order_multiple, N) == 1
    odd = order_multiple
    while odd % 2 == 0:
        odd //= 2
    z = pow(2, odd, N)
    while z * z % N != 1:
        z = z * z % N
    factor = gcd(z - 1, N)
    ranks = base_ranks(rows)
    assert ranks[0] == m and ranks[1] >= m - 2
    assert 1 < factor < N
    print(f"table: {len(labels)} residues, {len(representatives)} cosets, q period {d}")
    print(f"table: {len(rows)-1} edges, base ranks Z/F2 {ranks}, smoothness tests 0")
    print(f"base 2: coset cycle {t}, order multiple {order_multiple}, gcd {factor}")


def check_form():
    primes = base[:5]
    points = [q] + [q**3 * b for b in primes]
    t = pow(prod(-a for a in points), -1, N)

    def g(x):
        return 1 + (x - t) * prod(x - a for a in points)

    assert g(0) != 0 and g(0) % N == 0
    rows = []
    divisors = []
    residuals = []
    for i, a in enumerate(points):
        assert g(a) == 1
        cofactor = a * g(a)
        power = 0
        while cofactor % q == 0:
            cofactor //= q
            power += 1
        assert power == (3 if i else 1)
        assert cofactor == (primes[i-1] if i else 1)
        divisors.append(q**power)
        residuals.append(cofactor)
        rows.append([int(b == cofactor) for b in primes] + [power])
    assert base_ranks(rows) == (5, 5) and rank(rows, True) == len(rows)
    print(f"form: degree {len(points)+1}, t {t}, g(0) mod N {g(0)%N}")
    print(f"form: divisors {divisors}, remaining values {residuals}")
    print(f"form: base ranks Z/F2 {base_ranks(rows)}, binary nullity 0")


def check_constant_residual():
    p = 11
    large = 101
    points = [large**(2*i+1) for i in range(5)]
    t = p * pow(prod(-a for a in points), -1, N)

    def g(x):
        return p + (x - t) * prod(x - a for a in points)

    roots = [a % p for a in points]
    assert t % p == 0 and g(0) % N == 0
    assert 0 not in roots and len(set(roots)) == len(roots)
    rows = []
    for i, a in enumerate(points):
        deriv = (a - t) * prod(a - b for b in points if b != a)
        assert deriv % p and g(a) == p
        assert a * g(a) // large**(2*i+1) == p
        rows.append([int(r == a % p) for r in roots] + [2*i+1])
    assert base_ranks(rows) == (4, 4) and rank(rows, True) == 5
    print(f"constant residual: all {p}, roots {roots}, base ranks Z/F2 {base_ranks(rows)}")


if __name__ == "__main__":
    check_table()
    check_form()
    check_constant_residual()
