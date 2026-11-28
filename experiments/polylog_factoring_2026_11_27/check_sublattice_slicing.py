"""Check the interval between a sublattice's coordinate slice and projection."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from random import Random

from sympy import Matrix, eye
from sympy.matrices.normalforms import hermite_normal_form

if __package__:
    from .sublattice_sampler import lattice_from_types
else:
    from sublattice_sampler import lattice_from_types


def sliced_interval(lattice):
    """Return column bases of slice and projection, cyclic generator, and index."""
    columns = lattice[:-1, :]
    lower = columns[:, :-1]
    upper = hermite_normal_form(columns)
    index = int(lower.det()) // int(upper.det())
    return lower, upper, columns[:, -1], index


def points(generators, n, d):
    vectors = {tuple([0]*d)}
    for v in generators:
        vectors = {tuple((x+a*y) % n for x,y in zip(u,v))
                   for u in vectors for a in range(n)}
    return frozenset(vectors)


def check_space(generators, n):
    d = len(generators[0])
    subspace = points(generators, n, d)
    slice_points = frozenset(x[:-1] for x in subspace if x[-1] == 0)
    projected = frozenset(x[:-1] for x in subspace)
    lattice = hermite_normal_form(Matrix(generators).T.row_join(n*eye(d)))
    lower, upper, b, index = sliced_interval(lattice)
    assert index * len(slice_points) == len(projected) and n % index == 0
    # Check the cyclic generator directly in the quotient.
    h = tuple(map(int,b))
    first_return = next(i for i in range(1,n+1)
                        if tuple(i*x % n for x in h) in slice_points)
    assert first_return == index
    return index, (slice_points, projected), (lower, upper, b)


def field_bases(p, d, k):
    """The unique reduced row echelon bases of all k-subspaces over F_p."""
    for pivots in combinations(range(d), k):
        free = [(i,j) for i,col in enumerate(pivots)
                for j in range(col+1,d) if j not in pivots]
        for values in product(range(p), repeat=len(free)):
            rows = [[int(col == j) for j in range(d)] for col in pivots]
            for (i,j),value in zip(free, values):
                rows[i][j] = value
            yield rows


def field_check(p, k):
    total = horizontal = vertical = 0
    fibers = Counter()
    for basis in field_bases(p, 2*k, k):
        total += 1
        index, flag, _ = check_space(basis, p)
        under = len(flag[0]) == p**k
        contains = len(flag[1]) < p**k
        assert not (under and contains)
        horizontal += under
        vertical += contains
        assert (index == p) == (not under and not contains)
        if index == p:
            fibers[flag] += 1
    assert Fraction(horizontal, total) == Fraction(vertical, total) == Fraction(1, p**k+1)
    assert all(count == p-1 for count in fibers.values())
    return sum(fibers.values()), total


def ring_lines(n, primes):
    lines = {points((v,), n, 2): v for v in product(range(n), repeat=2) if gcd(*v,n) == 1}
    indices = Counter(check_space((v,), n)[0] for v in lines.values())
    probability = Fraction(1)
    for p in primes:
        probability *= Fraction(p-1,p+1)
    assert Fraction(indices[n], len(lines)) == probability
    return dict(sorted(indices.items()))


def mixed_interval():
    n = 15
    index, _, (lower, upper, b) = check_space(((1,0,0,1), (0,1,0,0)), n)
    assert index == n
    intermediate = {}
    for t in range(n):
        middle = hermite_normal_form(lower.row_join(t*b))
        label = int(middle.det()) // int(upper.det())
        assert label == gcd(t,n)
        intermediate.setdefault(label, set()).add(tuple(middle))
    assert set(intermediate) == {1,3,5,15}
    assert all(len(bases) == 1 for bases in intermediate.values())
    return sorted(intermediate)


if __name__ == "__main__":
    for p,k in ((2,1),(3,1),(5,1),(2,2),(3,2),(2,3)):
        good, total = field_check(p,k)
        print(f"F_{p}, dimension {2*k}: {good} good / {total} subspaces")
    for n, primes in ((9,(3,)), (15,(3,5)), (25,(5,))):
        print(f"N={n}, line index counts:", ring_lines(n,primes))
    print("N=15, rank two good interval indices:", mixed_interval())
    # Exercise the HNF returned by the public sampler at larger moduli.
    rng = Random(7109)
    for n,k in ((105,2),(3**4*5**2,2),(1000003*1000033,3)):
        for _ in range(6):
            _, _, _, index = sliced_interval(lattice_from_types(n,2*k,k,{},rng))
            assert n % index == 0
