"""Compare exact Hecke residue counts with independent coordinate enumeration."""

from collections import Counter
from itertools import product
from math import gcd

from sympy import factorint

if __package__:
    from .check_fixed_field_packing import local_vectors, pack, sphere_progression
    from .factored_ternary_hecke import FactoredTernaryHecke, prime_rotations, primitive_residue_genus
    from .fixed_field_packing import packing_params, recover_factor
else:
    from check_fixed_field_packing import local_vectors, pack, sphere_progression
    from factored_ternary_hecke import FactoredTernaryHecke, prime_rotations, primitive_residue_genus
    from fixed_field_packing import packing_params, recover_factor


def check_all_classes():
    checked = 0
    for p in (3, 5, 7, 11, 13, 17, 19):
        hs = list(prime_rotations(p))
        assert len(set(hs)) == len(hs) == p+1
        assert all(sum(a*a for a in h) == p for h in hs)
    for x, p, q, moduli in [(0, 3, 5, (1, 4, 7)), (1, 5, 13, (2, 7, 8)),
                            (2, 3, 7, (1, 5, 10, 16)), (7, 3, 7, (2, 5, 8, 10, 25))]:
        a = x*x+1
        target = a*(p*q)**2
        points = list(sphere_progression(target, (0, 0, 0), 1))
        for m, primitive in product(moduli, (False, True)):
            reader = FactoredTernaryHecke(x, p, q, m, primitive=primitive)
            actual = Counter(tuple(t % m for t in v) for v in points if not primitive or gcd(*v) == 1)
            signs = [1 if pow(-a, (r-1)//2, r) == 1 else -1 for r in (p,q)]
            assert actual.total() == reader.base_count*(p+1-signs[0]-primitive)*(q+1-signs[1]-primitive)
            size, average = primitive_residue_genus(x, p, q, factorint(m), primitive=primitive)
            assert size*average == sum(n for c, n in actual.items() if gcd(m, *c) == 1)
            two_part = m & -m
            local_classes = 0
            for c in product(range(m), repeat=3):
                assert reader.count(c) == actual[c], (x, p, q, m, primitive, c)
                delta = sum(t*t for t in c)-target
                locally_liftable = gcd(m, *c) == 1 and delta % (m//two_part) == 0
                if two_part == 2:
                    locally_liftable &= sum(t % 2 for t in c) == target % 4
                elif two_part >= 4:
                    locally_liftable &= delta % (2*two_part) == 0
                local_classes += locally_liftable
                checked += 1
            assert local_classes == size
    return checked


def check_packed_classes():
    n, x, p, q = 77, 1, 257, 397
    s = packing_params(n, x, 310, 1)
    assert s.T == p*q
    ell, u = 13, (1, 1, 0)
    fine_modulus = s.Q*ell
    residues = {tuple(t % fine_modulus for t in pack(u, z, s.P, s.Q))
                for z in local_vectors(n, x, u, ell)}
    counts = {}
    for primitive in (False, True):
        reader = FactoredTernaryHecke(x, p, q, fine_modulus, primitive=primitive)
        counts = Counter()
        for c in residues:
            counts[c] = reader.count(c)
            expected = sum(not primitive or gcd(*w) == 1
                           for w in sphere_progression(s.M, c, fine_modulus))
            assert counts[c] == expected
        assert 0 in counts.values() and max(counts.values()) >= 1
        reversed_order = FactoredTernaryHecke(x, q, p, fine_modulus, primitive=primitive)
        assert all(reversed_order.count(c) == value for c, value in counts.items())
        coarse = FactoredTernaryHecke(x, p, q, s.Q, primitive=primitive)
        c = tuple(s.P*i % s.Q for i in u)
        direct = sum(not primitive or gcd(*w) == 1 for w in sphere_progression(s.M, c, s.Q))
        assert coarse.count(c) == direct == counts.total()
    print(f'N={n}, T={p}*{q}: {len(residues)} fine residues, {counts.total()} points; '
          f'{reader.forward_images} forward images, {len(reader.table)} occupied table keys, '
          f'{len(reader.inverse_rotations)} backward rotations per valid query')


def check_recovery():
    n, x, p, q = 15, 0, 29, 101
    params = packing_params(n, x, 52, 1)
    assert params.T == p*q
    readers = {}
    queries = 0

    def decide(norm, modulus, residue):
        nonlocal queries
        assert norm == params.M
        if modulus not in readers:
            readers[modulus] = FactoredTernaryHecke(x, p, q, modulus, primitive=True)
        queries += 1
        return readers[modulus].count(residue) > 0

    assert recover_factor(n, x, decide, gaussian=(52, 1)) == 3
    print(f'N={n}, T={p}*{q}: recovered 3 in {queries} individual primitive decisions '
          f'at {len(readers)} moduli')


if __name__ == '__main__':
    print(f'{check_all_classes()} primitive/ordinary residue counts agree with brute force')
    check_packed_classes()
    check_recovery()
