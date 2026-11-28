"""Exact controls for supplied congruence bounds, aggregate bounds and fixed weights."""

from math import gcd, isqrt, prod
from pathlib import Path
import sys

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.weighted_order_hint import (
    _profile_basis, combine_congruences, fixed_weight_divisors, weighted_congruence_divisors,
)


def check_crt_and_basis():
    modulus = 23*31
    hints = ((8, 7), (9, -4), (55, 1))
    moduli = (modulus, 8, 9, 55)
    product, residue = combine_congruences(hints)
    assert product == 8*9*55 and all(residue % e == a % e for e, a in hints)
    for r in range(-100, 101):
        assert gcd(product, r-residue) == prod(gcd(e, r-a) for e, a in hints)
    crt = combine_congruences(((modulus, 0), *hints))
    weights = (2, 3, 1, 2)
    center, radius, d = 25, 7, 5
    rows = _profile_basis(moduli, crt, center, radius, d, weights)
    det = radius**(d*(d-1)//2)*prod(e**(m*(m+1)//2) for e, m in zip(moduli, weights))
    assert prod(row[i] for i, row in enumerate(rows)) == det
    for r in (23, 31):
        gs = (r, *(gcd(e, r-a) for e, a in hints))
        divides = prod(g**m for g, m in zip(gs, weights))
        for row in rows:
            assert all(v % radius**j == 0 for j, v in enumerate(row))
            assert sum(v//radius**j*(r-center)**j for j, v in enumerate(row)) % divides == 0
    bounds = ((8, 7, 4), (9, -4, 3), (55, 1, 11))
    assert weighted_congruence_divisors(modulus, bounds, center=center, radius=radius,
                                        max_dimension=8) == (23,)
    combined = (*combine_congruences(hints), 400)
    expected = tuple(r for r in (23, 31) if gcd(product, r-residue) >= 400)
    assert expected == (23,)
    assert weighted_congruence_divisors(modulus, [combined], center=center, radius=radius,
                                        max_dimension=8) == expected


def check_separation():
    p = 340282366921072654229351831270160596993
    q = 425352958651465639377564999549989158913
    for prime in (p, q):
        # Proth primality certificate: the power of two in prime-1 exceeds sqrt(prime).
        power_two = (prime-1) & -(prime-1)
        assert power_two**2 > prime and pow(3, (prime-1)//2, prime) == prime-1
    modulus = p*q
    assert p < q < 2*p
    e1, e2 = 5 << 61, ((p+1)//2)*7**144
    pairs = ((e1, 1), (e2, -1))
    product, residue = combine_congruences(pairs)
    assert gcd(product, modulus) == 1
    g1, g2 = gcd(e1, p-1), gcd(e2, p+1)
    assert g1 == 2**61 < e1 and 2**127 < g2 < e2
    assert gcd(product, p-residue) == g1*g2
    radius = isqrt(2*modulus)
    lower = isqrt((modulus-1)//2)+1
    # In every dimension D, D(D-1)z²-2Dmz+m(m+1)>=0 for integral 0<=m<D.
    # Consequently log X >= sum z_i² log e_i, with bound_i <= e_i^z_i,
    # precludes even the determinant test without LLL penalties. Include N at z=1/2.
    assert lower**2 <= modulus
    assert radius**4 >= modulus*e1**4                         # first hint, z=1
    assert g2**4 <= e2 and radius**16 >= modulus**4*e2         # second hint, z=1/4
    assert (g1*g2)**25 <= product**8                           # aggregate hint, z=8/25
    assert radius**2500 >= modulus**625 * product**256
    stats = {}
    hints = ((e1, 1, 2**61), (e2, -1, 2**127))
    assert weighted_congruence_divisors(modulus, hints, max_dimension=7, stats=stats) == (p,)
    assert stats['dimension'] == 7 and stats['multiplicities'] == (3, 6, 1)
    assert stats['evaluation_bound_squared'] < (lower**3 * (2**61)**6 * 2**127)**2
    # This reader is given only the multiplicities and public residues. Its
    # a posteriori bound refers to the actual gcds, even if no G_i are supplied.
    blind = {}
    assert fixed_weight_divisors(modulus, pairs, 7, (3, 6, 1), stats=blind) == (p,)
    assert blind == stats
    print('257-bit planted N: recovered a factor with weights (3,6,1) at dimension 7')
    print('Every dimension fails the sufficient determinant test for either singleton '
          'and for their exact combined gcd under a common hint weight')


def check_boundary_cases():
    assert fixed_weight_divisors(23*31, [(8, 7)], 2, (0, 0), center=23, radius=0) == (23,)
    assert weighted_congruence_divisors(23*31, [(8, 0, 8)], center=23, radius=0) == ()
    for run in (lambda: combine_congruences([(6, 1), (9, 2)]),
                lambda: weighted_congruence_divisors(35, [(15, 1, 2)]),
                lambda: weighted_congruence_divisors(35, [(4, 1, 5)]),
                lambda: fixed_weight_divisors(35, [(4, 1)], 3, (1,)),
                lambda: fixed_weight_divisors(35, [(4, 1)], 3, (1, 3))):
        try:
            run()
        except ValueError:
            pass
        else:
            raise AssertionError('expected invalid congruences, bounds or multiplicities to fail')


if __name__ == '__main__':
    check_crt_and_basis()
    check_separation()
    check_boundary_cases()
    print('multi-hint controls passed')
