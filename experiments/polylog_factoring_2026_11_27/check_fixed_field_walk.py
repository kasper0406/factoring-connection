"""Check every primitive quaternion fiber on small ternary shells."""

from collections import Counter
from fractions import Fraction
from math import gcd, prod
from random import Random

from sympy import jacobi_symbol

if __package__:
    from .check_fixed_field_packing import sphere_progression
    from .check_four_square_incidence import prime_powers, quadruples
    from .fixed_field_packing import fixed_field_parameters, packing_params
    from .fixed_field_walk import (
        base_points, conjugate_vector, packed_primitive_point, packed_sphere_point,
        primitive_base_points, primitive_scaled_sphere, scaled_sphere,
    )
else:
    from check_fixed_field_packing import sphere_progression
    from check_four_square_incidence import prime_powers, quadruples
    from fixed_field_packing import fixed_field_parameters, packing_params
    from fixed_field_walk import (
        base_points, conjugate_vector, packed_primitive_point, packed_sphere_point,
        primitive_base_points, primitive_scaled_sphere, scaled_sphere,
    )


def shell(norm):
    return {v for v in sphere_progression(norm, (0, 0, 0), 1) if gcd(*v) == 1}


def check_fibers(x, conductor):
    a = x*x+1
    base = primitive_base_points(x)
    assert len(base) == len(set(base)) and set(base) == shell(a)
    assert set(base_points(x, primitive=False)) == set(sphere_progression(a, (0, 0, 0), 1))
    quaternions = [q for q in quadruples(conductor) if gcd(*q) == 1]
    norm = a*conductor*conductor
    expected = shell(norm)
    counts = Counter()
    success = prod(Fraction(p-int(jacobi_symbol(-a, p)), p+1)
                   for p, _ in prime_powers(conductor))
    for s in base:
        passed = 0
        for q in quaternions:
            w = conjugate_vector(q, s)
            assert sum(c*c for c in w) == norm
            if gcd(*w) == 1:
                counts[w] += 1
                passed += 1
        assert Fraction(passed, len(quaternions)) == success
    assert set(counts) == expected and set(counts.values()) == {8}
    rng = Random(16381*x + conductor)
    for _ in range(8):
        assert primitive_scaled_sphere(x, conductor, rng) in expected
    return len(base)*len(quaternions), len(counts)


def check_public_packing():
    rng = Random(81869)
    for n, x, gaussian in ((15, 0, None), (15, 0, (128, 3)), (77, 1, None)):
        params = fixed_field_parameters(n, x) if gaussian is None else packing_params(n, x, *gaussian)
        for _ in range(16):
            w = packed_primitive_point(n, x, rng, gaussian=gaussian)
            assert gcd(*w) == 1 and sum(i*i for i in w) == params.M
    n = (2**256-2**32-977) * (2**256-2**224+2**192+2**96-1)
    x = 2
    params = fixed_field_parameters(n, x)
    assert n.bit_length() == 512 and params.M.bit_length() == 2064
    w = packed_primitive_point(n, x, rng)
    assert gcd(*w) == 1 and sum(i*i for i in w) == params.M
    w = packed_sphere_point(n, x, rng)
    assert sum(i*i for i in w) == params.M
    return params.T.bit_length(), params.M.bit_length()


def check_all_fibers(x, conductor):
    a = x*x+1
    counts = Counter()
    expected = set(sphere_progression(a*conductor**2, (0, 0, 0), 1))
    for g in range(1, conductor+1):
        if conductor % g:
            continue
        m = conductor//g
        quaternions = [q for q in quadruples(m) if gcd(*q) == 1]
        psi = prod(p**(e-1)*(p+1) for p, e in prime_powers(m))
        assert len(quaternions) == 8*psi
        for s in base_points(x, primitive=False):
            for q in quaternions:
                w = conjugate_vector(q, s)
                if gcd(*w, m) == 1:
                    counts[tuple(g*i for i in w)] += Fraction(psi, len(quaternions))
    assert set(counts) == expected and set(counts.values()) == {Fraction(1)}
    rng = Random(39341*x + conductor)
    for _ in range(20):
        assert scaled_sphere(x, conductor, rng) in expected
    return len(expected), sum(gcd(*v) > 1 for v in expected)


if __name__ == '__main__':
    # Squarefree and repeated conductors, primes dividing A, and A with
    # imprimitive base points (x=7 gives the excluded triples (±5,±5,0)).
    fixtures = ((0, 1), (0, 3), (1, 5), (2, 5), (2, 15), (2, 25), (3, 45), (7, 9), (7, 15))
    counts = [check_fibers(*f) for f in fixtures]
    attempts, points = map(sum, zip(*counts))
    print(f"{len(fixtures)} exact shells, {attempts} proposals, {points} primitive points "
          f"each of multiplicity eight, {8*len(fixtures)} sampled points")
    all_fixtures = ((0, 1), (0, 9), (1, 15), (7, 9), (7, 21))
    sizes = [check_all_fibers(x, t) for x, t in all_fixtures]
    all_points, imprimitive = map(sum, zip(*sizes))
    print(f"{len(sizes)} unrestricted shells, {all_points} equal point masses, "
          f"including {imprimitive} imprimitive points")
    conductor_bits, norm_bits = check_public_packing()
    print(f"48 small public samples and one 512-bit input with {conductor_bits}-bit conductor "
          f"and {norm_bits}-bit norm")
