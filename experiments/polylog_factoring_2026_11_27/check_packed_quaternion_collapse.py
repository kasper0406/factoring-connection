"""Independent inverse-packing tests, including a norm/modulus with two hidden inputs."""

from collections import Counter
from itertools import product
from math import gcd, isqrt
from random import Random

if __package__:
    from .check_fixed_field_packing import local_vectors, pack, sphere_progression
    from .check_mixed_ternary_embeddings import _all_frames
    from .factored_ternary_hecke import prime_rotations
    from .fixed_field_packing import packing_params
    from .fixed_field_walk import base_points, conjugate_vector
    from .four_square_sampler import _multiply
    from .packed_quaternion_collapse import (
        collapse, gaussian_decompositions, parameter_candidates, restore, scalar_points, unpack_candidates,
    )
else:
    from check_fixed_field_packing import local_vectors, pack, sphere_progression
    from check_mixed_ternary_embeddings import _all_frames
    from factored_ternary_hecke import prime_rotations
    from fixed_field_packing import packing_params
    from fixed_field_walk import base_points, conjugate_vector
    from four_square_sampler import _multiply
    from packed_quaternion_collapse import (
        collapse, gaussian_decompositions, parameter_candidates, restore, scalar_points, unpack_candidates,
    )


def norm(v):
    return sum(a*a for a in v)


def check_gaussian():
    rng = Random(6711)
    for primes in [(5,), (13,), (1861,), (5, 13), (17, 29), (13, 577), (29, 101)]:
        t = primes[0] if len(primes) == 1 else primes[0]*primes[1]
        ref = {(b, isqrt(t-b*b)) for b in range(1, isqrt(t)+1)
               if isqrt(t-b*b)**2 == t-b*b}
        for _ in range(4):
            assert set(gaussian_decompositions(primes, rng)) == ref


def check_components(x, b, primes, ell):
    n, t = 15, 1
    params = packing_params(n, x, b, t)
    assert params.T == primes[0]*primes[1]
    m = params.Q*ell
    rng = Random(831408)
    frames = {pack(u, v, params.P, params.Q): (u, v) for u, v in _all_frames(n, x)}
    spinors = [_multiply(k, g) for g, k in product(prime_rotations(primes[0]), prime_rotations(primes[1]))]
    seen = set()
    counts = Counter()
    for base in base_points(x):
        for h in spinors:
            w = conjugate_vector(h, base)
            if w not in frames:
                continue
            u, v = frames[w]
            c = tuple(j % m for j in w)
            valid = unpack_candidates(params.M, m, c, primes, rng)
            (mark,) = (mark for mark in valid if mark.packing.n == n and mark.u == u)
            plus, minus = collapse(mark, h, base)
            assert norm(plus) == 4*(x*x+1) and norm(minus) == 4*(x*x+1)*n*n
            assert _multiply((0, *u), plus) == _multiply(plus, (0, *base))
            assert _multiply((0, *u), minus) == tuple(-i for i in _multiply(minus, (0, *base)))
            assert gcd(n, *minus) == gcd(n, *v)
            assert restore(mark, base, plus, minus) == h
            assert restore(mark, base, (0,0,0,0), minus) is None
            counts[gcd(n, *minus)] += 1
            seen.add(w)
    assert seen == set(frames) and counts == {n: 24, (5 if x else 3): 48}
    checks = 0
    for u in base_points(x):
        for z in local_vectors(n, x, u, ell):
            c = tuple(j % m for j in pack(u, z, params.P, params.Q))
            exact = tuple(sorted(sphere_progression(params.M, c, m)))
            assert scalar_points(params.M, m, c, primes, rng=rng) == exact
            checks += 1
    if x == 0:
        extra = tuple(-params.T*j for j in (1,0,0))
        assert all((j-params.P*i) % params.Q == 0 for j,i in zip(extra, (1,0,0)))
        assert gcd(*extra) == params.T
        assert not unpack_candidates(params.M, m, extra, primes, rng)
        assert not scalar_points(params.M, m, extra, primes, rng=rng)
    print(f'x={x}, T={params.T}: {checks} exact scalar marks, V gcd distribution {dict(counts)}')


def check_ambiguity():
    primes, m = (13, 577), 7548
    target = (13*577)**2
    rng = Random(1367)
    options = parameter_candidates(target, m, primes, rng)
    assert {(o.n, o.t, o.ell, o.params.B) for o in options} == {(15,3,17,74), (35,2,37,51)}
    for o in options:
        u, v = (0, 1, 0), (-o.n, 0, 0)
        w = pack(u, v, o.params.P, o.params.Q)
        c = tuple(i % m for i in w)
        kept = unpack_candidates(target, m, c, primes, rng)
        assert {mark.packing.n for mark in kept} == {o.n}
        assert scalar_points(target, m, c, primes, rng=rng) == (w,)
    # A fully factored packing can also use one split conductor prime.
    n, x, b, t, ell = 15, 0, 31, 2, 11
    p = packing_params(n, x, b, t)
    assert p.T == 1861
    w = pack((0,1,0), (-n,0,0), p.P, p.Q)
    m = p.Q*ell
    assert scalar_points(p.M, m, w, (p.T,), rng=rng) == (w,)
    print('The same (M,m) encodes N=15 or 35. Each coordinate mark recovers its own input.')


if __name__ == '__main__':
    check_gaussian()
    check_components(0, 52, (29, 101), 11)
    check_components(1, 62, (13, 313), 11)
    check_ambiguity()
