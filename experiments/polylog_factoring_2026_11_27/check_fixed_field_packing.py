"""Exhaustive checks for the fixed-field sphere packing and its local markings."""

from collections import Counter
from math import gcd, isqrt

from sympy import jacobi_symbol, primerange

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .fixed_field_packing import fixed_field_parameters, packing_params, recover_factor
    from .short_jacobi_grams import select_x
else:
    from check_mixed_ternary_embeddings import _all_frames
    from fixed_field_packing import fixed_field_parameters, packing_params, recover_factor
    from short_jacobi_grams import select_x


def sphere_progression(norm, residues, step):
    """Enumerate exactly the signed sphere points in one coordinatewise class."""
    def interval(a, radius):
        return range(-radius + (a + radius) % step, radius + 1, step)

    r, s, t = residues
    for i in interval(r, isqrt(norm)):
        rem = norm - i*i
        for j in interval(s, isqrt(rem)):
            k2 = rem - j*j
            k = isqrt(k2)
            if k*k == k2:
                if (k-t) % step == 0:
                    yield i, j, k
                if k and (-k-t) % step == 0:
                    yield i, j, -k


def local_vectors(n, x, u, ell):
    """Solve the local plane and sphere by eliminating one coordinate."""
    pivot = next(i for i, a in enumerate(u) if a % ell)
    i, j = (s for s in range(3) if s != pivot)
    inverse = pow(u[pivot], -1, ell)
    for a in range(ell):
        for b in range(ell):
            c = (-n*x - u[i]*a - u[j]*b)*inverse % ell
            if (a*a + b*b + c*c - n*n) % ell == 0:
                z = [0]*3
                z[i], z[j], z[pivot] = a, b, c
                yield tuple(z)


def pack(u, v, p, q):
    return tuple(p*a-q*b for a, b in zip(u, v))


def recover_marked(markings):
    """CRT-reconstruct a bounded signed v from its coordinate labels."""
    current, modulus = [0]*3, 1
    for ell, z in markings:
        inv = pow(modulus, -1, ell)
        current = [a + modulus*((b-a)*inv % ell) for a, b in zip(current, z)]
        modulus *= ell
    return tuple(a-modulus if 2*a > modulus else a for a in current)


def check_fixture(n, x, factor):
    assert jacobi_symbol(-(x*x+1), n) == -1
    a = x*x+1
    b, p, q, t, m = fixed_field_parameters(n, x)
    assert b//2 <= (2*x+2)*n < b and b % (2*a) == 0
    primes = [ell for ell in primerange(11, 60) if gcd(ell, a*t*q*n*p) == 1]
    first, second = primes[:2]
    assert first*second > 2*n
    frames = set(_all_frames(n, x))
    contents, coarse_extras = Counter(), 0
    seen = set()
    for u in sphere_progression(a, (0, 0, 0), 1):
        target = {pack(u, v, p, q): v for f, v in frames if f == u}
        coarse = set(sphere_progression(m, tuple(p*i for i in u), q))
        assert set(target) <= coarse
        if x == 0:
            extra = tuple(-t*i for i in u)
            assert coarse - set(target) == {extra}
            assert pack(u, tuple(b*i for i in u), p, q) == extra
            assert gcd(*extra) == t
        coarse_extras += len(coarse) - len(target)
        labels = {}
        for ell in (first, second):
            occupied = {}
            for z in local_vectors(n, x, u, ell):
                for w in sphere_progression(m, pack(u, z, p, q), q*ell):
                    assert w not in occupied
                    occupied[w] = z
            assert set(occupied) == set(target)
            labels[ell] = occupied
        for w, v in target.items():
            assert gcd(*w) == 1
            decoded = recover_marked((ell, labels[ell][w]) for ell in labels)
            assert decoded == v and pack(u, decoded, p, q) == w
            contents[gcd(n, *decoded)] += 1
            seen.add((u, w))
    assert len(seen) == len(frames) == 72
    assert contents == {n: 24, factor: 48}
    assert coarse_extras == (6 if x == 0 else 0)
    decisions = check_recovery(n, x, factor)
    print(f"N={n}, x={x}, primes={first},{second}, coarse extras={coarse_extras}, "
          f"primitive frames={len(seen)}, recovered factor {factor} on {contents[factor]}, "
          f"native oracle calls {sum(decisions.values())} ({decisions[True]} positive)")


def check_recovery(n, x, factor, gaussian=None, primitive=False):
    counts = Counter()

    def decide(norm, modulus, residue):
        points = sphere_progression(norm, residue, modulus)
        answer = any(not primitive or gcd(*w) == 1 for w in points)
        counts[answer] += 1
        return answer

    assert recover_factor(n, x, decide, gaussian=gaussian) == factor
    assert counts[True] and counts[False]
    return counts


def check_alternative_packing():
    n, x, gaussian = 15, 0, (128, 3)
    _, p, q, t, m = packing_params(n, x, *gaussian)
    frames = set(_all_frames(n, x))
    assert (gaussian[0] % gaussian[1]) != 0
    for u in sphere_progression(1, (0, 0, 0), 1):
        coarse = set(sphere_progression(m, tuple(p*i for i in u), q))
        assert coarse == {pack(u, v, p, q) for f, v in frames if f == u}
        assert tuple(-t*i for i in u) not in coarse
        assert all(gcd(*w) == 1 for w in coarse)
    counts = check_recovery(n, x, 3, gaussian=gaussian, primitive=True)
    print(f"N=15, x=0, (B,t)={gaussian}, no coarse extras, primitive recovery "
          f"uses {sum(counts.values())} calls ({counts[True]} positive)")


def check_large():
    # The two standardized secp256k1 and secp256r1 prime field moduli.
    n = (2**256-2**32-977) * (2**256-2**224+2**192+2**96-1)
    x, divisor = select_x(n, 128)
    assert n.bit_length() == 512 and divisor is None and x is not None
    assert jacobi_symbol(-(x*x+1), n) == -1
    a = x*x+1
    b, p, q, t, m = fixed_field_parameters(n, x)
    assert (b//(2*a)).bit_count() == 1 and b//2 <= (2*x+2)*n < b
    assert (p, q, t, m) == (b*b-n*n-2*x*b*n, 2*a*b, b*b+n*n, a*(b*b+n*n)**2)
    u, v = (x, 1, 0), (-n, 0, 0)
    w = pack(u, v, p, q)
    assert sum(c*c for c in w) == m and gcd(*w) == gcd(m, n) == 1
    print(f"512-bit public construction, x={x}, packed norm {m.bit_length()} bits")


if __name__ == '__main__':
    for fixture in ((15, 0, 3), (35, 0, 7), (77, 1, 7), (65, 4, 5)):
        check_fixture(*fixture)
    check_alternative_packing()
    check_large()
