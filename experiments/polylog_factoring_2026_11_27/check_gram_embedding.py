"""Check exact primitive embeddings, their prescribed roots, and large inputs."""

from math import gcd, isqrt
from random import Random

if __package__:
    from .gram_embedding import embed_gram, orthonormal_coordinates
else:
    from gram_embedding import embed_gram, orthonormal_coordinates


def _bezout(a, b):
    g, h = a, b
    x, y, u, v = 1, 0, 0, 1
    while h:
        q = g // h
        g, h = h, g - q * h
        x, u = u, x - q * u
        y, v = v, y - q * v
    sign = -1 if g < 0 else 1
    return abs(g), sign * x, sign * y


def _check(a, b, c, root):
    d = a*c - b*b
    n = isqrt(d)
    u, v = embed_gram(a, b, c, root)
    assert sum(x*x for x in u) == a and sum(y*y for y in v) == c
    assert sum(x*y for x, y in zip(u, v)) == b
    f = (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0])
    assert gcd(*f) == 1 and sum(x*x for x in f) == d

    # Independently recover the unique lifted scalar from the output basis.
    h, x, y = _bezout(f[0], f[1])
    _, w, z = _bezout(h, f[2])
    bezout = (w*x, w*y, z)
    numerators = [b*x - a*y for x, y in zip(u, v)]
    lift = sum(x*y for x, y in zip(numerators, bezout)) % d
    assert (lift - root) % n == 0 and (lift*lift + a) % d == 0
    assert all((value - lift*x) % d == 0 for value, x in zip(numerators, f))
    return lift % n


def check_small():
    checked = 0
    for n in range(1, 36, 2):
        for a in range(1, 36):
            if gcd(a, n) != 1:
                continue
            roots = [i for i in range(n) if (i*i + a) % n == 0]
            for b in range(-(a // 2), a // 2 + 1):
                if (n*n + b*b) % a:
                    continue
                c = (n*n + b*b) // a
                for root in roots:
                    _check(a, b, c, root)
                    checked += 1
    return checked


def check_large():
    rng = Random(710942)
    checked = 0
    for bits in (100, 160, 256, 384, 512):
        for _ in range(3):
            # a=m²+1 gives a root m of −1 mod a; n=s²+a supplies the
            # modular embedding witness. No prime factors enter the check.
            m = rng.getrandbits(bits // 2 - 3) * 2
            a = m*m + 1
            s = isqrt(1 << (bits - 1)) + m
            s += s % 2
            while gcd(s, a) > 1:
                s += 2
            n = s*s + a
            assert n.bit_length() == bits
            b = n*m % a
            c = (n*n + b*b) // a
            shear = rng.getrandbits(2 * bits) * rng.choice((-1, 1))
            for root in (s, -s):
                _check(a, b + shear*a, c + 2*shear*b + shear*shear*a, root)
                checked += 1
    return checked


def check_invalid():
    cases = ((0,0,1,0), (1,2,1,0), (1,0,3,0), (1,0,4,1),
             (3,0,3,0), (1,0,9,1), (2,1,5,0))
    for case in cases:
        try:
            embed_gram(*case)
        except ValueError:
            continue
        raise AssertionError(f"accepted invalid input {case}")
    return len(cases)


def check_rank_four():
    rng = Random(493864)
    for _ in range(12):
        basis = [[int(i == j) for j in range(4)] for i in range(4)]
        for _ in range(12):
            i, j = rng.sample(range(4), 2)
            k = rng.randrange(-(1 << 30), 1 << 30)
            basis[i] = [x + k*y for x, y in zip(basis[i], basis[j])]
        gram = tuple(tuple(sum(x*y for x, y in zip(left, right)) for right in basis)
                     for left in basis)
        u, v = orthonormal_coordinates(gram)
        assert len(u) == len(v) == 4
        assert sum(x*x for x in u) == gram[0][0]
        assert sum(x*y for x, y in zip(u, v)) == gram[0][1]
        assert sum(y*y for y in v) == gram[1][1]
    for gram in ([[1]*3 for _ in range(3)], [[1,0,0],[0,1,0],[0,0,2]],
                 [[1,0,0],[0,0,1],[0,1,0]], [[1,2,0],[0,1,0],[0,0,1]]):
        try:
            orthonormal_coordinates(gram)
        except ValueError:
            continue
        raise AssertionError(f"accepted invalid Gram {gram}")
    return 12


if __name__ == "__main__":
    print("Small prescribed-root embeddings:", check_small())
    print("Large sheared embeddings (100–512 radius bits):", check_large())
    print("Rank-four unimodular controls:", check_rank_four())
    roots = [_check(26, -3, 9, s) for s in (2, 7, 8, 13)]
    assert gcd(roots[1] - roots[0], 15) == 5 and gcd(roots[2] - roots[0], 15) == 3
    print("Oriented roots recovered at radius 15:", roots)
    print("Rejected invalid inputs:", check_invalid())
