"""Reduce scalar square roots and factoring to prescribed perpendicular lattices."""

from math import gcd, isqrt
from itertools import permutations
from pathlib import Path
import sys

from sympy.polys.domains import ZZ

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.genusleak.forms import reduce_form

if __package__:
    from .four_square_sampler import gaussian_norm_pair, primitive_four_squares
    from .four_square_incidence import sample_content
    from .gram_embedding import orthonormal_coordinates
    from .quaternion_gram import embed_four_quaternion
    from .ternary_sphere import hopf_point
    from .ternary_embedding import _mul, _norm, residue_matching_norm_minus_one
else:
    from four_square_sampler import gaussian_norm_pair, primitive_four_squares
    from four_square_incidence import sample_content
    from gram_embedding import orthonormal_coordinates
    from quaternion_gram import embed_four_quaternion
    from ternary_sphere import hopf_point
    from ternary_embedding import _mul, _norm, residue_matching_norm_minus_one


def _dot(u, v):
    return sum(x*y for x, y in zip(u, v))


def _cross(u, v):
    x, y, z = u
    a, b, c = v
    return y*c-z*b, z*a-x*c, x*b-y*a


def line_gram(n, z):
    """Gram (A,B,C) for Z*z + n Z[i], plus h with Norm(z)=h² A."""
    if n <= 1 or n % 2 == 0 or gcd(n, _norm(z)) != 1:
        raise ValueError("expected an odd modulus >1 and a unit Gaussian residue")
    x, y = z
    h = gcd(x, y)
    x, y = x//h, y//h
    s, t, _ = map(int, ZZ.gcdex(x, y))
    # (x,y), (-t,s) have determinant 1. The first vector generates the line mod n.
    a = x*x+y*y
    b = n*(y*s-x*t)
    c = n*n*(s*s+t*t)
    return (a, b, c), h


def gram_for_square_root(n, target, rng):
    """Return a Gram and h; a primitive embedding gives a root of target scaled by h."""
    if n <= 1 or n % 2 == 0 or gcd(target, n) != 1:
        raise ValueError("expected an odd modulus >1 and a unit target")
    z = gaussian_norm_pair(n, -target, rng)
    gram, h = line_gram(n, z)
    return gram, h


def root_from_embedding(gram, u, v):
    """Read a square root of -A modulo det(gram) from a marked primitive embedding."""
    a, b, c = gram
    determinant = a*c-b*b
    if a <= 0 or determinant <= 1 or (_dot(u,u), _dot(u,v), _dot(v,v)) != gram:
        raise ValueError("vectors must have the specified positive definite Gram matrix")
    normal = _cross(u, v)
    s, t, g = map(int, ZZ.gcdex(normal[0], normal[1]))
    j, k, content = map(int, ZZ.gcdex(g, normal[2]))
    if content != 1:
        raise ValueError("the embedding must have primitive normal")
    coefficient = _dot((j*s, j*t, k), _cross(normal, u)) % determinant
    if (coefficient*coefficient + a) % determinant:
        raise ArithmeticError("embedding did not yield the required modular root")
    return coefficient


def _random_norm_fiber(n, fixed, rng):
    while True:
        w = rng.randrange(n), rng.randrange(n)
        d = _norm(w)
        if gcd(d, n) == 1:
            inverse = pow(d, -1, n)
            return tuple(x*inverse % n for x in _mul(fixed, _mul(w, w)))


def root_from_shape_oracle(n, target, oracle, rng):
    """One call on a random reduced form; return a root of target or None.

    For every unit square target, the queried oriented form class is uniform
    on the same occupied genus. The oracle receives only n and that form.
    """
    if n <= 1 or n % 2 == 0 or gcd(n,target) != 1:
        raise ValueError("expected an odd modulus >1 and a unit target")
    z = _random_norm_fiber(n, gaussian_norm_pair(n, -target, rng), rng)
    (a,b,c), _ = line_gram(n, z)
    point = oracle(n, reduce_form((a, 2*b, c)))
    if point is None or len(point) != 3 or _dot(point,point) != n*n or gcd(*point) != 1:
        return None
    return residue_matching_norm_minus_one(z, point)


def factor_from_shape_oracle(n, oracle, rng):
    """One random-square factoring trial, using a prescribed-shape oracle."""
    if n <= 1 or n % 2 == 0:
        raise ValueError("expected an odd modulus >1")
    while True:
        secret = rng.randrange(n)
        if gcd(secret, n) == 1:
            break
    root = root_from_shape_oracle(n, secret*secret % n, oracle, rng)
    if root is None:
        return None
    divisor = gcd(root-secret, n)
    return divisor if 1 < divisor < n else None


def _quaternary_completion(gram, z):
    """Glue the given binary lattice to the Gaussian line of z, with Norm(z)=-A mod D."""
    a, b, c = gram
    d = a*c-b*b
    n = isqrt(d)
    if (_norm(z)+a) % d:
        raise ValueError("complementary Gaussian norm must be −A modulo the determinant")
    (aa, bb, cc), h = line_gram(n, z)
    scale = pow(h*aa, -1, d)
    r1, r2, r3 = (-b*scale % d, a*scale % d, -h*bb*scale % d)
    # Basis e1,e2,e3,(r1*e1+r2*e2+r3*e3+e4)/D in the orthogonal sum.
    x = (a*r1+b*r2)//d
    y = (b*r1+c*r2)//d
    v = (aa*r3+bb)//d
    w = (a*r1*r1+2*b*r1*r2+c*r2*r2+aa*r3*r3+2*bb*r3+cc)//(d*d)
    return orthonormal_coordinates([[a,b,0,x], [b,c,0,y], [0,0,aa,v], [x,y,v,w]])


def embed_four(gram, rng):
    """Sample uniformly among all ordered primitive embeddings of this Gram in Z⁴."""
    a, b, c = gram
    d = a*c-b*b
    if a <= 0 or d <= 1:
        raise ValueError("expected a positive Gram of determinant greater than one")
    n = isqrt(d)
    if n*n != d or n % 2 == 0 or gcd(a,n) != 1:
        raise ValueError("expected determinant n² for odd n coprime to A")
    z = _random_norm_fiber(n, gaussian_norm_pair(n, -a, rng), rng)
    norm = _norm(z)
    # Multiplication by 1+k*n lifts the required norm without changing the line.
    k = ((-a-norm)//n)*pow(2*norm, -1, n) % n
    pair = _quaternary_completion(gram, tuple((1+k*n)*v for v in z))
    perm = tuple(permutations(range(4)))[rng.randrange(24)]
    signs = rng.randrange(16)
    return tuple(tuple((-1 if signs >> i & 1 else 1)*w[j] for i,j in enumerate(perm))
                 for w in pair)


def all_four_line_embeddings(n, x, rng):
    """Uniform over ALL ordered 4D frames of Gram (x²+1,−nx,n²), including imprimitive."""
    a = x*x+1
    if n <= 1 or n % 2 == 0 or gcd(a,n) != 1:
        raise ValueError("expected an odd modulus >1 coprime to x²+1")
    while True:
        content = sample_content(n, 1, rng)
        m = n//content
        if gcd(*hopf_point(primitive_four_squares(m,rng))) == 1:
            break
    u, v = embed_four_quaternion((a,-m*x,m*m),rng)
    return u, tuple(content*y for y in v)
