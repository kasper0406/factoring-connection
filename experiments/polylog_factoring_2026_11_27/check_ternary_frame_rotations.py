"""Verify the rational rotations and fixed-axis test for mixed ternary Grams."""

from collections import Counter, defaultdict
from math import gcd, isqrt
from random import Random

from sympy import jacobi_symbol

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .four_square_sampler import _multiply, primitive_four_squares
else:
    from check_mixed_ternary_embeddings import _all_frames
    from four_square_sampler import _multiply, primitive_four_squares


def _cross(u, v):
    a,b,c = u
    x,y,z = v
    return b*z-c*y, c*x-a*z, a*y-b*x


def _conjugate_vector(h, u):
    w, x, y, z = h
    return _multiply(_multiply(h, (0,*u)), (w,-x,-y,-z))[1:]


def quaternion_split_trial(n, x, rng, *, split=None):
    """One public rotation-numerator gcd, drawing h uniformly at norm n."""
    h = primitive_four_squares(n, rng, split=split)
    d = gcd(n, *_conjugate_vector(h, (x,1,0)))
    return d if 1 < d < n else None


def _norm_shell(n):
    pairs = defaultdict(list)
    for a in range(-isqrt(n), isqrt(n)+1):
        bmax = isqrt(n-a*a)
        for b in range(-bmax, bmax+1):
            pairs[a*a+b*b].append((a,b))
    return [pair + other for k in pairs for pair in pairs[k] for other in pairs.get(n-k, ())]


def _fixed_axis_frames(n, x, r):
    """Independent norm-equation generation of all frames fixing (x,1,0)."""
    frames = set()
    for norm in (r,2*r):
        for b in range(1, isqrt(norm//(x*x+1))+1):
            square = norm - (x*x+1)*b*b
            a = isqrt(square)
            if a*a != square or gcd(a,b) != 1:
                continue
            for real in {a, -a}:
                h = (real,b*x,b,0)
                column = _conjugate_vector(h, (-n,0,0))
                assert all(c % norm == 0 for c in column)
                frames.add(tuple(c//norm for c in column))
    return frames


def check_family():
    targets = frames = quaternion_checks = 0
    fixed_orbits = Counter()
    for p,q in ((3,5),(3,7),(5,7),(5,13),(7,11),(7,13)):
        n = p*q
        shell = _norm_shell(n)
        assert len(shell) == 8*(p+1)*(q+1)
        for x in range(n):
            a = x*x+1
            if jacobi_symbol(-a,n) != -1:
                continue
            r = p if jacobi_symbol(-a,p) == 1 else q
            denominators = Counter()
            same_orbit = 0
            fixed = set()
            for u,v in _all_frames(n,x):
                normal = _cross(u,v)
                columns = [tuple(-c for c in v), tuple(n*i+x*j for i,j in zip(u,v)), normal]
                assert _cross(columns[0],columns[1]) == tuple(n*c for c in columns[2])
                assert all(sum(c*c for c in col) == n*n for col in columns)
                denominator = n//gcd(n, *(c for col in columns for c in col))
                assert denominator == n//gcd(*v)
                denominators[denominator] += 1
                if denominator != 1:
                    if sorted(map(abs,u)) == sorted((x,1,0)):
                        same_orbit += 1
                    if u == (x,1,0):
                        fixed.add(v)
            assert denominators == {1:24,r:48}
            from_equation = _fixed_axis_frames(n,x,r)
            assert from_equation == fixed
            assert same_orbit == (48 if fixed else 0)
            assert not fixed or a <= 2*r
            fixed_orbits[bool(fixed)] += 1
            hits = Counter(gcd(n, *_conjugate_vector(h,(x,1,0))) for h in shell)
            count = len(shell)*2//(r+1)
            assert hits == {1:len(shell)-count,r:count}
            targets += 1
            frames += sum(denominators.values())
            quaternion_checks += len(shell)
    return targets, frames, quaternion_checks, dict(fixed_orbits)


if __name__ == "__main__":
    print("Targets, frames, quaternion checks, all-48 axis outcomes:", check_family())
    # Public trial smoke test. The preceding exhaustive count establishes its law.
    rng = Random(540725)
    runs = Counter(quaternion_split_trial(77,1,rng) for _ in range(150))
    assert set(runs) == {None,11}
    print("150 native trials at N=77, x=1:", dict(runs))
