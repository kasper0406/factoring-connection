"""Verify the prescribed-shape reduction and the native four-dimensional construction."""

from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import gcd, isqrt
from random import Random
from sympy.polys.domains import ZZ

if __package__:
    from .check_ternary_sphere import primitive_points
    from .gram_embedding import embed_gram
    from . import ternary_form_reduction as reduction
    from .ternary_sphere import perpendicular_form
else:
    from check_ternary_sphere import primitive_points
    from gram_embedding import embed_gram
    import ternary_form_reduction as reduction
    from ternary_sphere import perpendicular_form


def check_primitive(gram, pair):
    u, v = pair
    dot = reduction._dot
    assert (dot(u,u), dot(u,v), dot(v,v)) == gram
    assert gcd(*(u[i]*v[j]-u[j]*v[i] for i,j in combinations(range(len(u)), 2))) == 1


def check_public_queries(rng):
    count = 0
    for n in range(3, 36, 2):
        for r in range(1,n):
            if gcd(r,n) != 1:
                continue
            gram, h = reduction.gram_for_square_root(n, r*r, rng)
            supplied = r*pow(h,-1,n) % n
            triple = embed_gram(*gram, supplied)
            check_primitive(gram, triple)
            root = reduction.root_from_embedding(gram, *triple)
            assert (h*h*root*root-r*r) % n == 0
            assert (root+supplied) % n == 0
            count += 1
    return count


def check_complements(rng):
    cases = embeddings = axes = 0
    for n in (3,5,9,15,21,35):
        for target in (1, n-1):
            gram, _ = reduction.gram_for_square_root(n, target, rng)
            a = gram[0]
            observed = 0
            for x,y in product(range(n), repeat=2):
                error = -(x*x+y*y+a)
                if error % n:
                    continue
                # Lift this individual residue to the required norm modulo n².
                s,t,h = map(int, ZZ.gcdex(x,y))
                step = (error//n)*pow(2*h,-1,n) % n
                z = x+n*step*s, y+n*step*t
                pair = reduction._quaternary_completion(gram,z)
                check_primitive(gram,pair)
                dropped = any(i == j == 0 for i,j in zip(*pair))
                assert dropped == (x == 0 or y == 0)
                observed += dropped
                embeddings += 1
            assert observed == 2*sum((s*s+a) % n == 0 for s in range(n))
            axes += observed
            cases += 1
    return cases, embeddings, axes


def shell(norm):
    pairs = defaultdict(list)
    radius = isqrt(norm)
    for x in range(-radius,radius+1):
        for y in range(-isqrt(norm-x*x),isqrt(norm-x*x)+1):
            pairs[x*x+y*y].append((x,y))
    return [u+v for s,left in tuple(pairs.items())
            for u in left for v in pairs.get(norm-s,())]


def check_full_distribution():
    total = 0
    for gram in ((1,0,9),(2,1,5),(1,0,25),(26,-3,9)):
        a,b,c = gram
        n = isqrt(a*c-b*b)
        all_pairs = set()
        for u,v in product(shell(a),shell(c)):
            if reduction._dot(u,v) == b and gcd(*(u[i]*v[j]-u[j]*v[i]
                    for i,j in combinations(range(4),2))) == 1:
                all_pairs.add((u,v))
        measured = Counter()
        for x,y in product(range(n),repeat=2):
            norm = x*x+y*y
            if (norm+a)%n:
                continue
            k = ((-a-norm)//n)*pow(2*norm,-1,n) % n
            z = (1+n*k)*x, (1+n*k)*y
            pair = reduction._quaternary_completion(gram,z)
            for perm in permutations(range(4)):
                inversions = sum(perm[i]>perm[j] for i,j in combinations(range(4),2))
                for signs in product((-1,1),repeat=4):
                    if (inversions + signs.count(-1)) % 2:
                        continue
                    transformed = tuple(tuple(s*w[j] for s,j in zip(signs,perm)) for w in pair)
                    measured[transformed] += 1
        assert set(measured) == all_pairs and set(measured.values()) == {4}
        total += len(all_pairs)
    return total


def check_all_frames(rng):
    totals = []
    for n,x in ((9,1),(15,1),(21,0),(35,1),(77,1)):
        a = x*x+1
        histogram = Counter()
        for u,v in product(shell(a),shell(n*n)):
            if reduction._dot(u,v) == -n*x:
                histogram[gcd(*v)] += 1
        for content in range(1,n+1):
            if n % content == 0:
                m = n//content
                norm_fiber = sum((i*i+j*j-1)%m == 0
                                 for i,j in product(range(m),repeat=2))
                assert histogram[content] == 48*norm_fiber
        for _ in range(16):
            u,v = reduction.all_four_line_embeddings(n,x,rng)
            assert (reduction._dot(u,u),reduction._dot(u,v),reduction._dot(v,v)) == (a,-n*x,n*n)
            assert gcd(*v) in histogram
        totals.append((n, sum(histogram.values()), sum(histogram[g] for g in histogram if 1<g<n)))
    return totals


def check_shape_oracle(n, prime_count, rng):
    points = {}
    for point in primitive_points(n):
        key = reduction.reduce_form(perpendicular_form(point))
        points.setdefault(key, point)
    roots = defaultdict(list)
    for s in range(n):
        if gcd(s,n) == 1:
            roots[s*s % n].append(s)
    class_counts = Counter()
    by_target = defaultdict(Counter)
    splits = challenges = 0
    for z in product(range(n), repeat=2):
        target = -reduction._norm(z) % n
        if target not in roots:
            continue
        (a,b,c), _ = reduction.line_gram(n,z)
        key = reduction.reduce_form((a,2*b,c))
        recovered = reduction.residue_matching_norm_minus_one(z, points[key])
        assert recovered*recovered % n == target
        class_counts[key] += 1
        by_target[target][key] += 1
        for s in roots[target]:
            splits += 1 < gcd(s-recovered,n) < n
            challenges += 1
    assert set(class_counts) == set(points) and len(set(class_counts.values())) == 1
    assert all(set(row) == set(points) and len(set(row.values())) == 1 for row in by_target.values())
    assert splits*2**prime_count == challenges*(2**prime_count-2)
    for target in range(1,n):
        if gcd(n,target) == 1:
            root = reduction.root_from_shape_oracle(n,target,lambda modulus,key: points.get(key),rng)
            assert (root is not None) == (target in roots)
            if root is not None:
                assert root*root % n == target
    # Running the real randomized wrapper queries only the canonical form.
    trials = [reduction.factor_from_shape_oracle(n, lambda modulus, key: points[key], rng)
              for _ in range(40)]
    assert any(d is not None for d in trials)
    assert all(d is None or (1 < d < n and n % d == 0) for d in trials)
    return len(points), challenges, splits


if __name__ == "__main__":
    rng = Random(27601)
    print("Public square-root queries and prescribed-root embeddings:", check_public_queries(rng))
    print("Complement cases, primitive 4D embeddings, deletable axes:", *check_complements(rng))
    print("Distinct 4D embeddings with exact uniform mass:", check_full_distribution())
    print("All 4D embeddings and proper contents:", check_all_frames(rng))
    for n, omega in ((15,2),(35,2),(45,2),(105,3)):
        print("N, attained classes, exhaustive challenges, splits:", n,
              *check_shape_oracle(n,omega,rng))
    for n in (3, (2**61-1)*(2**31-1), (2**127-1)*(2**89-1)):
        gram = 1, 0, n*n
        pair = reduction.embed_four(gram,rng)
        check_primitive(gram,pair)
        print("Native 4D embedding at radius bits:", n.bit_length())
