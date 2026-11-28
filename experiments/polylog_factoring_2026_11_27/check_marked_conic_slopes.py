"""Compare exact marked slope searches with direct ternary frames and scalar queries."""

from collections import defaultdict
from itertools import product
from math import gcd
from random import Random

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .factored_gaussian_packing import certified_prime
    from .fixed_field_packing import packing_params
    from .marked_conic_slopes import marked_slope_partners, prime_conductor_points
    from .packed_quaternion_collapse import parameter_candidates
    from .packed_ternary_decision import _residue_candidates, brute_sphere_count
    from .short_jacobi_grams import first_vector_orbits
else:
    from check_mixed_ternary_embeddings import _all_frames
    from factored_gaussian_packing import certified_prime
    from fixed_field_packing import packing_params
    from marked_conic_slopes import marked_slope_partners, prime_conductor_points
    from packed_quaternion_collapse import parameter_candidates
    from packed_ternary_decision import _residue_candidates, brute_sphere_count
    from short_jacobi_grams import first_vector_orbits


def check_small_marks():
    cases = positives = 0
    for n, x in ((15, 0), (35, 4), (77, 1), (77, 8), (143, 12)):
        direct = defaultdict(set)
        for u, v in _all_frames(n, x):
            direct[u].add(v)
        for u in sorted(set(first_vector_orbits(x)) | set(direct)):
            for ell in (3, 7, 11):
                if gcd(ell, 2*n*(x*x+1)) != 1:
                    continue
                for z in _residue_candidates(u, -n*x, n*n, ell):
                    expected = sorted(v for v in direct[u] if tuple(i % ell for i in v) == z)
                    for method in ("enumerate", "pooled"):
                        assert list(marked_slope_partners(n, x, u, z, {ell: 1}, method=method)) == expected
                    cases += 1
                    positives += bool(expected)
    return cases, positives


def check_composite_digits():
    n, x, u = 15, 16, (6, 10, 11)
    rows = {v for first, v in _all_frames(n, x) if first == u}
    factors = {151: 2, 11: 2}
    modulus = 151**2 * 11**2
    lost = set()
    for v in rows:
        work = {}
        assert marked_slope_partners(n, x, u, v, factors, method="enumerate", stats=work) == (v,)
        assert marked_slope_partners(n, x, u, v, factors) == (v,)
        lost.add(modulus // work["slope_modulus"])
    assert rows and lost == {121}
    n, x, u = 15, 15, (0, 1, 15)
    work = {}
    actual = marked_slope_partners(n, x, u, (0, 0, 0), {}, stats=work)
    assert work["transport_denominator"] % n == 0 and work["proper_factor"] in (3, 5)
    assert set(actual) == {v for first, v in _all_frames(n, x) if first == u}

    # Exhaust all triples, including locally impossible ones, over a prime square.
    n, x, u, modulus = 35, 0, (0, 1, 0), 9
    direct = {v for first, v in _all_frames(n, x) if first == u}
    for z in product(range(modulus), repeat=3):
        expected = {v for v in direct if tuple(i % modulus for i in v) == z}
        for method in ("enumerate", "pooled"):
            assert set(marked_slope_partners(n, x, u, z, {3: 2}, method=method)) == expected
    return len(rows), lost, modulus**3


def check_scalar():
    n, x, u, ell = 77, 1, (1, 1, 0), 13
    params = packing_params(n, x, 163, 34)
    modulus = params.Q*ell
    assert params.T == 6880493
    packings = parameter_candidates(params.M, modulus, (params.T,), Random(9302))
    assert len(packings) == 1 and packings[0].n == n
    hits = queries = 0
    for z in _residue_candidates(u, -n*x, n*n, ell):
        residue = tuple((params.P*i - params.Q*j) % modulus for i, j in zip(u, z))
        count = brute_sphere_count(params.M, modulus, residue)
        for method in ("enumerate", "pooled"):
            points = prime_conductor_points(params.M, modulus, residue, params.T,
                                            method=method, rng=Random(129))
            assert len(points) == count
            assert all(sum(a*a for a in w) == params.M and
                       tuple(a % modulus for a in w) == residue for w in points)
        queries += 1
        hits += count
    return queries, hits


def check_full_menu():
    n, x, u = 251*257, 0, (0, 1, 0)
    observed = {}
    complete = set(marked_slope_partners(n, x, u, (0, 0, 0), {}, method="enumerate"))
    for ell in (3, 11, 31, 127, 263, 521):
        marks = {(-n % ell, 0, 0)}
        for k in range(ell):
            den = (k*k+1) % ell
            if den:
                t = 2*n*pow(den, -1, ell) % ell
                marks.add(((t-n) % ell, 0, k*t % ell))
        disk = pooled = 0
        found = set()
        for z in marks:
            stats = {}
            a = marked_slope_partners(n, x, u, z, {ell: 1}, method="enumerate", stats=stats)
            disk += stats["disk_slots"]
            b = marked_slope_partners(n, x, u, z, {ell: 1}, stats=stats)
            work = stats["pool"]
            pooled += work["quadratic_factors"] + work["block_values"] + work["axis_gcds"]
            assert a == b
            found.update(b)
        assert found == complete and len(marks) == ell - (1 if ell % 4 == 1 else -1)
        observed[ell] = (len(marks), disk, pooled)
    return observed


def check_large_marks():
    d, r = 2**127-1, 162*2**120+1
    a, b = 7523787201323602272, 12598712690840899073
    n, u = d*r, (0, 1, 0)
    v = (d*(b*b-a*a), 0, d*2*a*b)
    assert a*a+b*b == r and certified_prime(r) and certified_prime(d)
    costs = {}
    for fine in (1005*2**120+1, 2**107-1):
        assert certified_prime(fine)
        positive = tuple(i % fine for i in v)
        inv = pow(5, -1, fine)
        negative = (-3*n*inv % fine, 0, 4*n*inv % fine)
        for z, expected in ((positive, (v,)), (negative, ())):
            assert sum(i*i for i in z) % fine == n*n % fine
            assert gcd(n, *z) == 1
            work = {}
            assert marked_slope_partners(n, 0, u, z, {fine: 1}, stats=work) == expected
            if fine > r:
                assert marked_slope_partners(n, 0, u, z, {fine: 1}, method="enumerate") == expected
            tasks = work["pool"]
            costs[fine.bit_length(), bool(expected)] = (
                tasks["method"], tasks["quadratic_factors"], tasks["block_values"])
    return n.bit_length(), costs


if __name__ == "__main__":
    print("Admissible fine marks and occupied marks:", check_small_marks())
    print("Composite positive marks, lost precision factors, prime-square controls:", check_composite_digits())
    print("Single-prime scalar queries and exact point count:", check_scalar())
    print("Full menus (number of marks, disk candidates, pooled factors/blocks/axes):", check_full_menu())
    print("Large positive/negative query controls:", check_large_marks())
