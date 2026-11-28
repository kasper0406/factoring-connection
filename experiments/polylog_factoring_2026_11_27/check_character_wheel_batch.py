"""Check projected character groups, product evaluation at exact equalities and local coverage."""

from itertools import combinations, product
from math import gcd, prod
from random import Random

from sympy import primerange

if __package__:
    from .cfrac_character_rows import admissible_residues, collect_character_rows
    from .character_wheel_batch import _group_factor, _half_classes, _local_blocks, character_wheel_factor
    from .hh_low_residue import hh_low_factor
else:
    from cfrac_character_rows import admissible_residues, collect_character_rows
    from character_wheel_batch import _group_factor, _half_classes, _local_blocks, character_wheel_factor
    from hh_low_residue import hh_low_factor


def empty_stats():
    return dict(polynomial_coefficients=0, evaluation_points=0, derivative_points=0, saturated=0)


def check_products():
    stats = empty_stats()
    # At 10 the global equality is removed by differentiation. The other roots
    # have differences 7 and 11 there, so the derivative saturates.
    assert _group_factor(77, [10,87,17,21], [10], stats) in (7,11)
    assert stats["derivative_points"] == stats["saturated"] == 1
    assert _group_factor(77, [10,87,18], [10], empty_stats()) is None
    rng = Random(82)
    for _ in range(500):
        p,q = rng.choice(((7,11),(17,19),(31,47)))
        n = p*q
        A,B = ([rng.randrange(2*n) for _ in range(rng.randrange(1,15))] for _ in range(2))
        direct = {gcd(a-b,n) for a in A for b in B} & {p,q}
        got = _group_factor(n,A,B,empty_stats())
        assert got in direct if direct else got is None
    print("500 randomized product checks and two exact-equality controls")


def check_projections():
    primes = (3,5,7)
    count = 0
    blocks = _local_blocks(primes)
    for p,q in combinations(primerange(11,105),2):
        if q >= 2*p:
            continue
        n = p*q
        for budget in (1,6,22):
            constraints = collect_character_rows(n,primes,budget)
            if constraints.divisor:
                continue
            rows = [mask for mask,_ in constraints.rows]
            expected = sorted(admissible_residues(constraints))
            m = constraints.conductor
            # Compare projections and all matching right halves under two distinct cuts.
            for cut in (1,2):
                L = _half_classes(blocks[:cut],rows,m)
                R = _half_classes(blocks[cut:],rows,m)
                actual = sorted((a+b)%m for (s,A),(t,B) in zip(L,R,strict=True) if s==t
                                for a in A for b in B)
                assert actual == expected
            assert character_wheel_factor(n,constraints) in (p,q)
            count += 1
    print(f"{count} complete projected class sets matched explicit CRT and split each input")


def check_lock_cover():
    # Prescribe the positive quadratic character at 5,13,17, independent of N.
    primes = (5,13,17)
    n = 19*31
    m = prod(primes)
    roots = [sorted({pow(i,2,p) for i in range(1,p)}) for p in primes]
    classes = {r for r in range(m) if all(r%p in good for p,good in zip(primes,roots))}
    covered = set()
    idempotents = [(m//p)*pow(m//p,-1,p) for p in primes]
    pairs = [list(zip(good[::2],good[1::2])) for good in roots]
    count = 0
    for selected in product(*pairs):
        a = sum(e*u*v*pow(n,-1,p) for e,p,(u,v) in zip(idempotents,primes,selected))%m
        c = sum(e*(u+v) for e,(u,v) in zip(idempotents,selected))%m
        local = {r for r in classes if (r*r-c*r+a*n)%m == 0}
        assert len(local)==2**len(primes) and not covered & local
        covered.update(local)
        count += 1
    assert covered==classes and count==12
    print(f"{count} primitive linear targets attain the exact root-count bound on {len(classes)} classes")


def check_benchmarks():
    for n, primes, budget, rank in ((762201824107,(3,5,7,11,13,17),1600,4),
                                    (13069950186793286509,(3,5,7),4096,2)):
        rows_work = {}
        constraints = collect_character_rows(n,primes,budget,stop_rank=rank,stats=rows_work)
        assert not constraints.divisor and len(constraints.rows)==rank
        work = {}
        factor = character_wheel_factor(n,constraints,stats=work)
        assert 1 < factor < n and n%factor == 0
        if n==762201824107:
            assert constraints.class_count==23040
            assert (work["projected_left"],work["projected_right"])==(192,120)
            assert (work["polynomial_coefficients"],work["evaluation_points"])==(193,240)
            baseline = {}
            assert n % hh_low_factor(n,2,1,stats=baseline) == 0
            assert (baseline["baby_count"],baseline["giant_count"])==(159,218)
        print(f"{n.bit_length()} bits: {rows_work['steps']} input-only convergents, "
              f"{constraints.class_count} compatible classes, {work['projected_left']} and "
              f"{work['projected_right']} half classes, {work['polynomial_coefficients']} "
              f"polynomial coefficients and {work['evaluation_points']} evaluated points")


if __name__=="__main__":
    check_products()
    check_projections()
    check_lock_cover()
    check_benchmarks()
