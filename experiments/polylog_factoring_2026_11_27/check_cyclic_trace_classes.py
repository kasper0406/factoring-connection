"""Check projected dyadic residue classes and mixed-conductor order recovery."""

from math import gcd
from random import Random

from sympy import nextprime, primerange

if __package__:
    from .cyclic_trace_classes import trace_shortlist
    from .hh_residue_lists import _search_with_base, factor_from_residue_list
else:
    from cyclic_trace_classes import trace_shortlist
    from hh_residue_lists import _search_with_base, factor_from_residue_list


def check_local():
    checked = 0
    for prime in primerange(3, 2000):
        values = list(range(prime))
        for level in range(1, 9):
            values = [(v*v-2) % prime for v in values]
            modulus = 1 << (level+2)
            roots = (1 << level) if prime % modulus in (1,modulus-1) else 0
            assert values.count(0) == roots
            checked += 1
    print("Prime/level root counts:", checked)


def check_inputs():
    fixtures = ((17,31,5,2), (2207,3167,18,3), (257,383,15,4), (5662847,6803327,18,5))
    for p,q,x,j in fixtures:
        stats = {}
        assert trace_shortlist(p*q, (x,), j,j,stats=stats) in (p,q)
        assert stats["status"] == "two_classes"
    rng = Random(70395)
    for size in (16,24,32):
        statuses = {"empty":0, "trace_gcd":0, "two_classes":0}
        for _ in range(64):
            p,q=sorted(int(nextprime(rng.randrange(1 << (size-1),1 << size))) for _ in range(2))
            if p == q:
                continue
            stats = {}
            factor = trace_shortlist(p*q, range(3,35), 1,size,stats=stats)
            assert factor in (None,p,q)
            statuses[stats["status"]] += 1
        print("Independent inputs:",size,statuses)


def check_moduli():
    p,q = 12289,15361
    for base, expected in ((133104511,"common_order"),(162032817,"component_order")):
        work = {}
        assert _search_with_base(p*q,192,[1],1,15,base,work,(2,3)) in (p,q)
        assert work["status"] == expected
    checked = 0
    for p in primerange(37,190):
        for q in primerange(p+1,2*p):
            for m, primes in ((3,(3,)),(9,(3,)),(30,(2,3,5)),(45,(3,5)),(140,(2,5,7))):
                rs = [r for r in (p%m,1) if gcd(r,m)==1]
                assert factor_from_residue_list(p*q,m,rs,modulus_primes=primes) in (p,q)
                checked += 1
    print("Mixed-modulus promised lists:", checked)


if __name__ == "__main__":
    check_local()
    check_inputs()
    check_moduli()
