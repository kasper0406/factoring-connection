"""Verify CFRAC sign certificates, exact CRT class counts and fixed-budget coverage."""

from collections import Counter
from itertools import combinations
from math import gcd, prod
from random import Random

from sympy import nextprime, primerange

if __package__:
    from .cfrac_character_rows import admissible_residues, character_word, collect_character_rows
    from .hh_residue_lists import factor_from_residue_list
else:
    from cfrac_character_rows import admissible_residues, character_word, collect_character_rows
    from hh_residue_lists import factor_from_residue_list


BASE = (3,5,7)


def verify_certificates(result, n, p, q):
    signed = [-1, 2]+[ell if ell % 4 == 1 else -ell for ell in result.primes]
    words = [character_word(x, result.primes) for x in (p,q,n)]
    assert words[0] ^ words[1] == words[2]
    for mask, root in result.rows:
        radicand = prod(d for i,d in enumerate(signed) if mask & (1 << i))
        assert gcd(root,n) == 1 and root*root % n == radicand % n
        assert all((mask & w).bit_count() % 2 == 0 for w in words)


def check_small():
    count = 0
    for p,q in combinations(primerange(11,115),2):
        if q >= 2*p:
            continue
        n = p*q
        for limit in (1,6,30):
            result = collect_character_rows(n, BASE, limit)
            verify_certificates(result, n, p, q)
            if result.divisor:
                assert result.divisor in (p,q)
                continue
            classes = sorted(admissible_residues(result))
            direct = [u for u in range(result.conductor) if gcd(u, result.conductor) == 1
                      and all((mask & character_word(u, BASE)).bit_count() % 2 == 0
                              for mask,_ in result.rows)]
            assert classes == direct and len(classes) == result.class_count
            assert p % result.conductor in classes and q % result.conductor in classes
            count += 1
    print(f"{count} small unsplit relation matrices agreed with direct CRT enumeration")


def balanced_pairs(bits):
    rng = Random(20261129+bits)
    low, high = 3*(1 << (bits//2-2)), (1 << (bits//2))-128
    yielded = 0
    while yielded < 64:
        p,q = sorted(int(nextprime(rng.randrange(low, high))) for _ in range(2))
        if p != q:
            assert (p*q).bit_length() == bits and q < 2*p
            yield p,q
            yielded += 1


def check_cohorts():
    for bits, expected in ((32,(0,0,26,38,6886)), (48,(12,9,19,24,85879)),
                           (64,(48,12,2,2,252308)), (80,(64,0,0,0,409600))):
        counts = Counter()
        steps = 0
        for p,q in balanced_pairs(bits):
            n = p*q
            stats = {}
            result = collect_character_rows(n, BASE, bits*bits, stop_rank=2, stats=stats)
            verify_certificates(result, n, p, q)
            if result.divisor:
                assert result.divisor in (p,q)
                counts["split"] += 1
            else:
                counts[len(result.rows)] += 1
            steps += stats["steps"]
        record = counts[0], counts[1], counts[2], counts["split"], steps
        assert record == expected
        print(f"{bits} bits: ranks 0,1,2 = {record[:3]}, direct splits {record[3]}, "
              f"continued-fraction steps {steps}")


def check_completion():
    n = 13069950186793286509
    stats = {}
    result = collect_character_rows(n, BASE, n.bit_length()**2, stop_rank=2, stats=stats)
    assert result.divisor is None and stats["steps"] == 2465 and stats["dependencies"] == 0
    assert [mask for mask,_ in result.rows] == [5,9]
    for (mask,root), d in zip(result.rows,(3,-5)):
        assert root*root % n == d % n
    lists = (sorted(admissible_residues(result)),
             [u for u in range(840) if gcd(u,840) == 1])
    for candidates, expected in zip(lists, ((48,1189,34),(192,2559,21))):
        work = {}
        divisor = factor_from_residue_list(n, result.conductor, candidates,
                                          modulus_primes=(2,*BASE), stats=work)
        assert 1 < divisor < n and n % divisor == 0
        observed = len(candidates), work["baby_count"],work["giant_count"]
        assert observed == expected
        print(f"N={n}: {observed[0]} classes, {observed[1]} babies, {observed[2]} giants, "
              f"verified divisor {divisor}")


if __name__ == "__main__":
    check_small()
    check_cohorts()
    check_completion()
