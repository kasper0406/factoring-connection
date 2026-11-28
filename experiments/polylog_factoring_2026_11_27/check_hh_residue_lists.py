"""Check joint residue searches, common-order recovery, and supplied singleton answers."""

from itertools import combinations
from random import Random

from sympy import primerange

if __package__:
    from .check_cofactor_zero_queries import direct_zero
    from .hh_residue_lists import _search_with_base, _shared_prefix_list, factor_from_residue_list, factor_from_shared_bits
else:
    from check_cofactor_zero_queries import direct_zero
    from hh_residue_lists import _search_with_base, _shared_prefix_list, factor_from_residue_list, factor_from_shared_bits


def check_small():
    rng = Random(621004)
    promised = other = transcripts = 0
    for p, q in combinations(primerange(11, 110), 2):
        if q >= 2*p:
            continue
        n = p*q
        for m in (4, 16, 32):
            residues = rng.sample(range(1, m, 2), m//4)
            assert factor_from_residue_list(n, m, residues+[p % m, p % m + m]) in (p, q)
            promised += 1
            other_list = set(residues)-{p % m}
            assert factor_from_residue_list(n, m, other_list) in (None, p, q)
            other += 1
        k = max(2, (2*n.bit_length())//9)
        mask = 1+sum(1 << j for j in range(1, k) if direct_zero(n, 1 << j))
        assert mask == (p & q) % (1 << k)
        assert factor_from_shared_bits(n, k, mask) in (p, q)
        transcripts += 1
    print("True lists, other lists, direct-binomial transcripts:", promised, other, transcripts)


def check_branch_budgets():
    total = 0
    for k in range(2, 8):
        m = 1 << k
        for n in range(1, m, 2):
            by_shared = {}
            for x in range(1, m, 2):
                y = n*pow(x, -1, m) % m
                by_shared.setdefault(x & y, []).append((x, (x^y).bit_count()))
            for shared, states in by_shared.items():
                for budget in range(k):
                    residues, _ = _shared_prefix_list(n, k, shared, budget)
                    assert set(residues) == {r for r, d in states if d <= budget}
                    assert len(residues) <= 1 << budget
                    total += 1
    print("Exact branch-budget lists:", total)


def check_common_order():
    p, q = 12289, 15361
    n, base = p*q, 51095675
    for prime in (p, q):
        assert pow(base, 1024, prime) == 1 and pow(base, 512, prime) != 1
    stats = {}
    factor = _search_with_base(n, 64, [p % 64], 1, 15, base, stats)
    assert factor in (p, q) and stats["status"] == "common_order"
    assert stats["common_order"] == 1024 and stats["affine_intervals"] == 1
    print("Recovered common order:", stats)


def check_larger():
    p, q = 17179869209, 25769803799
    n = p*q
    for residues, m in (([1], 2), ([9, 15, 25, 31, 45, 63], 128)):
        assert p % m in residues
        stats = {}
        assert factor_from_residue_list(n, m, residues, stats=stats) in (p, q)
        assert stats["status"] in ("integer_sum", "collision", "saturated")
        print("Specified list:", m, residues, stats)
    # Separate fixed inputs; the larger transcripts below use the known factors
    # only to furnish the AND bits. The algorithm is given n, k, and the mask.
    for p, q in ((2639587,3890503), (17179869209,25769803799)):
        n = p*q
        for k in sorted({(2*n.bit_length())//9, n.bit_length()//4}):
            shared = (p&q) % (1 << k)
            stats = {}
            assert factor_from_shared_bits(n, k, shared, stats=stats) in (p,q)
            assert stats["branch_budget"] <= ((p^q) % (1 << k)).bit_count()
            print("Supplied AND:", n.bit_length(), k, shared, stats)


if __name__ == "__main__":
    check_small()
    check_branch_budgets()
    check_common_order()
    check_larger()
