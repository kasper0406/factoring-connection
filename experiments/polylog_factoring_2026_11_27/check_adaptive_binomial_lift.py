"""Check multiplier alignment, central-binomial-only lifting, and local singleton ambiguity."""

from itertools import combinations
from math import comb
from random import Random

from sympy import isprime, nextprime, primerange

if __package__:
    from .adaptive_binomial_lift import aligned_multiplier, factor_by_single_binomials
else:
    from adaptive_binomial_lift import aligned_multiplier, factor_by_single_binomials


def check_selector():
    count = 0
    for d in (8, 16, 32, 64, 128, 256):
        for a in range(1, d, 2):
            for b in range(1, d, 2):
                h = aligned_multiplier(a, b, d)
                exceptional = (a+b) % d == 0 or (a-b) % d == d//2
                assert (h is None) == exceptional
                if exceptional:
                    assert all((u*a % d < d//2) != (u*b % d < d//2) for u in range(1, d, 2))
                else:
                    assert 0 < h < d//4 and h % 2
                    assert (h*a % d < d//2) == (h*b % d < d//2)
                count += 1
    return count


def _carries(t, p):
    while t:
        if 2*(t % p) >= p:
            return True
        t //= p
    return False


def check_small():
    cases = calls = unbalanced = singular = 0
    primes = list(primerange(3, 300))
    largest_states = 0
    balanced_cases = shifted_calls = 0
    for p, q in combinations(primes, 2):
        n = p*q
        if n > 17000:
            continue
        def exact(modulus, t):
            nonlocal calls
            assert modulus == n
            answer = comb(2*t, t) % n == 0
            assert answer == (_carries(t, p) and _carries(t, q))
            calls += 1
            return answer
        stats = {}
        assert factor_by_single_binomials(n, exact, stats=stats) in (p, q)
        cases += 1
        unbalanced += q > 2*p
        singular += stats["opposite"]
        largest_states = max(largest_states, stats["peak_states"])
        if q < 2*p:
            assert factor_by_single_binomials(n, exact, balanced=True, stats=stats) in (p, q)
            assert stats["queries"] <= 4*stats["levels"] and stats["peak_states"] <= 4
            balanced_cases += 1
            shifted_calls += stats["shifted"]
    return cases, unbalanced, calls, singular, largest_states, balanced_cases, shifted_calls


def check_large():
    rng = Random(9302)
    fixtures = []
    for bits in (24, 40, 64, 96, 128):
        p = int(nextprime(2**(bits-1)+rng.randrange(2**(bits-2))))
        q = int(nextprime(p+rng.randrange(1, 2**(bits-2))))
        fixtures.append((p, q))
        fixtures.append((p, int(nextprime(2**(bits*2)+rng.randrange(2**bits)))))
    # Force many consecutive opposite-prefix stages.
    p = int(nextprime(3*2**45))
    while not isprime(2**48-p):
        p = int(nextprime(p))
    fixtures.append((p, 2**48-p))
    total = 0
    for p, q in fixtures:
        n = p*q
        stats = {}
        def oracle(modulus, t):
            assert modulus == n
            return _carries(t, p) and _carries(t, q)
        assert factor_by_single_binomials(n, oracle, stats=stats) in (p, q)
        levels = stats["levels"]
        assert stats["peak_states"] <= 2+4*levels
        assert stats["queries"] <= 2*levels+2*levels*(levels-1)
        if n % 8 != 7:
            assert stats["opposite"] == 0 and stats["peak_states"] <= 2
        total += stats["queries"]
    old = stats.copy()
    p, q = fixtures[-1]
    assert q < 2*p
    assert factor_by_single_binomials(p*q, oracle, balanced=True, stats=stats) in (p, q)
    assert stats["peak_states"] <= 4 and 0 < stats["shifted"] <= stats["queries"]
    return len(fixtures), total, old, stats


def check_affine_finish():
    results = []
    for p, q in ((4001, 5003), (60013, 88001), (802007, 1364009)):
        n = p*q
        assert isprime(p) and isprime(q) and p < q < 2*p
        calls = []
        def oracle(modulus, t):
            assert modulus == n
            calls.append(t)
            return _carries(t, p) and _carries(t, q)
        ordinary = {}
        assert factor_by_single_binomials(n, oracle, balanced=True, stats=ordinary) in (p, q)
        full_calls, calls = calls, []
        early = {}
        assert factor_by_single_binomials(n, oracle, balanced=True,
                                          affine_finish=True, stats=early) in (p, q)
        assert calls == full_calls[:len(calls)] and len(calls) < len(full_calls)
        assert early["levels"] < ordinary["levels"]
        assert early["completion_modulus"] == 4 << early["levels"]
        assert 1 <= early["affine_calls"] <= 2*early["peak_states"]
        assert early["maximum_dimension"] >= 3
        results.append((n, len(full_calls), len(calls), early["completion_modulus"],
                        early["affine_calls"], early["maximum_dimension"]))
    return results


def check_fibers():
    # Size-biased statistics: choose an actual pair before inspecting its cell.
    results = {}
    for k in (5, 7, 9, 11):
        m = 1 << k
        cells = {}
        for a in range(1, m, 2):
            for b in range(1, m, 2):
                key = a*b % m, a & b
                cells[key] = cells.get(key, 0)+1
        total = (m//2)**2
        assert sum(cells.values()) == total
        running = 0
        for width in sorted(cells.values()):
            running += width
            if running*2 >= total:
                median = width
                break
        results[k] = median, max(cells.values())
    return results


if __name__ == "__main__":
    print("Exhaustive multiplier pairs:", check_selector())
    print("Direct exact (cases, unbalanced, Boolean calls, opposite states, peak, balanced reruns, shifted states):", check_small())
    print("Factor-supplied controls (cases, calls, power-of-two-sum fixture):", check_large())
    print("Affine completion (N, full calls, early calls, modulus, windows, dimension):", check_affine_finish())
    print("Singleton low-prefix cells (pair-weighted median, maximum):", check_fibers())
