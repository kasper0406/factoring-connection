"""Exact controls for stabilised cyclotomic values and the first canonical digit."""

from collections import Counter, defaultdict
from functools import cache
from itertools import combinations
from math import gcd

from sympy import Poly, Symbol, cyclotomic_poly, factorint, n_order, primerange

if __package__:
    from .stationary_dual import prime_power_base
else:
    from stationary_dual import prime_power_base


def evaluate_power_lift(n, a, precision):
    """Construct b and return (b, Φ_n(b) mod n^precision) from public inputs."""
    if n <= 1 or precision < 1:
        raise ValueError("use an index greater than one and positive precision")
    modulus = n ** precision
    depth = max(precision, n.bit_length())
    b = pow(a, n ** depth, modulus)
    if gcd(a - 1, n) == 1:
        return b, 1
    prime = prime_power_base(n)
    value = prime if prime is not None and a % prime == 1 else 1
    return b, value % modulus


@cache
def coefficients(n):
    x = Symbol("x")
    return tuple(map(int, Poly(cyclotomic_poly(n, x), x).all_coeffs()))


def horner(coeffs, x, modulus):
    value = 0
    for coefficient in coeffs:
        value = (value*x + coefficient) % modulus
    return value


def check_stabilised():
    counts = Counter()
    for n in range(2, 151):
        factors = factorint(n)
        for k in (1, 2, 3):
            modulus = n**k
            for a in (-3, -1, 0, 1, 2, 3, 7):
                for mode, depth in (("precision", k), ("full", max(k, n.bit_length()))):
                    z = pow(a, n**depth, modulus)
                    value = horner(coefficients(n), z, modulus)
                    for p, e in factors.items():
                        local_modulus = p**(e*k)
                        assert pow(z, p, local_modulus) == z % local_modulus
                        exceptional = z % p != 0 and n_order(z % p, p) == n // p**e
                        expected = p if exceptional else 1
                        assert value % local_modulus == expected % local_modulus
                        counts[f"{mode}_local"] += 1
                        if exceptional and len(factors) > 1:
                            counts[f"{mode}_mixed_exception"] += 1
                    if mode == "full" and len(factors) > 1:
                        assert value == 1
                    if mode == "full":
                        assert evaluate_power_lift(n, a, k) == (z, value)
                    counts[mode] += 1
    return dict(counts)


def local_table(n, r, s):
    """Map a generic residue to the canonical digit at that lift and its slope."""
    table = {}
    roots = {b for b in range(1, r)
             if ((s-1)*pow(b, s, r) - s*pow(b, s-1, r) + 1) % r == 0}
    ratios = [u*pow(v, -1, r) % r for u in roots for v in roots if u != v]
    assert len(set(ratios)) == len(ratios) <= r-2
    for b in range(1, r):
        b_to_n = pow(b, n, r)
        if b_to_n == 1:
            continue
        value = horner(coefficients(n), b, r*r)
        assert value % r == 1
        digit = (value-1)//r * pow(s, -1, r) % r
        derivative = (pow(b-1, -1, r)
                      - s*pow(b, n-1, r)*pow(b_to_n-1, -1, r)) % r
        assert (derivative == 0) == (b in roots)
        table[b] = digit, derivative * pow(s, -1, r) % r
    return table, len(ratios)


def check_canonical():
    counts = Counter()
    for p, q in combinations(primerange(3, 38), 2):
        if q >= 2*p:
            continue
        n = p*q
        p_table, p_ratios = local_table(n, p, q)
        q_table, q_ratios = local_table(n, q, p)
        counts["balanced_moduli"] += 1
        counts["distinct_ratios"] += p_ratios + q_ratios
        histograms = {p: defaultdict(Counter), q: defaultdict(Counter)}
        local_predictions = {p: {}, q: {}}
        cells = []
        for x in range(1, n):
            if gcd(x, n) != 1 or gcd(pow(x, n, n)-1, n) != 1:
                continue
            value = horner(coefficients(n), x, n*n)
            assert value % n == 1
            digit = (value-1)//n
            cells.append((x, digit))
            for r, table in ((p, p_table), (q, q_table)):
                b = x % r
                initial, slope = table[b]
                assert digit % r == (initial + slope*((x-b)//r)) % r
                histograms[r][b][digit % r] += 1
                counts["local_taylor"] += 1
            if p_table[x % p][1] == 0 and q_table[x % q][1] == 0:
                counts["joint_critical"] += 1
        # On each noncritical fibre even a factor-informed local function
        # can agree on at most ceil(s/r) canonical lifts.
        bound = 0
        for r, s, table in ((p, q, p_table), (q, p, q_table)):
            for b, (_, slope) in table.items():
                prediction, multiplicity = histograms[r][b].most_common(1)[0]
                local_predictions[r][b] = prediction
                if slope:
                    assert multiplicity <= (s+r-1)//r
                    bound += multiplicity
                    counts["noncritical_fibres"] += 1
        critical_p = {b for b, (_, slope) in p_table.items() if slope == 0}
        critical_q = {b for b, (_, slope) in q_table.items() if slope == 0}
        bound += len(critical_p)*len(critical_q)
        if (p, q) == (11, 17):
            assert critical_p == {6, 9} and critical_q == {3, 7}
        hits = sum(all(local_predictions[r][x % r] == digit % r for r in (p, q))
                   for x, digit in cells)
        assert hits <= bound
        counts["marginal_mode_hits"] += hits
        counts["generic_canonical"] += len(cells)
    return dict(counts)


if __name__ == "__main__":
    print("Stabilised values", check_stabilised())
    print("Canonical digits", check_canonical())
