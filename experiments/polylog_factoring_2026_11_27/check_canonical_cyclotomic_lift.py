"""Seeded control for factoring from Φ_n(x) mod n² at canonical units 1 ≤ x < n."""

from collections import Counter
from itertools import combinations
from math import gcd, isqrt
from random import Random

from sympy import Poly, Symbol, cyclotomic_poly

if __package__:
    from .check_cyclotomic_at_index import geometric_sum
else:
    from check_cyclotomic_at_index import geometric_sum


class Nonunit(Exception):
    """A denominator whose gcd either factors n or makes the sample singular."""


class InvalidResidue(Exception):
    """A generic Φ_n(x) residue must reduce to 1 modulo n."""


def inverse(a, n):
    d = gcd(a, n)
    if d != 1:
        raise Nonunit(d)
    return pow(a, -1, n)


def uniform_unit(n, rng):
    while True:
        a = rng.randrange(1, n)
        if gcd(a, n) == 1:
            return a


def canonical_row(n, u, z, query):
    w = u * pow(z, n, n*n) % (n*n)
    carry, x = divmod(w, n)
    h = carry * inverse(x, n) % n
    v = pow(x, n, n)
    gamma = x * inverse(x - 1, n) % n
    alpha = v * inverse(v - 1, n) % n
    alpha_inv = inverse(alpha, n)
    residue = query(x)
    if residue % n != 1:
        raise InvalidResidue
    r = gamma * alpha_inv % n
    y = ((residue - 1) // n * alpha_inv + h*r) % n
    return r, h, y


def candidate_sum(n, rows):
    if len(rows) == 1:  # u = 1, so the nuisance coefficient θ vanishes.
        _, h, y = rows[0]
        return y * inverse(h, n) % n
    differences = [(r1-r0, h1-h0, y1-y0)
                   for (r0, h0, y0), (r1, h1, y1) in combinations(rows, 2)]
    for dr, _, _ in differences:
        if dr % n:
            inverse(dr, n)  # A coincidence in just one CRT component splits n.
    if all(dr % n == 0 for dr, _, _ in differences):
        for _, dh, dy in differences:
            if dh % n:
                return dy * inverse(dh, n) % n
        raise Nonunit(n)
    dr1, dh1, dy1 = differences[0]  # row 1 minus row 0
    dr2, dh2, dy2 = differences[1]  # row 2 minus row 0
    return (dr1*dy2 - dr2*dy1) * inverse(dr1*dh2 - dr2*dh1, n) % n


def recover(n, rng, query, principal):
    """The sampler and recovery receive only n, randomness, and a canonical oracle."""
    u = 1 if principal else uniform_unit(n*n, rng)
    try:
        rows = [canonical_row(n, u, uniform_unit(n, rng), query)
                for _ in range(1 if principal else 3)]
        s = candidate_sum(n, rows)
    except InvalidResidue:
        return "rejected", None
    except Nonunit as exc:
        d = exc.args[0]
        return ("gcd", d) if d != n else ("singular", None)
    delta = s*s - 4*n
    if delta >= 0:
        root = isqrt(delta)
        a, b = (s-root)//2, (s+root)//2
        if root*root == delta and 1 < a <= b and a*b == n:
            return "sum", a
    return "rejected", None


def factored_value(p, q, x):
    n = p*q
    assert 1 <= x < n and gcd(pow(x, n, n)-1, n) == 1
    modulus = n*n
    numerator = geometric_sum(pow(x, p, modulus), q, modulus)
    denominator = geometric_sum(x, q, modulus)
    return numerator * pow(denominator, -1, modulus) % modulus


def check_exact_values():
    X = Symbol("X")
    checked = 0
    for p, q in ((3, 5), (3, 7), (5, 11), (7, 13)):
        n = p*q
        polynomial = Poly(cyclotomic_poly(n, X), X)
        for x in range(1, n):
            if gcd(x, n) == gcd(pow(x, n, n) - 1, n) == 1:
                assert factored_value(p, q, x) == int(polynomial.eval(x)) % (n*n)
                checked += 1
    return checked


def check_equal_ratios():
    # All three ratios agree in both components, including when θ is nonzero.
    p, q = 11, 67
    n = p*q
    rows = [canonical_row(n, 2, z, lambda x: factored_value(p, q, x))
            for z in (694, 716, 729)]
    assert [r for r, _, _ in rows] == [679, 679, 679]
    assert candidate_sum(n, rows) == p+q


def trials(p, q, principal, mode, count=128):
    n = p*q
    # Reset for matched exact/noisy rows. A separate stream fixes errors by input.
    draws = Random(271828 + n + int(principal))
    errors = Random(314159 + n)
    offsets = {}
    stats = Counter()
    first_sum = None
    for trial in range(count):
        bases = []

        def query(x):
            bases.append(x)
            stats["queries"] += 1
            if mode == "miss-principal":
                # For n=11·67 the principal image modulo n is unrestricted modulo 11
                # and the subgroup of order 6 modulo 67. Give wrong low digits there.
                assert (p, q) == (11, 67)
                if pow(x, 6, q) == 1:
                    stats["corrupt"] += 1
                    return 0
            value = factored_value(p, q, x)
            if mode == "noisy":
                if x not in offsets:
                    offsets[x] = errors.randrange(1, n) if errors.randrange(4) == 0 else 0
                stats["corrupt"] += int(offsets[x] != 0)
                value = (value + n*offsets[x]) % (n*n)
            return value

        kind, d = recover(n, draws, query, principal)
        stats[kind] += 1
        if d is not None:
            assert d in (p, q)
        if kind == "sum" and first_sum is None:
            first_sum = trial, tuple(bases), d + n//d
    assert sum(stats[key] for key in ("sum", "gcd", "singular", "rejected")) == count
    if mode == "miss-principal" and principal:
        assert not stats["sum"] and stats["queries"] == stats["corrupt"] > 0
    else:
        assert stats["sum"] > 0
    if mode == "exact":
        assert not stats["rejected"]
    return stats, first_sum


if __name__ == "__main__":
    print("exact polynomial checks", check_exact_values())
    check_equal_ratios()
    print("p q orbit oracle trials sum gcd singular rejected queries corrupt")
    for p, q in ((101, 127), (1009, 1013), (11, 67)):
        modes = ("exact", "noisy", "miss-principal") if (p, q) == (11, 67) else ("exact", "noisy")
        for principal in (False, True):
            for mode in modes:
                stats, first = trials(p, q, principal, mode)
                print(p, q, "u=1" if principal else "random", mode, 128,
                      *(stats[k] for k in ("sum", "gcd", "singular", "rejected", "queries", "corrupt")))
                if (q-1) % p == 0 and not principal and mode != "noisy":
                    print("non-surjective sum (trial from zero, canonical bases, decoded sum)", first)
