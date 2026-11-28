"""Check the prescribed-index HNF sampler and reductions of power-index lattices."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd, prod
from random import Random

from sympy import Matrix

if __package__:
    from .determinant_lattices import sample_determinant_lattice
else:
    from determinant_lattices import sample_determinant_lattice


def diagonals(n, d):
    if d == 1:
        yield (n,)
        return
    for a in range(1, n+1):
        if n % a == 0:
            yield from ((a, *rest) for rest in diagonals(n//a, d-1))


def entries(diagonal):
    d = len(diagonal)
    slots = [(i, j) for i in range(d) for j in range(i+1, d)]
    for filling in product(*(range(diagonal[i]) for i, _ in slots)):
        rows = [[x if i == j else 0 for j in range(d)] for i, x in enumerate(diagonal)]
        for (i, j), value in zip(slots, filling):
            rows[i][j] = value
        yield tuple(map(tuple, rows))


def powers(n):
    factors = []
    for p in range(2, n+1):
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            factors.append((p, e))
    return factors


def index_count(p, m, d):
    return prod(Fraction(p**(m+j)-1, p**j-1) for j in range(1, d))


class Replay:
    def __init__(self, values):
        self.values = iter(values)

    def randrange(self, n):
        a = next(self.values)
        assert 0 <= a < n
        return a


def gcd_laws(n, d):
    modulus = n*n
    counts, laws = Counter({modulus: 1}), []
    for s in range(1, d):
        following = Counter()
        for previous, weight in counts.items():
            for x in range(modulus):
                following[gcd(previous, x)] += weight
        counts = following
        laws.append({a: Fraction(weight, modulus**s) for a, weight in counts.items()})
    return laws


def check_hnfs(n, d):
    laws = gcd_laws(n, d)
    masses, count, accepted = set(), 0, Fraction()
    for diagonal in diagonals(n, d):
        mass = prod(law[a] for law, a in zip(laws, diagonal[1:]))
        accepted += mass
        multiplicity = prod(a**(d-i-1) for i, a in enumerate(diagonal))
        masses.add(Fraction(mass, multiplicity))
        for rows in entries(diagonal):
            replay = []
            for s, a in enumerate(diagonal[1:], 1):
                replay.append(a % (n*n))
                if a != 1:
                    replay.extend([0]*(s-1))
            replay.extend(rows[i][j] for i in range(d) for j in range(i+1, d) if diagonal[i]>1)
            rng = Replay(replay)
            assert sample_determinant_lattice(n, d, rng) == Matrix(rows)
            assert next(rng.values, None) is None
            count += 1
    factors = powers(n)
    assert len(masses) == 1 and masses.pop()*count == accepted
    assert count == prod(index_count(p, e, d) for p, e in factors)
    assert accepted == prod(1-Fraction(1, p**(e+j)) for p, e in factors for j in range(1, d))
    return count


def rank_mod(rows, p):
    matrix = [[x % p for x in row] for row in rows]
    position = 0
    for j in range(len(rows[0])):
        pivot = next((i for i in range(position, len(matrix)) if matrix[i][j]), None)
        if pivot is None:
            continue
        matrix[pivot], matrix[position] = matrix[position], matrix[pivot]
        inv = pow(matrix[position][j], -1, p)
        for i in range(position+1, len(matrix)):
            multiplier = inv*matrix[i][j] % p
            matrix[i] = [(a-multiplier*b) % p for a, b in zip(matrix[i], matrix[position])]
        position += 1
    return position


def cyclic_probability(p, m, d):
    return prod(1-Fraction(1, p**j) for j in range(2, d+1)) / Fraction(
        prod(1-Fraction(1, p**(m+j)) for j in range(1, d))
    )


def check_power_cohort(n, k, d):
    primes = [p for p, _ in powers(n)]
    counts = Counter()
    first = 0
    split_diagonal = 0
    for diagonal in diagonals(n**k, d):
        multiplicity = prod(a**(d-i-1) for i, a in enumerate(diagonal))
        if diagonal[0] == n**k:
            first += multiplicity
        if any(1 < gcd(a, n) < n for a in diagonal):
            split_diagonal += multiplicity
        for rows in entries(diagonal):
            counts[tuple(rank_mod(rows, p) for p in primes)] += 1
    total = sum(counts.values())
    assert first == Fraction(total)*prod(
        Fraction(p**(k*(d-1)), index_count(p, k, d)) for p in primes
    )
    if k == 1 and len(primes) == 2:
        assert split_diagonal == total-sum(n**j for j in range(d))
    expected = prod(cyclic_probability(p, k, d) for p in primes)
    assert counts[(d-1,)*len(primes)] == expected*total
    for i, p in enumerate(primes):
        local = sum(weight for ranks, weight in counts.items() if ranks[i] == d-1)
        assert Fraction(local, total) == cyclic_probability(p, k, d)
    return counts


def check_layers(rows, p, k):
    """Directly compute the mod-p images of p^-j (L intersect p^j Z^d)."""
    d = len(rows)
    for j in range(k+1):
        divisor = p**j
        modulus = p*divisor
        vectors = set()
        for coeffs in product(range(modulus), repeat=d):
            point = [sum(a*b for a, b in zip(row, coeffs)) % modulus for row in rows]
            if all(x % divisor == 0 for x in point):
                vectors.add(tuple(x//divisor for x in point))
        assert len(vectors) == p**(d-1 if j < k else d)


def run():
    cohorts = [(n, 2) for n in range(1, 21)] + [(n, 3) for n in range(1, 11)]
    cohorts += [(1, 1), (35, 1), (2, 4), (3, 4), (4, 4), (6, 4)]
    count = sum(check_hnfs(n, d) for n, d in cohorts)
    print(f"{len(cohorts)} determinant/dimension pairs; replayed {count} equiprobable HNFs")

    for p, k, d in ((2, 3, 3), (3, 2, 4), (5, 2, 3)):
        frequencies = check_power_cohort(p, k, d)
        print(f"d={d}, determinant {p}^{k}: {dict(sorted(frequencies.items()))}")
    frequencies = check_power_cohort(6, 2, 3)
    split = sum(weight for (r, s), weight in frequencies.items() if r != s)
    print(f"d=3, determinant 36: {split}/{sum(frequencies.values())} mixed mod-2/mod-3 ranks")
    assert split == 1183 and sum(frequencies.values()) == 4550
    assert check_power_cohort(15, 1, 3) == {(2, 2): 403}
    assert sample_determinant_lattice(91, 3, Random(0)).tolist() == [
        [7, 2, 4], [0, 1, 0], [0, 0, 13]
    ]

    rng = Random(61729)
    lattice = sample_determinant_lattice(36, 3, rng)
    while any(rank_mod(lattice.tolist(), p) != 2 for p in (2, 3)):
        lattice = sample_determinant_lattice(36, 3, rng)
    for p in (2, 3):
        check_layers(lattice.tolist(), p, 2)
    for n, d in ((1, 12), (275, 9), (6**64, 16), (2**521-1, 32), ((2**127-1)**3, 40)):
        lattice = sample_determinant_lattice(n, d, rng)
        assert lattice.shape == (d, d) and prod(lattice[i, i] for i in range(d)) == n
        assert all(lattice[i,j] == 0 for i in range(d) for j in range(i))
        assert all(0 <= lattice[i,j] < lattice[i,i] for i in range(d) for j in range(i+1, d))
    print("Checked all six congruence layers of a cyclic sample and five public draws through "
          "dimension 40 and 521-bit determinant.")


if __name__ == "__main__":
    run()
