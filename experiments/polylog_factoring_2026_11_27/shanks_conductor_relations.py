"""Detect small-conductor combinations of public cyclic cubic splitting fields."""

from itertools import combinations, product
from math import gcd, prod
from random import Random

from sympy import GF, integer_nthroot, primerange
from sympy.polys.matrices import DomainMatrix

from check_shanks_split_seeds import discriminant_root, split_seed


def coprime_basis(integers):
    """Express all inputs as products of powers of pairwise coprime integers."""
    atoms = {n for n in integers if n > 1}
    while True:
        for a, b in combinations(sorted(atoms), 2):
            divisor = gcd(a, b)
            if divisor > 1:
                atoms.difference_update((a, b))
                atoms.update(c for c in (divisor, a//divisor, b//divisor) if c > 1)
                break
        else:
            return sorted(atoms)


def exponent(n, base):
    count = 0
    while n % base == 0:
        n //= base
        count += 1
    return count


def oriented_row(pairs, norms, atom):
    orders = [exponent(n, atom) for n in norms]
    ref = next((pair for pair, v in zip(pairs, orders) if v), None)
    if ref is None:
        return None
    c, d = ref
    return tuple((v if (a*d-b*c) % atom == 0 else -v) % 3
                 for (a,b), v in zip(pairs, orders))


def wild_coordinate(t):
    """Exponent of omega, up to sign, in gamma modulo (1-omega)^3."""
    residue = t % 27
    if residue in (0, 9, 18, 21):
        return 2
    if residue in (3, 6, 15, 24):
        return 1
    return 0


def relation_system(parameters, limit):
    return rational_relation_system([(t,1) for t in parameters],limit)


def rational_wild_coordinate(pair):
    a, b = pair
    return wild_coordinate(a*pow(b,-1,27)) if b % 3 else 0


def rational_relation_system(pairs, limit):
    """Return the rows requiring cancellation and those for the small conductor."""
    if limit < 3:
        raise ValueError("the small-prime support must include 3")
    if any(b <= 0 or gcd(a,b) != 1 for a,b in pairs):
        raise ValueError("expected reduced rational parameters with positive denominator")
    norms = [binary_norm(pair) for pair in pairs]
    tails = norms.copy()
    small = {3: tuple(rational_wild_coordinate(pair) for pair in pairs)}
    for prime in primerange(2, limit+1):
        if prime != 3:
            row = oriented_row(pairs, norms, prime)
            if row:
                small[prime] = row
        for i, tail in enumerate(tails):
            tails[i] //= prime**exponent(tail, prime)

    differences = [abs(a*d-b*c) for (a,b),(c,d) in combinations(pairs, 2)]
    required = []
    for atom in coprime_basis(tails + differences):
        if integer_nthroot(atom, 3)[1]:
            continue
        row = oriented_row(pairs, tails, atom)
        if row and any(row):
            required.append(row)
    return required, small


def kernel_basis(rows, dimension):
    if not rows:
        return [tuple(int(i == j) for j in range(dimension)) for i in range(dimension)]
    matrix = DomainMatrix.from_list(rows, GF(3))
    return [tuple(int(x) % 3 for x in row) for row in matrix.nullspace().to_Matrix().tolist()]


def dot_mod3(row, vector):
    return sum(x*y for x, y in zip(row, vector)) % 3


def supported_conductor(vector, small):
    """The vector must already satisfy all unbounded-prime cancellation rows."""
    return prod((9 if p == 3 else p) for p, row in small.items() if dot_mod3(row, vector))


def bounded_conductor(t, limit):
    """Return the exact conductor if <= limit, None otherwise; no large factorization."""
    if limit < 3:
        return None
    tail = discriminant_root(t)
    conductor = 9 if wild_coordinate(t) else 1
    for p in primerange(2, limit+1):
        v = exponent(tail, p)
        tail //= p**v
        if p != 3 and v % 3:
            conductor *= p
    if conductor <= limit and integer_nthroot(tail, 3)[1]:
        return conductor
    return None


def binary_norm(v):
    a, b = v
    return a*a + 3*a*b + 9*b*b


def shortest_rational_lift(n, t):
    """Gauss reduction in the Eisenstein norm; return (factor, primitive pair)."""
    u, v = (n, 0), (t, 1)
    while True:
        if binary_norm(v) < binary_norm(u):
            u, v = v, u
        a, b = u
        c, d = v
        weight = binary_norm(u)
        twice_inner = 2*a*c + 3*(a*d+b*c) + 18*b*d
        shift = (twice_inner+weight)//(2*weight)
        if shift == 0:
            break
        v = c-shift*a, d-shift*b
    if b < 0:
        a, b = -a, -b
    common = gcd(b, n)
    if common > 1:
        return common, None
    common = gcd(a, b)
    return None, (a//common, b//common)


def rational_trivial(pair):
    if rational_wild_coordinate(pair):
        return False
    delta = binary_norm(pair)
    return bool(integer_nthroot(delta//3**exponent(delta,3), 3)[1])


def eisenstein_product(x, y):
    a, b = x
    c, d = y
    return (a*c-b*d, a*d+b*c-b*d)


def wild_coordinate_direct(t):
    numerator, denominator = (t, -3), (t+3, 3)
    def divide_pi(x):
        a, b = x
        return ((2*a-b)//3, (a+b)//3)
    while sum(numerator) % 3 == 0:
        numerator, denominator = divide_pi(numerator), divide_pi(denominator)
    for j, root in enumerate(((1,0), (0,1), (-1,-1))):
        a, b = eisenstein_product(root, denominator)
        for sign in (-1,1):
            c, d = numerator[0]-sign*a, numerator[1]-sign*b
            if c % 3 == d % 3 == 0 and (c+d) % 9 == 0:
                return j
    raise AssertionError("unexpected norm-one residue modulo pi^3")


def splits_direct(parameters, vector, prime):
    """Power numerator and denominator independently in the finite Eisenstein algebra."""
    def power(x, e):
        out = (1,0)
        while e:
            if e % 2:
                out = tuple(a % prime for a in eisenstein_product(out, x))
            x = tuple(a % prime for a in eisenstein_product(x, x))
            e //= 2
        return out
    numerator = denominator = (1,0)
    for t, e in zip(parameters, vector):
        numerator = eisenstein_product(numerator, power((t,-3),e))
        denominator = eisenstein_product(denominator, power((t+3,3),e))
    degree = 1 if prime % 3 == 1 else 2
    e = (prime**degree-1)//3
    return power(numerator,e) == power(denominator,e)


def check():
    # The two inputs to the constructor are N and the chosen public seeds.
    n = 1247
    seeds = (63, 66)
    parameters = tuple(split_seed(n, x)[1] for x in seeds)
    assert parameters == (-7, 41) and [discriminant_root(t) for t in parameters] == [37, 1813]
    rows, small = relation_system(parameters, 7)
    basis = kernel_basis(rows, 2)
    assert basis == [(1, 1)]
    assert supported_conductor(basis[0], small) == 7
    assert [bounded_conductor(t, 260) for t in parameters] == [37, 259]
    assert [bounded_conductor(t, 7) for t in parameters] == [None, None]
    assert 29*43 == n and (29 % 7, 43 % 7) == (1, 1)
    assert bounded_conductor(1259, 7) == 7  # ignores the unsupplied cube 61^3
    print(f"N={n}, seeds={seeds}, parameters={parameters}, input conductors 37 and 259, "
          f"product conductor {supported_conductor(basis[0], small)}")
    for n, seeds, expected, relation in [(1247,(63,66),(-7,41),(1,1)),
                                         (2813,(78,111),(-44,-7),(2,1))]:
        ts = tuple(split_seed(n,x)[1] for x in seeds)
        assert ts == expected
        rows, small = relation_system(ts, 7)
        assert kernel_basis(rows,2) == [relation] and supported_conductor(relation,small) == 7
        for prime in primerange(2,250):
            if prime not in (3, 7, 37):
                assert splits_direct(ts,relation,prime) == (prime % 7 in (1,6))
    assert all(wild_coordinate_direct(t) == wild_coordinate(t) for t in range(-500,501))
    counts = [0,0,0]
    for n in (713,1247,2813):
        for x in range(n):
            _, t = split_seed(n,x)
            if t is None:
                continue
            divisor, pair = shortest_rational_lift(n,t)
            if divisor:
                assert 1 < divisor < n and n % divisor == 0
                counts[0] += 1
            else:
                a, b = pair
                assert b > 0 and gcd(a,b) == gcd(b,n) == 1 and (a-b*t) % n == 0
                assert binary_norm(pair) <= 3*n
                counts[1 if rational_trivial(pair) else 2] += 1
    print(f"Rational lifts: {counts[0]} gcd splits, {counts[1]} trivial characters, "
          f"{counts[2]} nontrivial characters within the 3N norm bound")

    from sympy import factorint

    # Fully factored, prime-by-prime control is confined to checking.
    rng = Random(217)
    cases = 0
    for _ in range(40):
        pairs = []
        while len(pairs) < 4:
            a, b = rng.randint(-200,200), rng.randint(1,20)
            if gcd(a,b) == 1:
                pairs.append((a,b))
        norms = [binary_norm(pair) for pair in pairs]
        primes = {p for norm in norms for p in factorint(norm)} - {3}
        cutoff = rng.choice([3, 7, 13, 31])
        required, local = rational_relation_system(pairs, cutoff)
        exact = {p: oriented_row(pairs, norms, p) for p in primes}
        basis = kernel_basis(required, len(pairs))
        accepted = 0
        for vector in product(range(3), repeat=len(pairs)):
            large = any(dot_mod3(row, vector) for p, row in exact.items() if p > cutoff)
            assert large == any(dot_mod3(row, vector) for row in required)
            if not large:
                accepted += 1
                tame = prod(p for p, row in exact.items() if dot_mod3(row, vector))
                wild = 9 if dot_mod3(local[3], vector) else 1
                assert supported_conductor(vector, local) == tame*wild
            cases += 1
        assert 3**len(basis) == accepted
    print(f"Checked {cases} combinations against fully factored tame valuation rows")


if __name__ == "__main__":
    check()
