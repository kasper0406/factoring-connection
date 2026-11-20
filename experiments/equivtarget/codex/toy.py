"""Diagnostics specified in prereg.txt. Run with Python and SymPy."""
import math
import random
from collections import Counter, defaultdict
from fractions import Fraction
from sympy import isprime, primerange, sqrt_mod

SEED = 20261120
SMALL = list(primerange(2, 98))
MULTIPLIERS = list(primerange(101, 200))
PAIRS = [(x, y) for x in range(1, 33) for y in range(1, 33)
         if math.gcd(x, y) == 1]


def random_prime(rng, lo, hi):
    while True:
        x = rng.randrange(lo | 1, hi + 1, 2)
        if isprime(x):
            return x


def modulus(rng, bits):
    lo, hi = 1 << (bits - 1), (1 << bits) - 1
    p = random_prime(rng, lo, hi)
    q = random_prime(rng, lo, hi)
    while p == q:
        q = random_prime(rng, lo, hi)
    return p * q


def score(d):
    total = 0
    for ell in SMALL[1:]:
        symbol = pow(d, (ell - 1) // 2, ell)
        total += (symbol if symbol <= 1 else -1) * math.log(ell) / ell
    return total


def is_smooth(value):
    for ell in SMALL:
        while value % ell == 0:
            value //= ell
    return value == 1


def form_yield(rng, m):
    scale = math.isqrt(m)
    while True:
        a = random_prime(rng, 55 * scale // 100, 75 * scale // 100)
        root = sqrt_mod(-m, a)
        if root is not None and root:
            break
    b = min(root, a - root) * rng.choice((-1, 1))
    c = (b*b + m) // a
    assert a*c - b*b == m and abs(2*b) <= a <= c
    return sum(is_smooth(a*x*x + 2*b*x*y + c*y*y) for x, y in PAIRS)


def forms():
    rng = random.Random(SEED)
    for bits in (16, 24):
        counts = [0, 0]
        wins = ties = 0
        for _ in range(24):
            n = modulus(rng, bits)
            best = max(MULTIPLIERS, key=lambda k: score(-k*n))
            other = rng.choice([k for k in MULTIPLIERS if k != best])
            hits = [sum(form_yield(rng, k*n) for _ in range(4))
                    for k in (best, other)]
            counts = [x+y for x, y in zip(counts, hits)]
            wins += hits[0] > hits[1]
            ties += hits[0] == hits[1]
        print("class", bits, "norms per arm", 24*4*len(PAIRS),
              "scored, control", *counts, "paired wins, ties", wins, ties)


def quotient_j(a, b, n):
    # Source shifted to y²=x³+a*x²+b*x. Quotient by (0,0)
    # has coefficients A=-2a and B=a²-4b in the same model.
    A, B = -2*a, a*a-4*b
    numer = 256 * pow(A*A-3*B, 3, n)
    denom = B*B*(A*A-4*B) % n
    assert math.gcd(denom, n) == 1
    return numer * pow(denom, -1, n) % n


def legendre_neighbours(lam, n):
    roots = (0, 1, lam)
    js = []
    for i, t in enumerate(roots):
        u, v = [val for j, val in enumerate(roots) if i != j]
        js.append(quotient_j(2*t-u-v, (t-u)*(t-v), n))
    j = 256 * pow(1-lam+lam*lam, 3, n) * pow(lam*lam*(lam-1)**2, -1, n) % n
    return j, js


def edges():
    rng = random.Random(SEED)
    for bits in (16, 24):
        bad = collisions = clean_splits = 0
        for _ in range(2000):
            n = modulus(rng, bits)
            lam = rng.randrange(n)
            g = math.gcd(lam*(lam-1), n)
            if g != 1:
                bad += (g < n or any(1 < math.gcd(x, n) < n for x in (lam, lam-1)))
                continue
            j, js = legendre_neighbours(lam, n)
            # Independently verify each quotient with the level-two modular polynomial.
            assert all((j**3+h**3-j*j*h*h + 1488*j*h*(j+h)
                        -162000*(j*j+h*h)+40773375*j*h
                        +8748000000*(j+h)-157464000000000) % n == 0 for h in js)
            diff = [math.gcd(js[0]-v, n) for v in js[1:]]
            split = any(1 < g < n for g in diff)
            delta = math.prod(js[i]-js[j] for i in range(3) for j in range(i))
            if math.gcd(delta, n) != 1:
                collisions += split
            else:
                clean_splits += split
        print("edges", bits, "bad-input, collision, other splits", bad, collisions, clean_splits)


def cm():
    for b in (16, 20, 24):
        lo, hi = 2**b, 2**b + 2**15
        primes = set(primerange(lo, hi))
        found = set()
        counts = []
        for d in range(1, 1024):
            if d % 8 == 3 and all(d % (i*i) for i in range(2, math.isqrt(d)+1)):
                for s in range(math.isqrt((4*lo-1)//d), math.isqrt((4*hi-1)//d)+1):
                    val = d*s*s + 1
                    if val % 4 == 0 and val//4 in primes:
                        found.add(val//4)
            if d in (63, 255, 1023):
                counts.append(len(found))
        print("CM", b, "primes", len(primes), "hits", *counts)


def masked_fibres():
    for primes in ((101, 103), (211, 223), (307, 311)):
        counts = []
        fibre_counts = []
        for p in primes:
            fibres = defaultdict(list)
            for lam in range(2, p):
                j, neighbours = legendre_neighbours(lam, p)
                if j not in (0, 1728 % p) and len(set(neighbours)) == 3:
                    fibres[j].append((neighbours[0], min(neighbours)))
            match = total = 0
            for observations in fibres.values():
                assert len(observations) == 6
                assert sorted(Counter(h for h, _ in observations).values()) == [2, 2, 2]
                assert len({selected for _, selected in observations}) == 1
                match += sum(hidden == selected for hidden, selected in observations)
                total += len(observations)
            counts.append((match, total))
            fibre_counts.append(len(fibres))
        (m, t), (n, u) = counts
        split = m*(u-n) + n*(t-m)
        assert Fraction(split, t*u) == Fraction(4, 9)
        print("private edges", primes, "j fibres", *fibre_counts, "splits/pairs", split, t*u)


if __name__ == "__main__":
    forms()
    edges()
    cm()
    masked_fibres()
