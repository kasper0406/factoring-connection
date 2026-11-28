"""Public ternary Gram challenges for odd distinct semiprimes and exact controls."""

from collections import Counter
from itertools import product
from math import gcd, isqrt
from random import Random

from sympy import jacobi_symbol

if __package__:
    from .gram_embedding import embed_gram
else:
    from gram_embedding import embed_gram


def _gram(n, x):
    return x*x + 1, -n*x, n*n


def sample_mixed_gram(n, rng):
    """Return (x,Gram) from the first uniform x with Jacobi(−x²−1,n)=−1.

    The expected-cost bound is for n=pq, p and q distinct odd primes.
    """
    if n <= 1 or n % 2 == 0 or isqrt(n)**2 == n:
        raise ValueError("expected an odd nonsquare n > 1")
    while True:
        x = rng.randrange(n)
        if jacobi_symbol(-(x*x + 1), n) == -1:
            return x, _gram(n, x)


def factor_from_frame(gram, u, v):
    """Verify the marked ternary Gram and return a proper gcd of v, if any."""
    if len(u) != 3 or len(v) != 3:
        raise ValueError("expected two ternary vectors")
    a, b, c = gram
    if (sum(x*x for x in u), sum(x*y for x, y in zip(u, v)), sum(y*y for y in v)) != (a, b, c):
        raise ValueError("the frame does not represent the supplied Gram")
    det = a*c - b*b
    n = isqrt(det)
    if n*n != det or n == 0:
        raise ValueError("expected a positive square Gram determinant")
    d = gcd(*v)
    return d if 1 < d < n and n % d == 0 else None


def _all_frames(n, x):
    """Enumerate both integer triples directly; no local signs or factors are used."""
    a, b, c = _gram(n, x)
    for i in range(-n, n + 1):
        radius = isqrt(c - i*i)
        for j in range(-radius, radius + 1):
            remainder = c - i*i - j*j
            k = isqrt(remainder)
            if k*k != remainder:
                continue
            for last in ((-k, k) if k else (0,)):
                v = (i, j, last)
                # ||u+(x/n)v||=1 bounds every coordinate in a closed interval.
                intervals = [range(-((x*t+n)//n), (-x*t+n)//n+1) for t in v]
                for u in product(*intervals):
                    if sum(t*t for t in u) == a and sum(s*t for s, t in zip(u, v)) == b:
                        yield u, v


def check_small():
    checked = frames = 0
    for p, q in ((3,5), (3,7), (3,11), (5,7), (5,13), (7,11), (7,13), (11,13)):
        n = p*q
        labels = [int(jacobi_symbol(-(x*x+1), n)) for x in range(n)]
        chi_p, chi_q = (1 if r % 4 == 1 else -1 for r in (p,q))
        expected_negatives = ((p-1-chi_p)*(q-1-chi_q) - chi_p*chi_q)//2
        assert labels.count(-1) == expected_negatives
        # Cover every sign pattern occurring at the two primes, plus a seeded
        # public challenge. Factor data here grades the independent enumeration.
        examples = {}
        for x in range(n):
            signs = tuple(int(jacobi_symbol(-(x*x+1), r)) for r in (p,q))
            if all(signs):
                examples.setdefault(signs, x)
        x_public, _ = sample_mixed_gram(n, Random(9201+n))
        for x in sorted({*examples.values(), x_public}):
            gram = _gram(n, x)
            observed = Counter()
            for u, v in _all_frames(n, x):
                factor = factor_from_frame(gram, u, v)
                minors = [u[i]*v[j] - u[j]*v[i] for i,j in ((0,1),(0,2),(1,2))]
                assert gcd(*minors) == gcd(*v)
                coordinate_plane = any(i == j == 0 for i, j in zip(u, v))
                observed[n if coordinate_plane else factor or 1] += 1
            p_good, q_good = [jacobi_symbol(-gram[0], r) == 1 for r in (p,q)]
            expected = {n: 24}
            if p_good:
                expected[q] = 48
            if q_good:
                expected[p] = 48
            if p_good and q_good:
                expected[1] = 96
            assert observed == expected, (n, x, observed, expected)
            if labels[x] == -1:
                assert sum(observed.values()) == 72 and observed[n] == 24
            checked += 1
            frames += sum(observed.values())
    return checked, frames


def check_large_with_factor():
    p, q = (1 << 107) - 1, (1 << 127) - 1  # Known Mersenne primes, both 3 mod 4.
    n = p*q
    rng = Random(71942)
    successes = Counter()
    for _ in range(16):
        x, gram = sample_mixed_gram(n, rng)
        a = gram[0]
        r = p if pow(-a, (p-1)//2, p) == 1 else q
        root = pow(-a, (r+1)//4, r)
        u, w = embed_gram(a, -r*x, r*r, root)
        v = tuple((n//r)*c for c in w)
        assert factor_from_frame(gram, u, v) == n//r
        successes[(n//r).bit_length()] += 1
    return dict(successes)


if __name__ == "__main__":
    print("Small sign controls and exact marked frames:", check_small())
    example = _gram(65, 4)
    assert factor_from_frame(example, (-4,-1,0), (60,20,15)) == 5
    assert factor_from_frame(example, (4,1,0), (-65,0,0)) is None
    assert factor_from_frame(_gram(77,1), (1,1,0), (-14,-63,42)) == 7
    print("Radius 65 and 77 explicit noncoordinate frames return 5 and 7")
    print("234-bit factor-supplied checks by returned bit length:", check_large_with_factor())
