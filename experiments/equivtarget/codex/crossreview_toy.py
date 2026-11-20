"""Exact-integer checks of the two predictions in crossreview-prereg.txt."""
from itertools import combinations
from math import gcd, isqrt


def prime(n):
    return n >= 2 and all(n % i for i in range(2, isqrt(n) + 1))


def mul(z, w):
    a, b = z
    c, d = w
    return a*c - b*d, a*d + b*c - b*d


def conj(z):
    a, b = z
    return a - b, -b


def norm(z):
    return mul(z, conj(z))[0]


def representations(n):
    bound = isqrt(4*n // 3)
    for b in range(-bound, bound + 1):
        disc = 4*n - 3*b*b
        r = isqrt(disc)
        if r*r == disc and (r + b) % 2 == 0:
            for a in sorted({(b-r)//2, (b+r)//2}):
                yield a, b


def mod(z, n):
    return tuple(v % n for v in z)


def kani(mu, other_mu, exclude_conjugate):
    m = norm(mu)
    assert norm(other_mu) == m and mod(mu, 3) == mod(other_mu, 3)
    primes = [p for p in range(7, 98) if p % 3 == 1 and m % p and prime(p)]
    for p, q in combinations(primes, 2):
        n = p*q
        for e in range((m*n).bit_length(), (m*n).bit_length() + 5):
            h = 2**e
            gap = h - m*n
            if e % 2 == 0 or gap % 3 != 1 or not prime(gap):
                continue
            ell = next(c*4*h-1 for c in range(3, 9999, 3) if prime(c*4*h-1))
            alpha = next(representations(gap))
            reps = list(representations(n))
            nu = next(z for z in reps if min(z) > 0)
            genuine, eta = mul(mu, nu), mul(other_mu, nu)
            divs = [mu, conj(mu)] if exclude_conjugate else [mu]
            assert all(mod(eta, h) != mod(mul(d, z), h) for d in divs for z in reps)
            assert all(mod(mul(conj(d), eta), m) != (0, 0) for d in divs)
            assert mod(eta, 3) == mod(genuine, 3)
            assert norm(eta) == m*n and norm(alpha) + norm(eta) == h
            # The off-diagonal of F^*F vanishes since its two products commute.
            assert mul(conj(alpha), conj(eta)) == mul(conj(eta), conj(alpha))
            print(f'Kani {m=} {p=} {q=} {n=} {e=} {ell=} {alpha=} {mu=} {other_mu=} '
                  f'{nu=} {genuine=} {eta=} all {len(reps)} norm-n representations checked')
            return
    raise AssertionError('No promised Kani example')


def acd():
    d = 2**14
    p = next(x for x in range(2**31+1, 2**32, d) if prime(x))
    q = next(x for x in range(2**32-d+1, 2**31, -d) if prime(x))
    n = p*q
    samples = [(n-1)//d+i for i in range(16)]
    weight = 2**32//d + 15
    for i, x in enumerate(samples):
        assert x % p == (p-1)//d+i and x % q == (q-1)//d+i
        assert x % p != x % q and max(x % p, x % q) <= weight

    def length_squared(t):
        def center(a):
            return (a + n//2) % n - n//2
        return (t*weight)**2 + sum(center(t*x)**2 for x in samples)

    assert gcd(d, n) == 1 and length_squared(d) < min(length_squared(p), length_squared(q))
    print(f'ACD {d=} {p=} {q=} {n=} {weight=} '
          f'first_residues={(samples[0]%p, samples[0]%q)} '
          f'floor_norms_d_p_q={[isqrt(length_squared(t)) for t in (d,p,q)]}')


if __name__ == '__main__':
    acd()
    kani((5, 3), conj((5, 3)), False)
    kani(mul((5, 3), (3, 1)), mul(conj((5, 3)), (3, 1)), True)
