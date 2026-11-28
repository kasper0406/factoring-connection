"""Invert a fixed-field packed query with known split conductor prime(s)."""

from itertools import permutations, product
from math import gcd, isqrt, prod
from pathlib import Path
from random import SystemRandom
import sys
from typing import NamedTuple

from sympy import factorint, isprime

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.thmA.arith import sqrt_mod_p

if __package__:
    from .fixed_field_packing import Parameters, packing_params
    from .fixed_field_walk import conjugate_vector
    from .four_square_sampler import _multiply
    from .gram_embedding import embed_gram
else:
    from fixed_field_packing import Parameters, packing_params
    from fixed_field_walk import conjugate_vector
    from four_square_sampler import _multiply
    from gram_embedding import embed_gram


def _norm(v):
    return sum(a*a for a in v)


def _dot(v, w):
    return sum(a*b for a, b in zip(v, w))


def _gaussian_prime(prime, rng):
    """Cornacchia at a supplied split prime. Each random root trial succeeds with probability 1/2."""
    while True:
        r = pow(rng.randrange(1, prime), (prime-1)//4, prime)
        if r*r % prime == prime-1:
            break
    a, b = prime, r
    while b*b > prime:
        a, b = b, a % b
    c = isqrt(prime-b*b)
    if c*c+b*b != prime:
        raise ArithmeticError("invalid prime sum-of-squares certificate")
    return b, c


def gaussian_decompositions(primes, rng=None):
    """Every ordered positive (B,C) of norm T for one or two known distinct split primes."""
    primes = tuple(map(int, primes))
    if (len(primes) not in (1, 2) or len(set(primes)) != len(primes)
            or any(p % 4 != 1 or not isprime(p) for p in primes)):
        raise ValueError("expected one or two distinct known primes congruent to 1 modulo 4")
    if rng is None:
        rng = SystemRandom()
    a, b = _gaussian_prime(primes[0], rng)
    pairs = {(a, b)}
    if len(primes) == 2:
        c, d = _gaussian_prime(primes[1], rng)
        pairs = {(abs(a*c-sign*b*d), abs(a*d+sign*b*c)) for sign in (-1, 1)}
    return tuple(sorted(pairs | {(imag, real) for real, imag in pairs}))


class Packing(NamedTuple):
    n: int
    x: int
    t: int
    ell: int
    params: Parameters


class Mark(NamedTuple):
    packing: Packing
    u: tuple
    z: tuple


def parameter_candidates(norm, modulus, conductor_primes, rng=None):
    """All arithmetically compatible packings; the original semiprime promise is not used to pick one."""
    primes = tuple(conductor_primes)
    conductor = prod(primes)
    if norm <= 0 or modulus <= 0 or conductor < 2 or norm % conductor**2:
        raise ValueError("invalid norm, modulus, or stated conductor")
    a = norm//conductor**2
    x = isqrt(a-1)
    if x*x+1 != a:
        raise ValueError("base norm is not x²+1")
    candidates = []
    for b, c in gaussian_decompositions(primes, rng):
        scale = 2*a*b
        if modulus % scale:
            continue
        quotient = modulus//scale
        t = gcd(c, quotient)
        n, ell = c//t, quotient//t
        if ell <= 1 or gcd(ell, 210*a*n*scale*t*conductor) != 1 or not isprime(ell):
            continue
        try:
            params = packing_params(n, x, b, t)
        except ValueError:
            continue
        candidates.append(Packing(n, x, t, ell, params))
    return tuple(candidates)


def unpack_candidates(norm, modulus, residue, conductor_primes, rng=None):
    """Recover every candidate with a valid first vector and fine Gram mark."""
    if len(residue) != 3:
        raise ValueError("expected three residue coordinates")
    output = []
    for candidate in parameter_candidates(norm, modulus, conductor_primes, rng):
        n, x, _, ell, (_, p, q, _, _) = candidate
        inv = pow(p, -1, q)
        u = tuple((inv*c+q//2) % q-q//2 for c in residue)
        if _norm(u) != x*x+1:
            continue
        z = tuple(((p*i-c)//q) % ell for i, c in zip(u, residue))
        if (_norm(z)-n*n) % ell == 0 and (_dot(u, z)+n*x) % ell == 0:
            output.append(Mark(candidate, u, z))
    return tuple(output)


def factored_partners(n, x, u, z, fine, divisor, *, rng=None):
    """Return all marked partners given a proper factor of the promised mixed semiprime."""
    if not 1 < divisor < n or n % divisor:
        raise ValueError("expected a proper divisor of the input")
    a = x*x + 1
    split = [r for r in (divisor, n // divisor) if pow(-a, (r-1)//2, r) == 1]
    if len(split) != 1:
        return ()
    r = split[0]
    root = sqrt_mod_p(-a, r, rng if rng is not None else SystemRandom())
    if root is None:
        return ()
    u1, v1 = embed_gram(a, -r*x, r*r, int(root))
    scaled = tuple((n//r)*v for v in v1)
    found = set()
    for first, second in (((x, 1, 0), (-n, 0, 0)), (u1, scaled)):
        for order in permutations(range(3)):
            for signs in product((-1, 1), repeat=3):
                if tuple(s*first[j] for j, s in zip(order, signs)) != tuple(u):
                    continue
                v = tuple(s*second[j] for j, s in zip(order, signs))
                if all((i-j) % fine == 0 for i, j in zip(v, z)):
                    found.add(v)
    return tuple(sorted(found))


def scalar_points(norm, modulus, residue, conductor_primes, *, factor=factorint, rng=None):
    """Exact primitive points for a promised packed scalar query; default baseline factors candidate N.

    The factor callback returns a complete prime-to-exponent mapping. A genuine
    candidate must be a product of two distinct odd primes. Unioning independently
    verified witnesses removes any dependence on uniqueness of the recovered N.
    """
    if modulus <= 0:
        raise ValueError("the residue modulus must be positive")
    c = tuple(t % modulus for t in residue)
    found = set()
    for packing, u, z in unpack_candidates(norm, modulus, c, conductor_primes, rng):
        n, x, _, ell, (_, p, q, _, _) = packing
        factors = factor(n)
        if len(factors) != 2 or any(e != 1 for e in factors.values()) or prod(factors) != n:
            continue
        for v in factored_partners(n, x, u, z, ell, int(next(iter(factors))), rng=rng):
            w = tuple(p*i-q*j for i, j in zip(u, v))
            if _norm(w) == norm and gcd(*w) == 1 and tuple(t % modulus for t in w) == c:
                found.add(w)
    return tuple(sorted(found))


def collapse(mark, h, base):
    """For a norm-T spinor realizing this mark, return its integral two quaternion components."""
    packing, u, z = mark
    n, x, t, ell, (b, p, q, conductor, _) = packing
    a = x*x+1
    if len(h) != 4 or len(base) != 3 or _norm(h) != conductor or _norm(base) != a:
        raise ValueError("incorrect quaternion or base norm")
    w = conjugate_vector(h, base)
    offsets = tuple(p*i-j for i, j in zip(u, w))
    if any(d % q for d in offsets):
        raise ValueError("rotation misses the coarse mark")
    v = tuple(d//q for d in offsets)
    if _norm(v) != n*n or _dot(u,v) != -n*x or tuple(j % ell for j in v) != z:
        raise ValueError("rotation misses the fine Gram mark")
    left = _multiply((0, *u), h)
    right = _multiply(h, (0, *base))
    plus = tuple(i+j for i, j in zip(left, right))
    minus = tuple(i-j for i, j in zip(left, right))
    if any(i % b for i in plus) or any(i % t for i in minus):
        raise ArithmeticError("quaternion components failed integral divisibility")
    return tuple(i//b for i in plus), tuple(i//t for i in minus)


def restore(mark, base, plus, minus):
    """Return the integral quaternion realizing these components, or None if they fail the mark."""
    packing, u, _ = mark
    n, x, t, _, params = packing
    a = x*x+1
    if len(plus) != 4 or len(minus) != 4 or _norm(plus) != 4*a or _norm(minus) != 4*a*n*n:
        return None
    total = tuple(params.B*i+t*j for i, j in zip(plus, minus))
    product = _multiply((0, *u), total)
    if any(d % (2*a) for d in product):
        return None
    h = tuple(-d//(2*a) for d in product)
    try:
        if collapse(mark, h, base) == (tuple(plus), tuple(minus)):
            return h
    except ValueError:
        pass
    return None
