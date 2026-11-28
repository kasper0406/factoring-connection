"""Recover a common squarefree polynomial basis, or a proper divisor of N."""

from dataclasses import dataclass
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.ptebarrier.rigidity import divmod_monic, mul, trim


class _Factor(Exception):
    pass


def squarefree_basis(n, h, fs):
    """Return (factor, []) or (None, [(monic block, multiplicities), ...]).

    Coefficients are low-degree first. All inputs are monic, including the
    constant 1. Each multiplicity vector follows the input order (h, *fs).
    """
    if n < 3 or n % 2 == 0:
        raise ValueError("expected an odd modulus of at least three")

    def reduce_poly(f):
        return trim([c % n for c in f])

    inputs = [reduce_poly(f) for f in (h, *fs)]
    if any(not f or f[-1] != 1 for f in inputs):
        raise ValueError("expected monic polynomials")
    degree = max(len(f) - 1 for f in inputs)
    for d in range(3, min(n, degree + 1), 2):
        if n % d == 0:
            return d, []

    def monic_gcd(a, b):
        while b:
            d = gcd(b[-1], n)
            if d > 1:
                raise _Factor(d)
            inverse = pow(b[-1], -1, n)
            b = reduce_poly([inverse * c for c in b])
            _, remainder = divmod_monic(a, b)
            a, b = b, reduce_poly(remainder)
        return a

    def quotient(a, b):
        q, r = divmod_monic(a, b)
        assert not reduce_poly(r)
        return reduce_poly(q)

    def layers(f):
        derivative = reduce_poly([i * f[i] for i in range(1, len(f))])
        repeated = monic_gcd(f, derivative)
        remaining = quotient(f, repeated)
        order = 0
        while remaining != [1]:
            order += 1
            shared = monic_gcd(remaining, repeated)
            block = quotient(remaining, shared)
            if block != [1]:
                yield block, order
            remaining = shared
            repeated = quotient(repeated, shared)
        if repeated != [1]:
            # Trial division left a prime modulus at most the input degree.
            assert n <= degree and all(c == 0 for i, c in enumerate(repeated) if i % n)
            for block, exponent in layers(repeated[::n]):
                yield block, n * exponent

    blocks = []
    try:
        for i, f in enumerate(inputs):
            for remaining, exponent in layers(f):
                refined = []
                for block, multiplicities in blocks:
                    shared = monic_gcd(block, remaining)
                    if shared == [1]:
                        refined.append((block, multiplicities))
                        continue
                    other = quotient(block, shared)
                    if other != [1]:
                        refined.append((other, multiplicities))
                    updated = multiplicities.copy()
                    updated[i] = exponent
                    refined.append((shared, updated))
                    remaining = quotient(remaining, shared)
                if remaining != [1]:
                    updated = [0] * len(inputs)
                    updated[i] = exponent
                    refined.append((remaining, updated))
                blocks = refined
    except _Factor as found:
        return found.args[0], []
    return None, [(tuple(block), tuple(v)) for block, v in blocks]


def squarefree_support(n, h):
    divisor, blocks = squarefree_basis(n, h, [])
    if divisor is not None:
        return divisor, None
    support = [1]
    for block, _ in blocks:
        support = [v % n for v in mul(support, block)]
    return None, tuple(support)


@dataclass(frozen=True)
class TargetClassification:
    factor: int | None = None
    scalar_exponents: tuple = ()
    norm_one_exponent: int | None = None
    universal_root: int | None = None


def classify_target(n, h, m, fs, exponents, unit=1):
    """Classify unit-shift tests (c-X)^m = unit * product(F_i(c)^e_i)."""
    if m < 1 or len(h) < 2 or len(fs) != len(exponents):
        raise ValueError("expected positive m, nonconstant H and one exponent per factor")
    divisor = gcd(unit, n)
    if divisor == n:
        raise ValueError("the target multiplier must be a unit")
    if divisor > 1:
        return TargetClassification(factor=divisor)
    divisor, blocks = squarefree_basis(n, h, fs)
    if divisor is not None:
        return TargetClassification(factor=divisor)
    weights = [(block, mult[0], sum(a*b for a, b in zip(mult[1:], exponents)))
               for block, mult in blocks]
    active = [(block, in_h, weight) for block, in_h, weight in weights if weight]
    if len(active) == 1:
        block, in_h, weight = active[0]
        if len(block) == 2 and in_h and weight == m:
            divisor = gcd(unit - 1, n)
            if 1 < divisor < n:
                return TargetClassification(factor=divisor)
            if divisor == n:
                return TargetClassification(universal_root=(-block[0]) % n)
    scalar = {m}
    for block, in_h, weight in active:
        scalar.add(abs(weight))
        if in_h and weight != m:
            scalar.add(abs(m - weight))
    return TargetClassification(scalar_exponents=tuple(sorted(scalar)), norm_one_exponent=m)


def check_image_orders():
    """Check the character bound using independent small field calculations."""
    if __package__:
        from .check_translated_batch import power
    else:
        from check_translated_batch import power
    factors = ((0, 1), (1, 1), (1, 0, 1))
    checks = 0
    # X²+1 is irreducible at all four primes; the other factors are distinct linear ones.
    for r in (7, 11, 19, 23):
        for h_index, h in enumerate(factors):
            d = len(h) - 1
            field_order = r**d - 1
            norm_exponent = field_order // (r-1)
            for m in (1, 12, 492):
                for weights in ((1,0,0), (0,-3,2**20+1), (0,0,12), (-17,9,0)):
                    own_order = field_order // gcd(field_order, m - weights[h_index]*norm_exponent)
                    image_order = own_order
                    for i, weight in enumerate(weights):
                        if i != h_index:
                            image_order = lcm(image_order, (r-1)//gcd(r-1, weight))
                    hits = 0
                    for c in range(1,r-1):
                        target = 1
                        for val, exponent in zip((c,c+1,c*c+1), weights):
                            target = target * pow(val, exponent, r) % r
                        hits += power((c,-1),m,h,r) == (target,) + (0,)*(d-1)
                    # Regular c excludes 0,-1. The full étale algebra has rank four.
                    assert (image_order*hits - (r-2))**2 <= (image_order-1)**2 * 9*r
                    checks += 1
    return checks


def check_basis():
    from sympy import Poly, Symbol

    x = Symbol("x")
    quadratic = [-2, 0, 1]
    assert squarefree_basis(77, quadratic, [mul(quadratic, quadratic)]) == (
        None, [((75, 0, 1), (1, 2))])
    assert squarefree_basis(77, [0, 1], [[7, 1]]) == (7, [])
    assert squarefree_basis(15, [0, 0, 0, 1], []) == (3, [])
    assert squarefree_basis(3, [1, 0, 0, 1], []) == (None, [((1, 1), (3,))])

    assert classify_target(77, (0,1), 5, [(0,1)], [5]).universal_root == 0
    assert classify_target(77, (0,1), 5, [(0,1)], [5], 8).factor == 7
    # One weight has hundreds of bits, while the polynomial support stays fixed.
    huge = 2**270 + 17
    test = classify_target(77, (-2,0,1), 492, [quadratic, (0,1)], [huge,-huge])
    assert set(test.scalar_exponents) == {huge, 492, huge - 492}
    assert test.norm_one_exponent == 492
    # At the inert prime 11, 492 is 12 modulo 120, the quadratic norm exponent.
    block_counts = {v: 0 for v in (1, 7, 11, 77)}
    if __package__:
        from .check_translated_batch import power
    else:
        from check_translated_batch import power
    for c in range(77):
        # For X²=2, norm(a+bX)=a²−2b².
        a, b = power((c,-1), 492, quadratic, 77)
        block_counts[gcd((a-(c*c-2))**2-2*b*b,77)] += 1
    assert block_counts == {1: 0, 7: 0, 11: 55, 77: 22}

    splits = bases = local_checks = 0
    for p, q in ((7, 11), (11, 13), (13, 17)):
        for a in range(6):
            for b in range(6):
                first, second, common = [-a, 1], [-b, 1], [1, 0, 1]
                h = mul(first, second)
                fs = [mul(mul(first, first), common), mul(second, mul(common, common))]
                divisor, blocks = squarefree_basis(p * q, h, fs)
                if divisor is not None:
                    assert divisor in (p, q)
                    splits += 1
                    continue
                bases += 1
                for r in (p, q):
                    polys = [Poly.from_list(g[::-1], x, modulus=r) for g, _ in blocks]
                    one = Poly(1, x, modulus=r)
                    assert all(g.gcd(g.diff()) == one for g in polys)
                    assert all(g.gcd(k) == one for g, k in combinations(polys, 2))
                    for i, f in enumerate([h, *fs]):
                        rebuilt = one
                        for g, (_, multiplicities) in zip(polys, blocks):
                            rebuilt *= g**multiplicities[i]
                        assert rebuilt == Poly.from_list(f[::-1], x, modulus=r)
                    local_checks += 1
    return bases, splits, local_checks


if __name__ == "__main__":
    print("bases, proper splits, prime-field reconstruction checks", check_basis())
    print("correlated-target image checks", check_image_orders())
