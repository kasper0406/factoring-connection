"""Read prime character constraints from an unfinished continued-fraction relation matrix."""

from math import gcd, isqrt, prod
from typing import NamedTuple

class CharacterRows(NamedTuple):
    primes: tuple
    rows: tuple
    divisor: int | None

    @property
    def conductor(self):
        return 8*prod(self.primes)

    @property
    def class_count(self):
        return (4*prod(p-1 for p in self.primes)) >> len(self.rows)


def character_word(value, primes):
    """The minus-sign indicators of (-1/r),(2/r), and (r/l) at each odd base prime l."""
    if value % 2 == 0:
        raise ValueError("character argument must be odd")
    word = int(value % 4 == 3) | (2 if value % 8 in (3,5) else 0)
    for i, prime in enumerate(primes, 2):
        unit = value % prime
        if not unit:
            raise ValueError("character argument must be a unit")
        if pow(unit, (prime-1)//2, prime) != 1:
            word |= 1 << i
    return word


def _signed_product(mask, basis):
    value = 1
    while mask:
        bit = mask & -mask
        value *= basis[bit.bit_length()-1]
        mask -= bit
    return value


def _extract_row(n, numerator, quotient, basis):
    """Return (root, factor mask, positive cofactor) with root²=cofactor*D(mask)."""
    mask = int(quotient < 0)
    remaining = abs(quotient)
    square = 1
    for i, signed in enumerate(basis[1:], 1):
        prime = abs(signed)
        count = 0
        while remaining % prime == 0:
            count += 1
            remaining //= prime
        square *= prime**(count//2)
        if count & 1:
            mask ^= 1 << i
            if signed < 0:
                mask ^= 1
    root = isqrt(remaining)
    if root*root == remaining:
        square *= root
        remaining = 1
    return numerator*pow(square, -1, n) % n, mask, remaining


def _combine(n, basis, x, first, y, second):
    shared = _signed_product(first & second, basis)
    return x*y*pow(shared, -1, n) % n, first ^ second


def collect_character_rows(n, primes, budget, *, stop_rank=None, stats=None):
    """Collect up to budget convergents of sqrt(n), using N alone and a chosen factor base.

    Return a proper gcd immediately, including any split from a dependent row.
    An optional stop_rank ends collection as soon as this many independent rows are available.
    """
    primes = tuple(primes)
    if (n < 3 or n % 2 == 0 or isqrt(n)**2 == n or budget < 0
            or (stop_rank is not None and stop_rank < 1)
            or tuple(sorted(set(primes))) != primes
            or any(p <= 2 or any(p % d == 0 for d in range(2, isqrt(p)+1)) for p in primes)):
        raise ValueError("expected an odd nonsquare, increasing odd base primes and valid budget")
    work = {} if stats is None else stats
    work.update(steps=0, complete=0, paired=0, dependencies=0, rank=0)
    for prime in primes:
        if n % prime == 0:
            if prime == n:
                raise ValueError("a factor-base prime equals the input")
            return CharacterRows(primes, (), prime)
    basis = (-1, 2) + tuple(p if p % 4 == 1 else -p for p in primes)
    pivots, partials = {}, {}
    root_n = isqrt(n)
    offset, denominator, term = 0, 1, root_n
    prev2, prev = 0, 1
    divisor = None
    for i in range(budget):
        work["steps"] += 1
        numerator = (term*prev+prev2) % n
        prev2, prev = prev, numerator
        offset = denominator*term-offset
        denominator = (n-offset*offset)//denominator
        term = (root_n+offset)//denominator
        signed = denominator if i % 2 else -denominator
        d = gcd(n, numerator)
        if d > 1:
            if d < n:
                divisor = d
                break
            continue
        x, mask, cofactor = _extract_row(n, numerator, signed, basis)
        if cofactor != 1:
            existing = partials.get(cofactor)
            if existing is None:
                partials[cofactor] = x, mask
                continue
            y, old_mask = existing
            x, mask = _combine(n, basis, x, mask, y, old_mask)
            x = x*pow(cofactor, -1, n) % n
            work["paired"] += 1
        else:
            work["complete"] += 1
        while mask:
            pivot = mask.bit_length()-1
            existing = pivots.get(pivot)
            if existing is None:
                pivots[pivot] = mask, x
                break
            old_mask, y = existing
            x, mask = _combine(n, basis, x, mask, y, old_mask)
        if not mask:
            work["dependencies"] += 1
            d = gcd(x-1, n)
            if 1 < d < n:
                divisor = d
                break
        if len(pivots) == stop_rank:
            break
    work["rank"] = len(pivots)
    return CharacterRows(primes, tuple(pivots[i] for i in sorted(pivots)), divisor)


def admissible_residues(result):
    """Generate the CRT classes satisfying the character rows, in output-polynomial time."""
    primes, rows, _ = result
    if result.divisor:
        raise ValueError("collection already split the input")
    pivots = {mask.bit_length()-1: mask for mask, _ in rows}
    dimension = 2+len(primes)
    free = [i for i in range(dimension) if i not in pivots]
    bins = []
    modulus = 8
    for p in primes:
        square = {k*k % p for k in range(1, (p+1)//2)}
        bins.append((p, pow(modulus, -1, p), [k for k in range(1,p) if k in square],
                     [k for k in range(1,p) if k not in square]))
        modulus *= p
    # For each free sign vector, solve the triangular character constraints from low pivot to high.
    units8 = (1, 7, 5, 3)
    for index in range(1 << len(free)):
        word = sum(((index >> j) & 1) << pos for j, pos in enumerate(free))
        for i in sorted(pivots):
            if (word & pivots[i]).bit_count() % 2:
                word |= 1 << i
        residues = [units8[word & 3]]
        modulus = 8
        for i, (p, inverse, plus, minus) in enumerate(bins, 2):
            choices = minus if word & (1 << i) else plus
            residues = [r+modulus*((s-r)*inverse % p) for r in residues for s in choices]
            modulus *= p
        yield from residues
