"""Exact unfactored reference evaluator for the quartic collision zero test."""

from math import gcd, isqrt

from sympy import integer_nthroot, isprime, primerange

try:
    from flint import fmpz, fmpz_mod_poly_ctx
except ImportError:
    fmpz = fmpz_mod_poly_ctx = None


def certified_prime(value):
    if not isprime(value):
        return False
    # SymPy uses deterministic primality testing below 2**64.
    if value < 2**64:
        return True
    if fmpz is not None:
        return bool(fmpz(value).is_prime())
    return all(value % d for d in range(3, isqrt(value) + 1, 2))


def progression_split(part, candidates, work, batched):
    def record(name, amount=1):
        if work is not None:
            work[name] = work.get(name, 0) + amount

    def checked_gcd(value):
        record("gcds")
        return gcd(value, part)

    size = max(0, (candidates.stop - candidates.start + candidates.step - 1) // candidates.step)
    if not batched or fmpz_mod_poly_ctx is None or size < 64:
        for value in candidates:
            divisor = checked_gcd(value)
            if divisor > 1:
                return divisor
        return 1

    width = isqrt(size)
    full, tail = divmod(size, width)
    ring = fmpz_mod_poly_ctx(part)

    def window(low, high):
        if high - low == 1:
            return ring([low * candidates.step, 1])
        mid = (low + high) // 2
        return window(low, mid) * window(mid, high)

    polynomial = window(0, width)
    starts = [candidates.start + k * width * candidates.step for k in range(full)]
    values = polynomial.multipoint_evaluate(starts)
    record("linear_factors", width)
    record("block_values", full)
    for start, value in zip(starts, values):
        divisor = checked_gcd(int(value))
        if 1 < divisor < part:
            return divisor
        if divisor == part:
            for index in range(width):
                divisor = checked_gcd(start + index * candidates.step)
                if divisor > 1:
                    return divisor
            raise ArithmeticError("Nonunit product with only unit factors")
    for index in range(size - tail, size):
        divisor = checked_gcd(candidates[index])
        if divisor > 1:
            return divisor
    return 1


def small_collision_witness(n, j, ell):
    """Return a verified (prime, valuation) implying divisibility, or None."""
    for prime in primerange(3, n.bit_length() + 1):
        if n % prime:
            continue
        exponent, rest = 0, n
        while rest % prime == 0:
            exponent += 1
            rest //= prime
        character = 1 if prime % 4 == 1 else -1
        if (5 * exponent * j * (prime - 1) - 2 * prime - character) % ell == 0:
            return prime, exponent
    return None


def divides_collision_power(n, j, ell, work=None, *, batched=True):
    """Return whether ell divides H(n**j), for odd n and prime ell >= 7 coprime to n."""
    if n < 1 or n % 2 == 0 or j < 1 or ell < 7 or ell % 2 == 0 or not certified_prime(ell) or n % ell == 0:
        raise ValueError("Expected odd n, positive j and a coprime prime ell >= 7")
    witness = small_collision_witness(n, j, ell)
    if witness is not None:
        if work is not None:
            work["small_witness"] = witness
        return True
    stride = 4 * ell
    pending = [n]
    while pending:
        part = pending.pop()
        if part == 1:
            continue
        if certified_prime(part):
            exponent, rest = 0, n
            while rest % part == 0:
                exponent += 1
                rest //= part
            character = 1 if part % 4 == 1 else -1
            if (5 * exponent * j * (part - 1) - 2 * part - character) % ell == 0:
                return True
            continue

        limits = {}
        for exponent in range(1, n.bit_length()):
            slope = 5 * exponent * j % ell
            if slope == 2:
                continue
            for character in (1, -1):
                residue = (slope + character) * pow((slope - 2) % ell, -1, ell) % ell
                if residue == 0:
                    continue
                target = residue + ell * ((character - residue) * pow(ell, -1, 4) % 4)
                cap = int(integer_nthroot(n, exponent)[0])
                direct = min(isqrt(part), cap)
                limits[target] = max(limits.get(target, 0), direct)
                if cap > isqrt(part):
                    other = part * pow(target, -1, stride) % stride
                    limits[other] = isqrt(part)
        split = 1
        for residue, limit in limits.items():
            start = residue if residue > 1 else stride + 1
            split = progression_split(part, range(start, limit + 1, stride), work, batched)
            if split > 1:
                break
        if split > 1:
            pending.extend([split, part // split])
    return False
