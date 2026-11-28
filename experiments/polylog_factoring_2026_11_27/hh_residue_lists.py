"""Split a balanced odd semiprime from a list containing its smaller prime modulo m."""

from math import gcd, isqrt

if __package__:
    from .hh_low_residue import (
        _affine_class, hh_collision, hh_giants, hh_order_base, hh_parameters,
    )
else:
    from hh_low_residue import (
        _affine_class, hh_collision, hh_giants, hh_order_base, hh_parameters,
    )
from src.thmA.arith import iroot


def _recover_common_order(n, m, base, powered_order, modulus_primes, stats):
    # At the first full repetition, both local orders of base^(m^2) equal i.
    # Reconstruct each component removed by the stride; the other prime-power
    # components of both local orders already agree with those in i.
    annihilator = powered_order*m*m
    order = powered_order
    for prime in modulus_primes:
        while order % prime == 0:
            order //= prime
        exponent = annihilator
        while exponent % prime == 0:
            exponent //= prime
        value = pow(base, exponent, n)
        while value != 1:
            divisor = gcd(value-1, n)
            if 1 < divisor < n:
                stats["status"] = "component_order"
                return divisor
            order *= prime
            value = pow(value, prime, n)
    divisor = gcd(order, n)
    if 1 < divisor < n:
        stats["status"] = "order_gcd"
        return divisor
    if divisor != 1 or order <= 8*isqrt(isqrt(n)):
        raise ValueError("common-order recovery needs the promised large-order base")
    lattice_work = {}
    divisor = _affine_class(n, order, 1, lattice_work)
    if divisor is None:
        raise ArithmeticError("missing prime in the common-order class")
    stats.update(status="common_order", common_order=order, affine_intervals=lattice_work["intervals"])
    return divisor


def _search_with_base(n, m, residues, buckets, lam, base, stats, modulus_primes=(2,)):
    # Undo the single-class exponent normalization: every class then uses
    # exactly the same baby table with stride m^2.
    def giants():
        for residue in residues:
            g = gcd(residue-1, m)
            for a, b, j, exponent in hh_giants(n, m, residue, buckets, lam):
                yield a, b, j, exponent*g

    divisor = hh_collision(n, m, base, 1, lam, giants(), stats=stats)
    if divisor is None and stats["status"] == "short_order":
        return _recover_common_order(n, m, base, stats["baby_count"], modulus_primes, stats)
    return divisor


def factor_from_residue_list(n, m, residues, *, modulus_primes=None, stats=None):
    """Return a divisor if the list contains the smaller prime modulo m.

    Promise n=p*q with distinct odd primes p<q<2p. For a non-dyadic m,
    supply its distinct prime divisors in modulus_primes. Incorrect classes
    may still produce a verified incidental divisor.
    """
    if n < 15 or n % 2 == 0 or m < 2:
        raise ValueError("expected an odd balanced semiprime and a modulus >= 2")
    if modulus_primes is None:
        if m & (m-1):
            raise ValueError("a non-dyadic modulus needs its distinct prime divisors")
        modulus_primes = (2,)
    else:
        from sympy import isprime
        modulus_primes = sorted(set(modulus_primes))
        remaining = m
        for prime in modulus_primes:
            if prime < 2 or remaining % prime or not isprime(prime):
                raise ValueError("incorrect prime divisors for the modulus")
            while remaining % prime == 0:
                remaining //= prime
        if remaining != 1:
            raise ValueError("incomplete prime divisors for the modulus")
    common = gcd(m, n)
    if common != 1:
        if common < n:
            if stats is not None:
                stats["status"] = "modulus_gcd"
            return common
        for prime in modulus_primes:
            if n % prime == 0:
                if stats is not None:
                    stats["status"] = "modulus_gcd"
                return prime
        raise ValueError("modulus factorization inconsistent with the input")
    residues = sorted({r % m for r in residues})
    if any(gcd(r, m) != 1 for r in residues):
        raise ValueError("factor residues must be units modulo m")
    work = {} if stats is None else stats
    count = len(residues)
    work.update(classes=count, status="empty", baby_count=0, giant_count=0)
    if not residues:
        return None
    if m**8 >= n*count**8:
        work.update(method="affine", affine_intervals=0)
        for residue in residues:
            lattice_work = {}
            divisor = _affine_class(n, m, residue, lattice_work)
            work["affine_intervals"] += lattice_work["intervals"]
            if divisor:
                work["status"] = "affine_factor"
                return divisor
        return None

    buckets = 4*(1+iroot(n//(m**3*count**2), 5))
    lam = hh_parameters(n, m, buckets)
    bound = 8*(isqrt(isqrt(n))+1)
    work.update(method="HH", buckets=buckets, lam=lam, order_bound=bound)
    divisor, base = hh_order_base(n, bound)
    if divisor:
        work["status"] = "order_setup"
        return divisor
    return _search_with_base(n, m, residues, buckets, lam, base, work, modulus_primes)


def _shared_prefix_list(n, bits, shared, budget):
    """Grow product-compatible residues; a pair of differing bits spends one branch."""
    modulus = 2
    states = [(1, 0)]
    tested = 0
    for _ in range(1, bits):
        next_modulus = 2*modulus
        want = shared & modulus
        following = []
        for low, spent in states:
            partner = n*pow(low, -1, next_modulus) % next_modulus
            tested += 1
            if partner & modulus:
                # The two lifts are (0,1) and (1,0) in their new bits.
                if not want and spent < budget:
                    following.extend(((low, spent+1), (low+modulus, spent+1)))
            else:
                # The two possibilities are (0,0) and (1,1).
                following.append((low+want, spent))
        states, modulus = following, next_modulus
        if not states:
            break
    return [a for a, _ in states], tested


def factor_from_shared_bits(n, bits, shared, *, complete_list=False, stats=None):
    """Factor using the supplied low bits of p & q on an odd balanced semiprime.

    Increasing the binary-branch budget keeps every run finite and eventually
    includes the true residue. complete_list=True grows the entire compatible
    tree once. The latter has the sharp density-one entropy bound.
    """
    if n < 15 or n % 2 == 0 or bits < 1 or not 0 < shared < (1 << bits) or shared % 2 == 0:
        raise ValueError("expected an odd balanced semiprime and a low AND with bit zero set")
    modulus = 1 << bits
    tested = 0
    work = {} if stats is None else stats
    last_budget = bits-shared.bit_count()
    budgets = (last_budget,) if complete_list else range(last_budget+1)
    for budget in budgets:
        residues, visits = _shared_prefix_list(n, bits, shared, budget)
        tested += visits
        if not residues:
            continue
        factor = factor_from_residue_list(n, modulus, residues, stats=work)
        if factor:
            work.update(bits_supplied=bits, branch_budget=budget, prefix_nodes_tested=tested)
            return factor
    raise ValueError("no balanced factor compatible with the supplied shared bits")
