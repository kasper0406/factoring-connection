"""Check short Chebyshev traces for a split or a two-class certificate."""

from math import gcd

if __package__:
    from .hh_residue_lists import factor_from_residue_list
else:
    from hh_residue_lists import factor_from_residue_list


def trace_shortlist(n, traces, min_level, max_level, *, stats=None):
    """At level j test C_(2^j)(x)=0 mod n, C_2(x)=x^2-2.

    A proper gcd is a divisor. A full gcd certifies that both prime divisors
    are ±1 modulo 2^(j+2) and the balanced class search can finish.
    """
    if n < 15 or n % 2 == 0 or not 1 <= min_level <= max_level:
        raise ValueError("expected an odd balanced semiprime and positive ordered levels")
    work = {} if stats is None else stats
    work.update(probes=0, status="empty")
    for x in traces:
        value = x % n
        for level in range(1, max_level+1):
            value = (value*value-2) % n
            if level < min_level:
                continue
            work["probes"] += 1
            divisor = gcd(value, n)
            if 1 < divisor < n:
                work.update(status="trace_gcd", trace=x, level=level)
                return divisor
            if divisor == n:
                m = 1 << (level+2)
                search = {}
                divisor = factor_from_residue_list(n, m, (1, m-1), stats=search)
                if divisor is None:
                    raise ArithmeticError("missing factor in certified signed classes")
                work.update(status="two_classes", trace=x, level=level, modulus=m, search=search)
                return divisor
    return None
