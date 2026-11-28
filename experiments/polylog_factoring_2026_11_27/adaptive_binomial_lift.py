"""Recover a divisor of an odd distinct semiprime using individual binomial zero tests."""

from math import comb, isqrt


def exact_central_zero(n, half):
    """Direct reference. Its integer cost grows with half."""
    return comb(2*half, half) % n == 0


def _lower(odd, modulus):
    return odd % modulus < modulus//2


def _fold(h, half):
    h %= 2*half
    h = min(h, 2*half-h)
    return min(h, half-h)


def aligned_multiplier(a, b, d):
    """Return odd h<d/4 with ha,hb in the same half modulo d, or None."""
    if d < 8 or d & (d-1) or not (a & b & 1):
        raise ValueError("expected odd residues and a power-of-two modulus at least eight")
    inverse = pow(a, -1, d)
    c = b*inverse % d
    half = d//2
    if c < half:
        x = 1
    else:
        r = d-c
        s = min(r, half-r)
        if s == 1:
            return None
        x = (half+s-1)//s
        if x % 2 == 0:
            x += 1
    return _fold(x*inverse, half)


def factor_by_single_binomials(n, zero=exact_central_zero, *, balanced=False,
                              affine_finish=False, stats=None):
    """zero(n,t) returns just the bit n | binom(2t,t). Promise n=pq, odd p<q.

    With the additional guarantee q<2p, balanced=True permits shifted
    center indices to resolve the opposite residue pairs at every d>=16.
    In that mode affine_finish=True completes a sufficiently long low
    residue with the existing affine-divisor window, avoiding later queries.
    """
    if n < 15 or n % 2 == 0:
        raise ValueError("expected a product of two distinct odd primes")
    m = 4
    states = {tuple(sorted((a, n*pow(a, -1, m) % m))) for a in (1, 3)}
    work = {} if stats is None else stats
    work.update(queries=0, levels=0, peak_states=len(states), opposite=0, equal=0, shifted=0)
    if affine_finish:
        if not balanced:
            raise ValueError("affine completion requires the balanced semiprime promise")
        if __package__:
            from .affine_prime_window import affine_divisors, window_radius
        else:
            from affine_prime_window import affine_divisors, window_radius
        lower, upper = isqrt((n-1)//2)+1, isqrt(n)
        radius_limit = max(2, window_radius(n))
        # For each residue the quotient interval has diameter <= (upper-lower)//m.
        first_window_modulus = (upper-lower)//(2*radius_limit+1)+1
        work.update(completion_modulus=0, affine_calls=0, maximum_dimension=0)
    while m <= 2*n:
        # Across all live states, this also finds the true smaller prime
        # before any oracle call at a level whose previous modulus exceeds it.
        for pair in states:
            for a in pair:
                if 1 < a < n and n % a == 0:
                    return a
        if affine_finish and m >= first_window_modulus:
            work["completion_modulus"] = m
            inverse = pow(m, -1, n)
            for residue in sorted({a for pair in states for a in pair}):
                lo = (lower-residue+m-1)//m
                hi = (upper-residue)//m
                if lo > hi:
                    continue
                middle = (lo+hi)//2
                window_stats = {}
                divisors = affine_divisors(n, (residue*inverse+middle) % n,
                                           hi-middle, lower, stats=window_stats)
                work["affine_calls"] += 1
                work["maximum_dimension"] = max(work["maximum_dimension"], window_stats["dimension"])
                for divisor in divisors:
                    if divisor < n:
                        return divisor
            raise ValueError("no balanced divisor compatible with the binomial answers")
        d = 2*m
        readings = {}
        following = set()
        work["levels"] += 1
        for a, _ in sorted(states):
            b = n*pow(a, -1, d) % d
            first = tuple(sorted((a, b)))
            other = tuple(sorted((a+m, (b+m) % d)))
            if first == other:
                work["equal"] += 1
                following.add(first)
                continue
            h = aligned_multiplier(a, b, d)
            shift = 0
            if h is None:
                if (a+b) % d:
                    raise ArithmeticError("unclassified multiplier symmetry")
                work["opposite"] += 1
                if not balanced or d < 16:
                    following.update((first, other))
                    continue
                # The folded multiplier places this pair at ±1 or m±1.
                # The positive shift carries both central digits, while the
                # low digit at residue 1 stays below its prime's midpoint.
                h = _fold((m-1)*pow(a, -1, d), m)
                shift = (isqrt(2*n)+1)//d+2
                predicts_zero = min(h*a % d, h*b % d) == m-1
                work["shifted"] += 1
            else:
                predicts_zero = not _lower(h*a, d)
            t = h*n//d + shift
            if t not in readings:
                readings[t] = bool(zero(n, t))
                work["queries"] += 1
            following.add(first if readings[t] == predicts_zero else other)
        states, m = following, d
        work["peak_states"] = max(work["peak_states"], len(states))
    raise ValueError("inconsistent binomial answers or semiprime promise")
