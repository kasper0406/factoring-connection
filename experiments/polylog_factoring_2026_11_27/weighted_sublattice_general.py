"""Sample intermediate lattices with relative weight sqrt(|L/nZ^d|), for every n."""

from random import SystemRandom

if __package__:
    from .projector_decompositions import _remove_prime
    from .sqrt_prime_marks import sqrt_mark_probability, sqrt_prime_marks
    from .sqrt_weighted_sublattice import _quadratic_sign
    from .sublattice_sampler import lattice_from_types, type_counts
else:
    from projector_decompositions import _remove_prime
    from sqrt_prime_marks import sqrt_mark_probability, sqrt_prime_marks
    from sqrt_weighted_sublattice import _quadratic_sign
    from sublattice_sampler import lattice_from_types, type_counts


def _step(value, weight, odd, p):
    """Multiply a+b√p by an integer weight and optionally by √p."""
    a,b = value
    return (weight*p*b, weight*a) if odd else (weight*a, weight*b)


def _tables(p, e, d):
    """Weighted Birkhoff transitions and totals over t future rank columns."""
    transitions, _ = type_counts(p,e,d)
    weights = [[v*p**(r//2) for r,v in enumerate(row)] for row in transitions]
    dp = [[(1,0)]*(d+1)]
    for _ in range(e):
        totals = []
        for s in range(d+1):
            a = b = 0
            for r in range(s,d+1):
                x,y = _step(dp[-1][r], weights[s][r], r%2, p)
                a += x
                b += y
            totals.append((a,b))
        dp.append(totals)

    # tails[t] is the mass of t central columns with a central predecessor.
    center = (d+1)//2
    tails = [(1,0)]
    for _ in range(e):
        tails.append(_step(tails[-1],weights[center][center],center%2,p))
    return weights, dp, tails


def _choose(masses, p, rng):
    """Choose an index with nonnegative real masses in Z[√p] by dyadic refinement."""
    cumulative = []
    active = []
    a = b = 0
    for index, (x,y) in enumerate(masses):
        sign = _quadratic_sign(x,y,p)
        if sign < 0:
            raise ArithmeticError("negative algebraic sampling mass")
        if sign:
            a += x
            b += y
            cumulative.append((a,b))
            active.append(index)
    if not active:
        raise ArithmeticError("empty conditional type distribution")
    if len(active) == 1:
        return active[0]

    ticket, scale = 0, 1
    while True:
        ticket = 256*ticket+rng.randrange(256)
        scale *= 256
        lo, hi = 0, len(active)-1
        while lo < hi:
            mid = (lo+hi)//2
            x,y = cumulative[mid]
            if _quadratic_sign(scale*x-ticket*a, scale*y-ticket*b,p) > 0:
                hi = mid
            else:
                lo = mid+1
        x,y = cumulative[lo]
        if _quadratic_sign(scale*x-(ticket+1)*a, scale*y-(ticket+1)*b,p) >= 0:
            return active[lo]


def _conditional_masses(p, d, weights, dp, tails, s, remaining, matched):
    center = (d+1)//2
    masses = []
    for r in range(s,d+1):
        a,b = _step(dp[remaining-1][r],weights[s][r],r%2,p)
        if matched and r == center:
            x,y = _step(tails[remaining-1],weights[s][r],r%2,p)
            a,b = a-x,b-y
        masses.append((a,b))
    return masses


def _marked_profile(p, e, d, rng):
    center = (d+1)//2
    weights, dp, tails = _tables(p,e,d)
    a,b = dp[e][0]
    x,y = _step(tails[e-1],weights[0][center],center%2,p)
    h = sqrt_mark_probability(p)
    exceptional = h.denominator*(a-x), h.denominator*(b-y)
    ordinary = h.numerator*a-exceptional[0], h.numerator*b-exceptional[1]
    if _choose((ordinary,exceptional),p,rng) == 0:
        return None

    reversed_profile = []
    s, matched = 0, True
    for remaining in range(e,0,-1):
        masses = _conditional_masses(p,d,weights,dp,tails,s,remaining,matched)
        s += _choose(masses,p,rng)
        reversed_profile.append(s)
        matched = matched and s == center
    return tuple(reversed(reversed_profile))


def sample_weighted_sublattice(n, dimension, rng=None, *, backend="auto"):
    """Return the column HNF of L with probability proportional to sqrt(|L/nZ^d|).

    The input n is any positive integer and d is any positive dimension. The
    prime-marker backend determines how the independent corrections are found.
    """
    if n < 1 or dimension < 1:
        raise ValueError("expected positive modulus and dimension")
    if rng is None:
        rng = SystemRandom()
    exceptions = {}
    if n > 1:
        for p in sqrt_prime_marks(n,rng,backend=backend):
            _, e = _remove_prime(n,p)
            profile = _marked_profile(p,e,dimension,rng)
            if profile is not None:
                exceptions[p] = profile
    return lattice_from_types(n,dimension,(dimension+1)//2,exceptions,rng)
