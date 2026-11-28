"""Evaluate the prescribed gap comparison given the two prime factors."""

from dataclasses import astuple, dataclass
from math import gcd
from random import Random

from check_single_quadratic_count import direct_contrast, stationary_correction


@dataclass(frozen=True)
class Query:
    n: int
    c: int
    j: int
    width: int
    offset: int
    alpha: int
    beta: int

    def weight(self, a):
        if (a+self.offset) % self.n >= self.width:
            return 0
        return 1 if (self.alpha*a+self.beta) % self.n <= self.n//2 else -1


def block_change(query, a):
    """Direct inner-lift formula on a specified unit lower digit."""
    n, c = query.n, query.c
    d = gcd(a*a-c*c, n)
    inverse = pow(a, -1, n*n)
    phase = (a+c*c*inverse) % (n*n)
    old = (phase//n) % d
    new = (old+2*c*inverse*query.j) % d
    return d*query.weight(a)*(int(new <= d//2)-int(old <= d//2))


def stratum_change(query, d, m, sign):
    """Exact contribution for gcd=d and a=sign*c mod d, in m iterations."""
    c, j = query.c, query.j
    base = (sign*c) % d
    b0 = (2*sign*c) % d
    inverse_m = pow(m, -1, d)
    inverse_d = pow(d, -1, m)
    intercept = ((2*sign*c-b0)//d) * inverse_m
    total = 0
    for index in range(m):
        a = base+d*index
        x = a % m
        if x in (0, c % m, (-c) % m):
            continue
        k = ((x+c*c*pow(x, -1, m)-b0)*inverse_d) % m
        old = (intercept-inverse_m*k) % d
        new = (old+2*sign*j) % d
        total += d*query.weight(a)*(int(new <= d//2)-int(old <= d//2))
    return total


def stationary_change(query, p, q):
    mixed = 1+p*((-2*pow(p, -1, q)) % q)
    return sum(block_change(query, (sign*query.c*mixed) % query.n)
               for sign in (-1, 1))


def estimate(query, p, q):
    """Integer estimate of corrected raw contrast, within N/4 on all inputs."""
    assert p*q == query.n and p != q and gcd(query.c, query.n) == 1
    cutoff = 2**22*(2+query.n.bit_length())**8
    answer = stationary_change(query, p, q)
    omitted = 0
    for d, m in ((p, q), (q, p)):
        if m >= cutoff:
            omitted += 2
        else:
            answer += sum(stratum_change(query, d, m, sign) for sign in (-1, 1))
    return answer, omitted


def comparison(query, p, q, polarity):
    value, _ = estimate(query, p, q)
    return 2*polarity*value >= query.n


def split_from_gap(n, oracle):
    """Dyadic-window stage after trial division, with only gap access to N."""
    lo, hi = 2, n//2
    width = 1 << ((hi-lo+1).bit_length()-1)
    calls = 0

    def ask(start, size, j, polarity):
        nonlocal calls
        calls += 1
        query = Query(n, 1, j % n, size, (-start) % n, 0, 0)
        return oracle(query, polarity)

    for i in range(n.bit_length()+1):
        for j in (2**i, -(2**i)):
            for polarity in (-1, 1):
                for start in (lo, hi-width+1):
                    if not ask(start, width, j, polarity):
                        continue
                    size = width
                    while size > 1:
                        size //= 2
                        if not ask(start, size, j, polarity):
                            start += size
                    factor = gcd(start-1, n)
                    assert 1 < factor < n
                    return factor, calls
    raise AssertionError("no changing stationary shift")


def grid_count(d, m, intercept, shift):
    inv = pow(m, -1, d)
    return sum(int((intercept+shift-inv*k) % d <= d//2)
               - int((intercept-inv*k) % d <= d//2) for k in range(m))


def check():
    rng = Random(264031)
    small = twisted = grid_cases = 0
    for p, q in ((3, 5), (5, 7), (3, 23), (11, 17), (31, 37), (7, 101)):
        n = p*q
        for _ in range(18):
            c = rng.randrange(1, n)
            while gcd(c, n) != 1:
                c = rng.randrange(1, n)
            args = [n, c, rng.randrange(n), rng.randrange(1, n+1)]
            query = Query(*args, *(rng.randrange(n) for _ in range(3)))
            exact = sum(block_change(query, a) for a in range(1, n)
                        if gcd(a, n) == 1 and a not in (c, n-c))
            predicted, omitted = estimate(query, p, q)
            assert predicted == exact and omitted == 0
            if n <= 187:
                original = astuple(query)
                assert exact == direct_contrast(*original)-stationary_correction(*original)
                small += 1
            for d, m in ((p, q), (q, p)):
                for sign in (-1, 1):
                    expected = sum(block_change(query, a) for a in range(1, n)
                                   if gcd(a, n) == 1 and gcd(a*a-c*c, n) == d
                                   and (a-sign*c) % d == 0)
                    assert expected == stratum_change(query, d, m, sign)
                    twisted += 1
                for j in range(d):
                    assert abs(grid_count(d, m, c, j)) <= 8
                    grid_cases += 1

    # Known Mersenne primes make the asymptotic omission branch reachable.
    # Narrow windows give a second certificate using only their lower digits.
    large = []
    for p, q in ((11, 2**127-1), (2**107-1, 2**127-1)):
        n = p*q
        mixed = 1+p*((-2*pow(p, -1, q)) % q)
        for i in range(20):
            query = Query(n, 1, rng.randrange(n), 17, (-mixed+i%17) % n,
                          rng.randrange(n), rng.randrange(n))
            exact = 0
            for index in range(query.width):
                a = (index-query.offset) % n
                if gcd(a, n) == 1 and a not in (1, n-1):
                    exact += block_change(query, a)
            predicted, omitted = estimate(query, p, q)
            assert omitted == (2 if p == 11 else 4)
            assert 4*abs(predicted-exact) < n
            for polarity in (-1, 1):
                signed = polarity*exact
                if 4*signed <= n or 4*signed >= 3*n:
                    assert comparison(query, p, q, polarity) == (4*signed >= 3*n)
            large.append(omitted)
    print(f"{small} literal contrasts, {twisted} exact strata and {grid_cases} grid sums passed.")
    print(f"{len(large)} narrow-window controls at 131 and 234 bits, omitting 2 or 4 strata.")
    calls = []
    for p, q in ((101, 103), (251, 257), (521, 607), (2**107-1, 2**127-1)):
        factor, query_count = split_from_gap(p*q, lambda query, sign: comparison(query, p, q, sign))
        assert factor in (p, q)
        calls.append(query_count)
    print(f"Four dyadic gap bisections returned a verified factor, using {calls} queries.")


if __name__ == "__main__":
    check()
