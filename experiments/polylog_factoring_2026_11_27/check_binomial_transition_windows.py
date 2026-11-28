"""Check consecutive-step products, shifts, and all opposite states on small semiprimes."""

from itertools import combinations
from math import comb, gcd, isqrt, prod

from sympy import primerange

if __package__:
    from .adaptive_binomial_lift import _fold
    from .binomial_transition_windows import opposite_window, transition_gcds
else:
    from adaptive_binomial_lift import _fold
    from binomial_transition_windows import opposite_window, transition_gcds


def check_transitions():
    count = plateaus = true_center = true_edge = wrong = 0
    for p, q in combinations(primerange(11, 110), 2):
        if q >= 2*p:
            continue
        n = p*q
        lower = isqrt((n-1)//2)+1
        for d in (8, 16, 32, 64, 128):
            half = d//2
            for h in range(1, d//4, 2):
                t = h*n//d
                length = (isqrt(2*n)+1)//d+2
                numerator = prod(2*(t+j)-1 for j in range(1, length+1))
                denominator = prod(t+j for j in range(1, length+1))
                assert transition_gcds(n, h, d, length) == (gcd(n, numerator), gcd(n, denominator))
                count += 1
                bound = (lower-half)//d
                if bound > 0:
                    base = comb(2*t, t) % n == 0
                    for offset in range(-bound, bound+1):
                        u = t+offset
                        assert (comb(2*u, u) % n == 0) == base
                        plateaus += 1
            if d < 16 or half >= lower:
                continue
            for a in range(1, half, 2):
                if (a*a+n) % d:
                    continue
                numerator, denominator = opposite_window(n, a, d)
                if a not in (p % half, q % half):
                    wrong += 1
                    continue
                h = _fold((half-1)*pow(a, -1, d), half)
                scaled = {h*p % d, h*q % d}
                if scaled == {half-1, half+1}:
                    assert numerator in (p, q) and denominator == 1
                    true_center += 1
                else:
                    assert scaled == {1, d-1} and numerator == 1 and denominator in (p, q)
                    true_edge += 1
    return count, plateaus, true_center, true_edge, wrong


def check_large_indices():
    total = 0
    for power in range(4, 17):
        d = 1 << power
        aligned = [h for h in range(1, d//4, 2) if (-3*h) % d < d//2]
        assert aligned and aligned == [h for h in range(1, d//4, 2) if 6*h > d]
        total += len(aligned)
    return total


if __name__ == "__main__":
    print("Interval gcd comparisons, plateau points, center states, edge states, false states:", check_transitions())
    print("Aligned large numerators across thirteen dyadic periods:", check_large_indices())
