"""Compare direct full-period L-values, noisy local controls and sine-phase recovery."""

from math import gcd, isqrt
from random import Random
from mpmath import mp
from sympy import kronecker_symbol, nextprime

if __package__:
    from .induced_euler_phase import factor_from_phase
else:
    from induced_euler_phase import factor_from_phase


def period_ratio(n, delta, t):
    """Hurwitz formula through the actual public period; no input factors."""
    s = mp.j*t
    f = abs(delta)
    table = [int(kronecker_symbol(delta, a)) for a in range(f)]
    def value(period, coprime):
        return (mp.fsum(table[a % f]*mp.zeta(s, mp.mpf(a)/period)
                        for a in range(1, period+1) if table[a % f] and gcd(a, coprime) == 1)
                / mp.power(period, s))
    return value(f*n, n)/value(f, 1)


def local_ratio(r, d, t):
    return -mp.expm1(-mp.j*t*mp.log(r))*(1+mp.expj(-t*mp.log(d)))


def verify():
    with mp.workprec(100):
        for n, delta, r in ((15, -4, 5), (35, -4, 5)):
            work = {}
            def ratio(t):
                result = period_ratio(n, delta, t)
                assert abs(result-local_ratio(r, n//r, t)) < mp.mpf("1e-18")
                return result
            assert factor_from_phase(n, ratio, stats=work, round_only=True) == r
            print(f"Input-only period reader at N={n}: {work}")

    rng = Random(271828)
    count = 0
    for bits in (16, 32, 48, 64, 96, 128):
        p = int(nextprime(2**(bits-1)+rng.randrange(2**(bits-3))))
        q = int(nextprime(3*2**(bits-2)+rng.randrange(2**(bits-3))))
        for r, d in ((p,q), (q,p)):
            for direction in (-1, 1):
                delta = next(-m for m in range(3, 512) if m % 4 == 3 and
                             int(kronecker_symbol(-m, r)) == 1 and
                             int(kronecker_symbol(-m, d)) == -1 and
                             all(m % (v*v) for v in range(2, isqrt(m)+1)))
                assert kronecker_symbol(delta, r*d) == -1
                def ratio(t):
                    # Error magnitude exactly 1/256 in the observed sine direction.
                    return local_ratio(r,d,t) + direction*mp.j*mp.expj(-t*mp.log(r*d)/2)/256
                work = {}
                assert factor_from_phase(r*d, ratio, stats=work, round_only=True) == r
                if bits <= 32:
                    assert factor_from_phase(r*d, ratio, stats=work) == r
                count += 1
        print(f"Local Euler-factor control, N bits={(p*q).bit_length()}, "
              f"largest queried t bits={work['max_height'].bit_length()}")
    p, q = 1000000007, 1000000009
    for r, d in ((p, q), (q, p)):
        assert factor_from_phase(p*q, lambda t: local_ratio(r, d, t), round_only=True) == r
        assert factor_from_phase(p*q, lambda t: local_ratio(r, d, t)) in (p, q)
    print(f"{count} noisy local controls and both close-prime orientations passed")


if __name__ == "__main__":
    verify()
