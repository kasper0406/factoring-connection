"""Numerical controls for the conditioning of quartic power means."""

from math import comb, factorial
from mpmath import mp
from sympy import nextprime


def profile_matrix(size):
    return mp.matrix([
        [mp.log1p(1 / (2*j - mp.mpf(2)/(5*e))) for e in range(1, size+1)]
        for j in range(1, size+1)
    ])


def inverse_bound(size):
    return size * (2*size+1)**(2*size-1) * (3*size**2)**(size-1)


def paired_profiles(order):
    base = order * (order - 1).bit_length()
    profiles = [[], []]
    for offset in range(order + 1):
        profiles[offset % 2].extend([base + offset] * comb(order, offset))
    # One common simple factor keeps both inputs from being perfect powers.
    for profile in profiles:
        profile.append(1)
    return profiles


def log_mean(primes, profile, power):
    return sum(mp.log(5*e*power - 2 - mp.mpf(5*e*power + (1 if p%4 == 1 else -1))/p)
               for p, e in zip(primes, profile))


if __name__ == "__main__":
    with mp.workdps(100):
        for size in range(1, 9):
            matrix = profile_matrix(size)
            inverse = matrix**-1
            row_norm = max(sum(abs(inverse[i, j]) for j in range(size)) for i in range(size))
            assert row_norm < inverse_bound(size)
            hist = mp.matrix([3*i+1 for i in range(size)])
            error = mp.matrix([(-1)**i / mp.mpf(8*inverse_bound(size)) for i in range(size)])
            recovered = mp.lu_solve(matrix, matrix*hist + error)
            assert all(mp.nint(x) == y for x, y in zip(recovered, hist))

        for order in (4, 6, 8):
            left, right = paired_profiles(order)
            prime = 2**(2**order)
            primes = []
            for _ in left:
                prime = int(nextprime(prime))
                primes.append(prime)
            base = min(left[:-1] + right[:-1])
            bound = mp.mpf(factorial(order-1))/(base-mp.mpf(2)/5)**order + 8*len(primes)/mp.mpf(primes[0])
            gaps = [abs(log_mean(primes, left, j) - log_mean(primes, right, j))
                    for j in (1, 2, 19, 10**6)]
            assert max(gaps) < bound
            print(order, max(left + right), mp.nstr(max(gaps), 6), mp.nstr(bound, 6))
