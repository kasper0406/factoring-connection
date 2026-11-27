"""Skip impossible factor sums using clipped CRT progressions."""

from math import isqrt


def sieve_interval(n, low, high, primes):
    """Return the even sums in [low, high] passing all local square tests."""
    classes = []
    for prime in primes:
        squares = {x*x % prime for x in range(prime)}
        classes.append((prime, {s for s in range(prime) if (s*s-4*n) % prime in squares}))
    starts = [low + low % 2] if low + low % 2 <= high else []
    step = 2
    probes = 0
    for index, (prime, allowed) in enumerate(classes):
        if step * prime > high - low:
            survivors = []
            for start in starts:
                for s in range(start, high+1, step):
                    probes += 1
                    if all(s % modulus in residues for modulus, residues in classes[index:]):
                        survivors.append(s)
            return sorted(survivors), probes
        inverse = pow(step, -1, prime)
        children = []
        for start in starts:
            for residue in allowed:
                first = start + step * ((residue-start)*inverse % prime)
                if first <= high:
                    children.append(first)
        starts = children
        step *= prime
    survivors = sorted(s for start in starts for s in range(start, high+1, step))
    return survivors, len(survivors)


if __name__ == '__main__':
    n = 2989601243
    primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
    low, high = isqrt(4*n-1)+1, isqrt(9*n//2)
    mid = (low+high)//2
    for left, right in ((low, high), (low, mid), (mid+1, high)):
        sums, probes = sieve_interval(n, left, right, primes)
        reference = [s for s in range(left+left%2, right+1, 2)
                     if all(pow((s*s-4*n) % p, (p-1)//2, p) in (0, 1) for p in primes)]
        assert sums == reference
        print(left, right, probes, sums)
