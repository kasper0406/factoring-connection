"""Probe AKS-style coefficients on the existing missed balanced controls."""

from math import comb, gcd

from check_lehman_window_baseline import cohorts, window_answer


def cyclic_power(n, r, a):
    # Integer multiplication gives exact convolution before cyclic folding.
    shift = (r * (n-1)**2).bit_length()
    mask = (1 << shift) - 1

    def multiply(left, right):
        packed = left * right
        coeffs = [0] * r
        for j in range(2*r-1):
            coeffs[j % r] += packed & mask
            packed >>= shift
        return sum((c % n) << (j * shift) for j, c in enumerate(coeffs))

    power = 1
    base = a % n + (1 << shift)
    for bit in bin(n)[2:]:
        power = multiply(power, power)
        if bit == '1':
            power = multiply(power, base)
    return [(power >> (j * shift)) & mask for j in range(r)]


def probe(n, shifts):
    rejected = coefficient_split = scalar_split = False
    for a in shifts:
        scalar_split |= 1 < gcd(pow(a, n-1, n)-1, n) < n
        for r in (3, 5, 7, 11):
            coeffs = cyclic_power(n, r, a)
            assert sum(coeffs) % n == pow(a+1, n, n)
            coeffs[0] -= a
            coeffs[n % r] -= 1
            rejected |= any(c % n for c in coeffs)
            coefficient_split |= any(1 < gcd(c, n) < n for c in coeffs)
    return rejected, coefficient_split, scalar_split


if __name__ == '__main__':
    for n in range(5, 32):
        for r in (3, 5, 7, 11):
            for a in (1, 2, 3):
                ref = [0] * r
                for j in range(n+1):
                    ref[j % r] += comb(n, j) * a**(n-j)
                assert cyclic_power(n, r, a) == [v % n for v in ref]
    for bits, planted, cases in cohorts():
        if planted or bits == 16:
            continue
        totals = [0] * 6
        for p, q in cases:
            n = p*q
            assert window_answer(n, n//4, n//2) is None
            fixed = probe(n, (1, 2, 3))
            window = probe(n, (n//4, 3*n//8, n//2))
            totals = [x+y for x, y in zip(totals, fixed + window)]
        print(bits, len(cases), *totals)
