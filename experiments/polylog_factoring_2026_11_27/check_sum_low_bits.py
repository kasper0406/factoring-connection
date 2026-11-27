"""Decode low bits of the factor sum into a few larger-factor progressions."""

from math import isqrt
from statistics import median

from check_lehman_window_baseline import cohorts
from check_sum_residue_lattice import query_classes


def factor_classes(n, residue, bits):
    if bits < 1 or n % 2 == 0:
        raise ValueError('Expected positive precision and odd modulus')
    modulus = 1 << bits
    if residue % 2:
        return modulus, []
    half = (residue % modulus)//2
    delta = (half*half-n) % modulus
    if delta == 0:
        stride = 1 << ((bits+1)//2)
        return stride, [half % stride]
    valuation = (delta & -delta).bit_length()-1
    if valuation % 2:
        return modulus, []
    shift = valuation//2
    unit = delta >> valuation
    roots = [1]
    for j in range(1, bits-valuation):
        mask = (2 << j)-1
        roots = [v for r in roots for v in (r,r+(1 << j)) if (v*v-unit) & mask == 0]
    stride = 1 << (bits-shift)
    return stride, [(half+(r << shift)) % stride for r in roots]


if __name__ == '__main__':
    for bits in range(1, 9):
        modulus = 1 << bits
        for n in (15, 21, 77, 143, 257):
            for s in range(0, modulus, 2):
                step, classes = factor_classes(n,s,bits)
                assert {x for c in classes for x in range(c,modulus,step)} == {
                    x for x in range(modulus) if (x*x-s*x+n)%modulus == 0}
    for size, planted, cases in cohorts():
        if size != 48 or planted:
            continue
        thresholds = []
        for p,q in cases:
            n = p*q
            low, high = isqrt(4*n-1)+1, isqrt(9*n//2)
            for bits in range(1, (p+q).bit_length()+1):
                modulus = 1 << bits
                residue = (p+q)%modulus
                stride, classes = factor_classes(n,residue,bits)
                try:
                    f = query_classes(n,low,high,stride,classes,modulus,residue)
                except ValueError as exc:
                    if str(exc) != 'Progression exceeds the certified lattice radius':
                        raise
                    continue
                assert f in (p,q)
                thresholds.append(bits)
                break
        assert len(thresholds)==len(cases)
        print('inputs, min/median/max supplied bits:',len(cases),
              min(thresholds),median(thresholds),max(thresholds))
