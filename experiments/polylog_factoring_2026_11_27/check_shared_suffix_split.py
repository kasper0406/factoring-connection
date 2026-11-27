"""Try a shared binary suffix of the two factors using only their product."""

from math import isqrt
from sympy import isprime

from check_lehman_window_baseline import cohorts, window_answer
from check_sum_low_bits import odd_square_roots
from check_sum_residue_lattice import _rigorous_ok, query_classes


def shared_suffix_split(n):
    """Assume n=pq with distinct odd primes p<q<=2p."""
    low, high = isqrt(4*n-1)+1, isqrt(9*n//2)
    lower_factor_bound = isqrt(n)+1
    # query_classes uses this same conservative larger-factor interval.
    first = max(lower_factor_bound, (low+isqrt(low*low-4*n))//2)
    last = (high+isqrt(high*high-4*n))//2
    for bits in range(1,n.bit_length()+1):
        stride = 1 << bits
        radius = max(1, ((last-first)//stride+1)//2)
        if _rigorous_ok(n,lower_factor_bound,radius,9,4):
            roots = odd_square_roots(n,bits)
            return bits, query_classes(n,low,high,stride,roots,1,0)
    raise ValueError('No certified lattice radius')


if __name__ == '__main__':
    for bits, planted, cases in cohorts():
        if bits != 48 or planted:
            continue
        hits = control_hits = overlap = 0
        precisions = []
        for p, control_q in cases:
            _, split = shared_suffix_split(p*control_q)
            if split is not None:
                assert split in (p,control_q)
                control_hits += 1
            # Preserve the ratio up to a short search in this progression.
            step = 1 << 32
            q = p + ((control_q-p)//step)*step
            while not isprime(q):
                q += step
            assert p < q <= 2*p
            n = p*q
            k, split = shared_suffix_split(n)
            assert split in (p,q)
            hits += 1
            precisions.append(k)
            overlap += window_answer(n,n//4,n//2) is not None
        print('planted, unplanted, Lehman overlap, tested bits:',hits, control_hits,
              overlap, min(precisions),max(precisions))
