"""HARNESS / DIAGNOSTIC ONLY: ground truth for sections of F_N using the factors.

Never import this from operational code (src/trackA/selectors.py etc.).

For squarefree N = prod of distinct primes l, by Lucas (notes/seed_lemmas.md L7):
    supp F_N = U_l { l*j : 0 < l*j < N, j digit-dominated by N/l in base l }.
`section_nonzero(N, factors, r, M)` decides whether S_{r,M}(F_N) != 0, i.e.
whether some support index is = r (mod M), without enumerating the support:
for each l, it enumerates the digit-dominated j only modulo the needed
congruence (feasible for the small-digit-count cases we use: balanced and
3-prime N; falls back to explicit enumeration otherwise).
"""
from __future__ import annotations

from itertools import product
from math import comb, gcd
from typing import Dict, Iterable, List, Sequence


def base_digits(m: int, b: int) -> List[int]:
    out = []
    while m:
        out.append(m % b)
        m //= b
    return out or [0]


def dominated_set_size(m: int, b: int) -> int:
    s = 1
    for d in base_digits(m, b):
        s *= d + 1
    return s


def support_explicit(N: int, factors: Sequence[int]) -> List[int]:
    """Explicit interior support of F_N for squarefree N (small cases)."""
    S = set()
    for l in set(factors):
        m = N // l
        digs = base_digits(m, l)
        for js in product(*[range(d + 1) for d in digs]):
            j = sum(d * l ** i for i, d in enumerate(js))
            if 0 < j < m:
                S.add(l * j)
    return sorted(S)


def _class_hit_for_prime(N: int, l: int, r: int, M: int) -> List[int]:
    """Indices l*j (0<j<N/l, j dominated by N/l in base l) with l*j = r mod M.

    Strategy: write j = j_hi * l + j_0 (j_0 = lowest base-l digit).  For each
    admissible high part j_hi (digit-dominated by m // l), the constraint on
    j_0 is a single congruence mod M/gcd(l,M) (l odd, M power of 2 -> l
    invertible mod M), solved directly.  Cost: number of admissible j_hi.
    """
    m = N // l
    digs = base_digits(m, l)
    d0, hi_digits = digs[0], digs[1:]
    if M > 1 and gcd(l, M) != 1:  # fallback: explicit enumeration
        return [i for i in support_explicit(N, [l]) if i % M == r % M]
    linv = pow(l, -1, M) if M > 1 else 0
    target = (r * linv) % M if M > 1 else 0          # need j = target (mod M)
    hits = []
    his = [0] if not hi_digits else None
    if his is None:
        his = []
        for js in product(*[range(d + 1) for d in hi_digits]):
            his.append(sum(d * l ** i for i, d in enumerate(js)))
    for jh in his:
        base = jh * l
        # j0 in [0, d0], j0 = target - base (mod M)
        j0 = (target - base) % M if M > 1 else 0
        while j0 <= d0:
            j = base + j0
            if 0 < j < m:
                hits.append(l * j)
            if M == 1:
                j0 += 1
            else:
                j0 += M
    return hits


def section_support(N: int, factors: Sequence[int], r: int, M: int) -> List[int]:
    """Support indices of F_N in the class r mod M (squarefree N)."""
    assert len(set(factors)) == len(factors), "squarefree only"
    S = set()
    for l in factors:
        S.update(_class_hit_for_prime(N, l, r % M, M))
    return sorted(S)


def section_nonzero(N: int, factors: Sequence[int], r: int, M: int) -> bool:
    return bool(section_support(N, factors, r, M))


def section_coeffs_explicit(N: int, r: int, M: int) -> Dict[int, int]:
    """{j: a_{r+jM}} nonzero coefficients of S_{r,M}(Y) for F_N (explicit, small N)."""
    out = {}
    i = r
    j = 0
    while i <= N:
        if 0 < i < N:
            a = comb(N, i) % N
            if a:
                out[j] = a
        i += M
        j += 1
    return out


def count_hit_classes(N: int, factors: Sequence[int], M: int) -> int:
    """Number of residues r mod M with S_{r,M} != 0 (explicit support; small N)."""
    sup = support_explicit(N, factors)
    return len({i % M for i in sup})
