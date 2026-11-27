"""Redei-matrix 4-ranks and PARI 2-Sylow ground truth (diagnostic: uses p, q).

fund_disc(m)            fundamental discriminant of Q(sqrt(-m)), m > 0 squarefree
prime_disc_factors(D)   list of (prime, prime discriminant) with product D
r4_redei(D)             4-rank via Redei matrix (needs the factorisation of D)
two_sylow(D)            PARI quadclassunit 2-Sylow invariants (list of 2-powers)
"""
import functools
import cypari2
from sympy import factorint, jacobi_symbol

pari = cypari2.Pari()
pari.allocatemem(2 * 10**9, silent=True)


def kron(a, n):
    """Kronecker symbol (a/n), n prime (incl. 2)."""
    return int(pari.kronecker(a, n))


def fund_disc(m):
    """Fundamental discriminant of Q(sqrt(-m)), m squarefree positive."""
    d = -m
    return d if d % 4 == 1 else 4 * d


def prime_disc_factors(D, fac=None):
    fac = fac or factorint(abs(D))
    out = []
    rest = D
    for ell in sorted(fac):
        if ell == 2:
            continue
        ds = ell if ell % 4 == 1 else -ell
        out.append((ell, ds))
        rest //= ds
    if rest != 1:
        assert rest in (-4, 8, -8), (D, rest)
        out.append((2, rest))
    return out


def rank_f2(rows):
    rows = [r for r in rows]
    rk = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = None
        for i in range(rk, len(rows)):
            if rows[i][c]:
                piv = i
                break
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                rows[i] = [x ^ y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def redei_matrix(D, fac=None):
    pf = prime_disc_factors(D, fac)
    t = len(pf)
    R = [[0] * t for _ in range(t)]
    for i, (pi, di) in enumerate(pf):
        for j, (pj, dj) in enumerate(pf):
            if i != j:
                R[i][j] = 0 if kron(dj, pi) == 1 else 1
        R[i][i] = sum(R[i][j] for j in range(t) if j != i) % 2
    return R, pf


def r4_redei(D, fac=None):
    R, pf = redei_matrix(D, fac)
    return len(pf) - 1 - rank_f2(R)


@functools.lru_cache(maxsize=None)
def two_sylow(D):
    cyc = pari.quadclassunit(D)[1]
    out = []
    for c in cyc:
        c = int(c)
        e = 1
        while c % 2 == 0:
            c //= 2
            e *= 2
        if e > 1:
            out.append(e)
    return sorted(out, reverse=True)


def rank_2k(syl, k):
    """2^k-rank: number of cyclic factors divisible by 2^k."""
    return sum(1 for e in syl if e % (2**k) == 0)
