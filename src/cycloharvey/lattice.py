"""Target lattices for one-sided collision searches (MODEL CODE; no factors).

Target functions
----------------
Every exponent F(p, q) in Z[p, q] restricted to pq = N is a Laurent
polynomial in p:  F(p, N/p) = sum_k a_k p^k  with a_k in Z for k >= 0 and
a_k in N^{|k|} Z for k < 0  (p^i q^j = N^{min(i,j)} p^{i-j}).  With the scale
P = sqrt(N) and t = p/P this is G(t) = sum_k b_k t^k with b_k in P^{|k|} Z:
the *target lattice* Lambda_I for a support I (a set of nonzero exponents).

Residue guesses
---------------
Guessing p == r (mod M) lets one use the *integer-valued lattice*
L_{r,M} = { g in Q-span : g(p, N/p) == const (mod 1) on the class }, i.e.
(F - c)/M^e with F(p,q) == c (mod M^e) for every p == r (mod M) (q = N p^{-1}
mod M^e).  Its index over Lambda_I is M^kappa; Lemma K (notes) says
kappa = K(K+1)/2 for M coprime to small factorials, and HH Prop. 3.3 is the
case I = {-1, 1}, kappa = 3.  `class_lattice_index` computes kappa exactly.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd, prod
from typing import Dict, List, Sequence, Tuple

import flint

ALL_EXPONENTS = (-3, -2, -1, 1, 2, 3)


# --------------------------------------------------------------------------
# volume heuristic (Section 3 of the results note)
# --------------------------------------------------------------------------
def support_stats(I: Sequence[int]) -> Dict:
    """K = |I|, S = sum |i|; heuristic exponent S/(K(K+3)) (in N);
    per-branch steps S_branch ~ N^{a} / T^{delta} with a = S/(2K),
    delta = (K+1)/2; effective delta at a = 1/2: K(K+3)/(2S) - 1."""
    K = len(I)
    S = sum(abs(i) for i in I)
    e = Fraction(S, K * (K + 3))
    return {"I": list(I), "K": K, "S": S, "exponent": e,
            "a": Fraction(S, 2 * K), "delta": Fraction(K + 1, 2),
            "delta_eff": Fraction(K * (K + 3), 2 * S) - 1}


def all_supports(exps: Sequence[int] = ALL_EXPONENTS) -> List[Tuple[int, ...]]:
    out = []
    for K in range(1, len(exps) + 1):
        out.extend(combinations(exps, K))
    return out


def exponent_with_guesses(I: Sequence[int], mu: Fraction, kappa: Fraction = None,
                          D_exponent: Fraction = Fraction(0)) -> Fraction:
    """Heuristic total exponent (in N) when the residue class of p mod M = N^mu
    is guessed (phi(M) ~ M runs), the target lattice gains index M^kappa, and
    every ring operation costs N^{D_exponent} (cyclotomic: D = phi(m) ~ m).
    Derivation in the notes: e = S/(K(K+3)) + mu (1 - 2(K+kappa)/(K(K+3))) + D."""
    st = support_stats(I)
    K, S = st["K"], st["S"]
    if kappa is None:
        kappa = Fraction(K * (K + 1), 2)
    return Fraction(S, K * (K + 3)) + mu * (1 - Fraction(2 * (K + kappa), K * (K + 3))) + D_exponent


# --------------------------------------------------------------------------
# Wronskian (nondegeneracy) -- used in the proof of Lemma C
# --------------------------------------------------------------------------
def falling(i: int, k: int) -> int:
    return prod(i - j for j in range(k))


def wronskian_det(I: Sequence[int]) -> int:
    """det[(i)_k]_{k=1..K, i in I}; equals prod(i) * prod_{i<i'} (i'-i)."""
    K = len(I)
    Mx = flint.fmpz_mat([[falling(i, k) for i in I] for k in range(1, K + 1)])
    return int(Mx.det())


def wronskian_formula(I: Sequence[int]) -> int:
    I = sorted(I)
    v = prod(I)
    for a in range(len(I)):
        for b in range(a + 1, len(I)):
            v *= I[b] - I[a]
    return v


# --------------------------------------------------------------------------
# integer-valued lattice on a residue class (exact HNF computation)
# --------------------------------------------------------------------------
def _basis_values_mod(I: Sequence[int], N: int, r: int, M: int, e: int) -> List[List[int]]:
    """rows: for x in [0, M^(e-1)), the vector (f_i(p) mod M^e)_i with
    p = r + M x, f_i = p^i (i>0) or (N p^{-1})^{|i|} (i<0)."""
    Me = M ** e
    rows = []
    for x in range(M ** (e - 1)):
        pp = (r + M * x) % Me
        qq = (N * pow(pp, -1, Me)) % Me
        rows.append([pow(pp, i, Me) if i > 0 else pow(qq, -i, Me) for i in I])
    return rows


def class_kernel(I: Sequence[int], N: int, r: int, M: int, e: int) -> flint.fmpz_mat:
    """Lattice {a in Z^K : sum a_i f_i(p) == const (mod M^e) on p == r (M)}
    (HNF basis, rows)."""
    K = len(I)
    Me = M ** e
    vals = _basis_values_mod(I, N, r, M, e)
    diffs = [[(v[i] - vals[0][i]) % Me for i in range(K)] for v in vals[1:]]
    # kernel of a -> (diffs . a) mod Me: compute via HNF of [A^T | I] trick
    # rows of [ C ; Me*I ] where C has columns = conditions -> use nullspace mod
    # Standard: lattice = { a : D a == 0 mod Me }. Build
    #   [ I_K   D^T ]
    #   [ 0   Me I_c]
    # and HNF; rows whose condition part is 0 form the kernel.
    if diffs:
        # replace the (many) conditions by an equivalent K-row basis:
        # row space of [diffs; Me*I] over Z is a full-rank lattice in Z^K
        Hc = flint.fmpz_mat(diffs + [[Me if i == j else 0 for j in range(K)]
                                     for i in range(K)]).hnf().tolist()
        diffs = [[int(v) for v in row] for row in Hc if any(int(v) for v in row)]
    c = len(diffs)
    if c == 0:
        return flint.fmpz_mat([[1 if i == j else 0 for j in range(K)] for i in range(K)])
    big = []
    for i in range(K):
        big.append([1 if i == j else 0 for j in range(K)] + [diffs[t][i] for t in range(c)])
    for t in range(c):
        big.append([0] * K + [Me if s == t else 0 for s in range(c)])
    # put condition columns first so HNF eliminates them
    perm = [[row[K + t] for t in range(c)] + row[:K] for row in big]
    H = flint.fmpz_mat(perm).hnf()
    rows = H.tolist()
    ker = [[int(v) for v in row[c:]] for row in rows if all(int(v) == 0 for v in row[:c])]
    ker = [row for row in ker if any(row)]
    return flint.fmpz_mat(ker)


def class_lattice_index(I: Sequence[int], N: int, r: int, M: int, emax: int = None) -> Dict:
    """Exact index [L_{r,M} : Z^K] where L = sum_e M^{-e} Ker_e (e <= emax),
    in coefficient coordinates a (target sum a_i f_i).  Returns
    {'log_index': kappa (as Fraction of log M if a pure power), 'index': int}.
    """
    K = len(I)
    if emax is None:
        emax = K
    Mtop = M ** emax
    gens = []
    for e in range(0, emax + 1):
        ker = class_kernel(I, N, r, M, e) if e >= 1 else flint.fmpz_mat(
            [[1 if i == j else 0 for j in range(K)] for i in range(K)])
        scale = M ** (emax - e)
        for row in ker.tolist():
            gens.append([int(v) * scale for v in row])
    H = flint.fmpz_mat(gens).hnf()
    rows = [r_ for r_ in H.tolist() if any(int(v) for v in r_)]
    det = abs(int(flint.fmpz_mat(rows).det()))
    # L = (1/Mtop) * span(gens); index over Z^K = Mtop^K / det
    num = Mtop ** K
    assert num % det == 0
    index = num // det
    kappa = None
    k, t = 0, index
    while t % M == 0:
        t //= M
        k += 1
    if t == 1:
        kappa = k
    return {"index": index, "kappa": kappa, "basis_scaled": rows, "scale": Mtop}
