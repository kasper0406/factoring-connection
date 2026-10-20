"""Exact binary-form tools for the generalized resultant barrier (agent resbarrier).

Everything here is exact integer arithmetic (python ints / python-flint). These are DIAGNOSTIC / harness
routines: they take polynomial pairs and N, never factor N, and implement no factoring algorithm.

Conventions. A univariate integer polynomial is a coefficient list c[0..d] (low to high). Its homogenisation
is F(a, b) = sum_i c[i] a^i b^(d-i), so a - b*m corresponds to g = [-m, 1]  (g(x) = x - m) and NFS candidates
a - b*alpha have norm F(a, b) (up to sign).
"""
from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction as Fr
from math import comb, gcd, isqrt
from typing import Dict, Iterable, List, Sequence, Tuple

import flint

Poly = Sequence[int]

# Bean's constant 3*B(1/3,1/3) = 15.8997...; we use the rational upper bound 16 (see notes, Lemma C2).
BEAN_UPPER = Fr(16)


def fpoly(c: Poly) -> flint.fmpz_poly:
    return flint.fmpz_poly([int(x) for x in c])


def polymul(f: Poly, g: Poly) -> List[int]:
    return [int(x) for x in (fpoly(f) * fpoly(g)).coeffs()]


def resultant(f: Poly, g: Poly) -> int:
    return int(fpoly(f).resultant(fpoly(g)))


def discriminant(f: Poly) -> int:
    return int(fpoly(f).discriminant())


def binary_disc(f: Poly) -> int:
    """Discriminant of the binary form F(a,b) = b^d f(a/b) of degree d = len(f)-1.

    Equals disc(f) when the leading coefficient is nonzero (the only case used here)."""
    if f[-1] == 0:
        raise ValueError("degree must be exact")
    return discriminant(f)


def height(f: Poly) -> int:
    return max(abs(int(x)) for x in f)


def hom_eval(f: Poly, a: int, b: int) -> int:
    """F(a,b) = sum c_i a^i b^(d-i), exact."""
    d = len(f) - 1
    return sum(int(f[i]) * a ** i * b ** (d - i) for i in range(d + 1))


def base_m_pair(N: int, d: int, m: int) -> Tuple[List[int], List[int]]:
    """Base-m expansion N = sum c_i m^i (0 <= c_i < m, top digit free): f = sum c_i x^i, g = x - m."""
    c = []
    r = N
    for _ in range(d):
        c.append(r % m)
        r //= m
    c.append(r)
    return c, [-m, 1]


def snfs_pair(e: int, d: int) -> Tuple[int, List[int], List[int]]:
    """N = 2^e - 1, f = x^d - 2^(d k - e), g = x - 2^k with k = ceil(e/d): f(2^k) = 2^(dk-e) N."""
    N = (1 << e) - 1
    k = -(-e // d)
    f = [-(1 << (d * k - e))] + [0] * (d - 1) + [1]
    return N, f, [-(1 << k), 1]


# ----------------------------------------------------------------------------- Lemma RES / C2 checks


def res_check(f: Poly, g: Poly, N: int) -> Dict[str, int]:
    """Lemma RES (snfsrig): N | Res(f,g) != 0 and |Res| <= ||f||_2^deg g ||g||_2^deg f (checked squared)."""
    R = resultant(f, g)
    df, dg = len(f) - 1, len(g) - 1
    n2f = sum(x * x for x in f)
    n2g = sum(x * x for x in g)
    return {
        "res": R,
        "divisible": int(R != 0 and R % N == 0),
        "hadamard_ok": int(R * R <= n2f ** dg * n2g ** df),
    }


def disc_product_floor(f: Poly, g: Poly, N: int) -> Dict[str, int]:
    """|disc(F1 F2)| = |disc F1| |disc F2| Res^2 >= N^2 when N | Res and both discriminants are nonzero."""
    fg = polymul(f, g)
    dF = binary_disc(f)
    # disc of a linear form is 1; general identity disc(PQ) = disc(P) disc(Q) Res(P,Q)^2
    dG = binary_disc(g) if len(g) > 2 else 1
    R = resultant(f, g)
    prod = dF * dG * R * R
    full = binary_disc(fg)
    return {"disc_prod": full, "identity_ok": int(abs(full) == abs(prod)), "ge_N2": int(abs(full) >= N * N)}


def box_values(F: Poly, E: int) -> List[int]:
    """Sorted |F(a,b)| over |a| <= E, 0 <= b <= E, excluding zeros (exact)."""
    d = len(F) - 1
    vals = []
    for b in range(0, E + 1):
        bpow = [1] * (d + 1)
        for k in range(1, d + 1):
            bpow[k] = bpow[k - 1] * b
        coeffs = [F[i] * bpow[d - i] for i in range(d + 1)]  # polynomial in a
        for a in range(-E, E + 1):
            s = 0
            for c in reversed(coeffs):  # Horner in a
                s = s * a + c
            if s:
                vals.append(abs(s))
    vals.sort()
    return vals


def count_le(vals: List[int], X: int) -> int:
    return bisect_right(vals, X)


def c2_bound_holds(count: int, D: int, E: int, X: int, disc_floor: int, c: Fr = BEAN_UPPER) -> bool:
    """Exact test of count <= c * disc_floor^(-1/(D(D-1))) * X^(2/D) + 2(D+1)(2E+1).

    With S = count - err: S <= 0, or (S/c)^(D(D-1)) * disc_floor <= X^(2(D-1))."""
    S = count - 2 * (D + 1) * (2 * E + 1)
    if S <= 0:
        return True
    k = D * (D - 1)
    lhs = (Fr(S) / c) ** k * disc_floor
    return lhs <= Fr(X) ** (2 * (D - 1))


def c2_saturation(count: int, D: int, E: int, X: int, disc_floor: int) -> float:
    """(count - err)_+ / main term, as a float (diagnostic summary only)."""
    S = count - 2 * (D + 1) * (2 * E + 1)
    if S <= 0:
        return 0.0
    from math import log, exp
    main_log = log(float(BEAN_UPPER)) - log(disc_floor) / (D * (D - 1)) + 2 * log(X) / D
    return exp(log(S) - main_log)


# ----------------------------------------------------------------------------- Lemma KL (rank-r Liouville)


def kernel_short_vector(v: Sequence[int]) -> List[int]:
    """A nonzero integer vector a with sum a_j v_j = 0 of small size (LLL on the kernel lattice)."""
    r = len(v)
    # kernel basis via the lattice [I | K v] reduced; flint fmpz_mat.lll
    K = 1
    for x in v:
        K = max(K, abs(x))
    K = K * (r + 1) * 4
    rows = [[1 if i == j else 0 for j in range(r)] + [K * v[i]] for i in range(r)]
    M = flint.fmpz_mat(rows).lll()
    best = None
    for i in range(M.nrows()):
        row = [int(M[i, j]) for j in range(r + 1)]
        if row[-1] == 0 and any(row[:r]):
            vec = row[:r]
            if best is None or max(map(abs, vec)) < max(map(abs, best)):
                best = vec
    if best is None:
        raise ValueError("no kernel vector found")
    return best


def lemma_kl_check(v: Sequence[int], W2_eval, N: int, d2: int, H2: int) -> Dict[str, int]:
    """Lemma KL: if the rational side is W1(a) = sum a_j v_j (primitive v) and W2 is a degree-d2 form in r
    variables with W2 = 0 mod N on ker W1 (linking), then for the shortest kernel vector a, W2(a) != 0
    gives N <= |W2(a)| <= C(d2+r-1, r-1) H2 |a|_inf^d2. Returns the exact quantities."""
    a = kernel_short_vector(v)
    w = W2_eval(a)
    r = len(v)
    amax = max(map(abs, a))
    return {
        "a_inf": amax,
        "W2a": w,
        "divisible": int(w % N == 0),
        "nonzero": int(w != 0),
        "ge_N": int(abs(w) >= N),
        "upper_ok": int(abs(w) <= comb(d2 + r - 1, r - 1) * H2 * amax ** d2),
    }
