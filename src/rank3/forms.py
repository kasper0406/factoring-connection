"""Exact tools for rank-r linked modules (agent rank3).

DIAGNOSTIC / harness code: it takes polynomials, module generators and N, never factors N inside any algorithm,
and implements no factoring method. Real/complex numbers appear only as certified `arb`/`acb` balls (python-flint),
used to certify integer values (Delta_F) or rigorous inequalities; every reported integer is exact.

Conventions.
- A polynomial is a coefficient list low -> high. f, g are squarefree integer polynomials defining the two sides
  (K = Q[x]/f, K' = Q[x]/g, possibly reducible / non-monic in E2, monic irreducible in E1).
- A linked rank-r module is given by polynomial lists P = [P_1..P_r] (K side) and Q = [Q_1..Q_r] (K' side):
  u_j = P_j(alpha), u'_j = Q_j(beta), with P_j(m) = Q_j(m) mod N for a common root m of f, g mod N
  (the linking v_j = P_j(m) mod N). Power basis: P_j = Q_j = x^(j-1).
- Candidate a in Z^r: A_K = sum a_j P_j, A_K' = sum a_j Q_j; W(a) = Res(f, A_K), W'(a) = Res(g, A_K'); for monic f
  this is N_K(sum a_j u_j). F = W W' is the decomposable form whose linear factors are
  L_sigma = (P_j(alpha_sigma))_j (times lc(f)^{deg/d} normalisations for non-monic f; E1 uses monic f, g).
- Delta_F = prod over ORDERED r-tuples of distinct factors of det(L_{i1},...,L_{ir}) (Bean-Thunder 1996).
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction as Fr
from typing import Dict, List, Sequence, Tuple

import flint
from flint import acb, arb, ctx, fmpz_poly

Poly = Sequence[int]


def fpoly(c: Poly) -> fmpz_poly:
    return fmpz_poly([int(x) for x in c])


def resultant(f: Poly, g: Poly) -> int:
    return int(fpoly(f).resultant(fpoly(g)))


def discriminant(f: Poly) -> int:
    return int(fpoly(f).discriminant())


def poly_comb(polys: Sequence[Poly], a: Sequence[int]) -> List[int]:
    n = max(len(p) for p in polys)
    out = [0] * n
    for aj, p in zip(a, polys):
        if aj:
            for k, c in enumerate(p):
                out[k] += aj * c
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def predicted_valuation(r: int, D: int) -> int:
    """Number of ordered r-tuples of distinct indices from D containing two given indices: r(r-1)(D-2)!/(D-r)!."""
    return r * (r - 1) * math.factorial(D - 2) // math.factorial(D - r)


def kappa(r: int, D: int) -> Fr:
    """N-exponent of Theorem R1: P(r,D) * (D-r)!/D! = r(r-1)/(D(D-1))."""
    return Fr(r * (r - 1), D * (D - 1))


# ----------------------------------------------------------------------------- certified constants


def bean_thunder_C(r: int, prec: int = 256) -> arb:
    """C_r = (2/r) prod_{k=1}^{r-1} [B(1/(r+1),k/(r+1)) + B((r-k)/(r+1),k/(r+1)) + B((r-k)/(r+1),1/(r+1))],
    as an arb ball (Beta via Gamma). C_2 = 3B(1/3,1/3) = 15.8997 (Bean)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        def B(x: Fr, y: Fr) -> arb:
            ax, ay = arb(x.numerator) / x.denominator, arb(y.numerator) / y.denominator
            return ax.gamma() * ay.gamma() / (ax + ay).gamma()
        s = arb(2) / r
        for k in range(1, r):
            s *= (B(Fr(1, r + 1), Fr(k, r + 1)) + B(Fr(r - k, r + 1), Fr(k, r + 1))
                  + B(Fr(r - k, r + 1), Fr(1, r + 1)))
        return s
    finally:
        ctx.prec = old


def log_C_upper_elementary(r: int) -> float:
    """Elementary bound log C_r <= (r-1) log(12(r+1)) + log(2/r) (from B(x,y) <= 2/x + 2/y on (0,1])."""
    return (r - 1) * math.log(12 * (r + 1)) + math.log(2 / r)


# ----------------------------------------------------------------------------- linear forms and Delta_F


def acb_roots(f: Poly) -> List[acb]:
    rts = []
    for z, mult in fpoly(f).complex_roots():
        if mult != 1:
            raise ValueError("polynomial must be squarefree")
        rts.append(z)
    return rts


def _ev(p: Poly, z: acb) -> acb:
    s = acb(0)
    for c in reversed(list(p)):
        s = s * z + int(c)
    return s


def linear_forms(f: Poly, P: Sequence[Poly], g: Poly, Q: Sequence[Poly]) -> Tuple[List[List[acb]], int]:
    """Rows L_sigma = (P_j(alpha_sigma))_j for roots of f, then (Q_j(beta_tau))_j for roots of g.
    Requires monic f, g (then prod_i L_i(a) = Res(f, A_K) Res(g, A_K') exactly)."""
    if f[-1] != 1 or g[-1] != 1:
        raise ValueError("linear_forms requires monic f, g")
    rows = [[_ev(p, z) for p in P] for z in acb_roots(f)] + [[_ev(q, z) for q in Q] for z in acb_roots(g)]
    return rows, len(acb_roots(f))


def _det(M: List[List[acb]]) -> acb:
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    s = acb(0)
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = M[0][j] * _det(minor)
        s = s + term if j % 2 == 0 else s - term
    return s


def delta_F(f: Poly, P: Sequence[Poly], g: Poly, Q: Sequence[Poly], r: int,
            prec0: int = 512, prec_max: int = 1 << 17) -> Dict[str, object]:
    """Certified exact Delta_F.

    Delta_F = prod over ordered tuples = sign * prod_{|S|=r} det(S)^{r!}. We certify the rational integer
    Q2 := (prod_S det S)^2 (Galois-invariant: Galois permutes subsets and changes each det by a sign), then
    |Delta_F| = |Q2|^{r!/2}. Returns Q2 (exact int, possibly 0), r!/2, and the precision used."""
    prec = prec0
    while prec <= prec_max:
        old = ctx.prec
        ctx.prec = prec
        try:
            rows, _ = linear_forms(f, P, g, Q)
            D = len(rows)
            prod = acb(1)
            for S in itertools.combinations(range(D), r):
                prod *= _det([rows[i] for i in S])
            q2 = prod * prod
            z = q2.real.unique_fmpz()
            im_ok = q2.imag.contains(0)
        finally:
            ctx.prec = old
        if z is not None and im_ok:
            return {"Q2": int(z), "exp": math.factorial(r) // 2, "prec": prec, "D": D}
        prec *= 2
    raise RuntimeError("could not certify Delta_F")


def vp(n: int, p: int) -> int:
    if n == 0:
        return 10 ** 9
    c = 0
    while n % p == 0:
        n //= p
        c += 1
    return c


# ----------------------------------------------------------------------------- counts and checks


def box_points(r: int, E: int):
    return itertools.product(range(-E, E + 1), repeat=r)


def F_value(f: Poly, P: Sequence[Poly], g: Poly, Q: Sequence[Poly], a: Sequence[int]) -> Tuple[int, int]:
    AK = poly_comb(P, a)
    AK2 = poly_comb(Q, a)
    if all(c == 0 for c in AK) or all(c == 0 for c in AK2):
        return 0, 0
    return resultant(f, AK), resultant(g, AK2)


def box_abs_values(f, P, g, Q, r: int, E: int) -> List[int]:
    """Sorted |F(a)| over the box |a_j| <= E, zeros excluded (exact)."""
    vals = []
    for a in box_points(r, E):
        w, w2 = F_value(f, P, g, Q, a)
        v = abs(w * w2)
        if v:
            vals.append(v)
    vals.sort()
    return vals


def count_bound_terms(r: int, D: int, E: int, X: int, info: Dict[str, object], prec: int = 256) -> Tuple[arb, int]:
    """Theorem-R3 bound pieces: main = C_r |Delta_F|^{-(D-r)!/D!} X^{r/D} (arb ball), err = r(D+1)(2E+1)^{r-1}."""
    old = ctx.prec
    ctx.prec = prec
    try:
        C = bean_thunder_C(r, prec)
        e_delta = arb(info["exp"] * math.factorial(D - r)) / math.factorial(D)
        logmain = C.log() - e_delta * arb(abs(int(info["Q2"]))).log() + arb(r) / D * arb(X).log()
        main = logmain.exp()
    finally:
        ctx.prec = old
    return main, r * (D + 1) * (2 * E + 1) ** (r - 1)


def liouville_check(f, P, g, Q, r: int, info: Dict[str, object], prec: int = 512) -> Dict[str, float]:
    """Proposition R4: M := prod_i max_j |L_i(e_j)| >= (2^r/(r! C_r))^{D/r} |Delta_F|^{(D-r)!/(r (D-1)!)}.
    Returns certified lower bound of log M, upper bound of log RHS, and whether lower >= upper."""
    old = ctx.prec
    ctx.prec = prec
    try:
        rows, _ = linear_forms(f, P, g, Q)
        D = len(rows)
        logM = arb(0)
        for row in rows:
            mags = [abs(z) for z in row]
            mx = mags[0]
            for m_ in mags[1:]:
                mx = mx.max(m_)
            logM += mx.log()
        C = bean_thunder_C(r, prec)
        lhs_const = arb(D) / r * (arb(2) ** r / math.factorial(r) / C).log()
        e = arb(info["exp"] * math.factorial(D - r)) / (r * math.factorial(D - 1))
        rhs = lhs_const + e * arb(abs(int(info["Q2"]))).log()
        lo = logM.lower()
        hi = rhs.upper()
        return {"logM_lower": float(lo), "logRHS_upper": float(hi), "holds": bool(lo >= hi)}
    finally:
        ctx.prec = old
