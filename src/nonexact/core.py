"""Exact helpers for the nonexact agent.  Exact integer / rational arithmetic only.

Polynomials are python-flint fmpz_poly or lists of ints (lowest degree first).
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt
from typing import Iterable, Sequence

import numpy as np
from flint import fmpz_mat, fmpz_poly, nmod_poly


# --------------------------------------------------------------------------- polynomials
def from_roots(roots: Iterable[int]) -> fmpz_poly:
    p = fmpz_poly([1])
    for a in roots:
        p *= fmpz_poly([-a, 1])
    return p


def coeffs(p: fmpz_poly, n: int | None = None) -> list[int]:
    c = [int(x) for x in p.coeffs()]
    if n is not None:
        c = c + [0] * (n - len(c))
    return c


def height(c: Sequence[int]) -> int:
    return max((abs(int(x)) for x in c), default=0)


def centred(x: int, N: int) -> int:
    x %= N
    return x - N if x > N // 2 else x


def red_mod(P: fmpz_poly, F: fmpz_poly) -> fmpz_poly:
    """Remainder of P modulo the MONIC polynomial F, in Z[X]."""
    assert int(F.coeffs()[-1]) == 1
    return P % F


def red_mod_N_centred(P: fmpz_poly, F: fmpz_poly, N: int) -> list[int]:
    d = F.degree()
    return [centred(c, N) for c in coeffs(red_mod(P, F), d)]


def split_over_Z_small(c: Sequence[int], B: int) -> bool:
    """True iff the nonzero integer polynomial c equals (const) * prod (X - a) with all a integral, |a| <= B.
    Constants count as split (empty product)."""
    p = fmpz_poly([int(x) for x in c])
    if p == 0:
        return False
    if p.degree() == 0:
        return True
    _, facs = p.factor()
    for f, _e in facs:
        if f.degree() != 1:
            return False
        a, b = int(f.coeffs()[0]), int(f.coeffs()[1])
        if b not in (1, -1):
            return False
        root = -a * b  # f = b X + a, root -a/b
        if abs(root) > B:
            return False
    return True


def int_roots_count(c: Sequence[int]) -> int:
    """Number of distinct integer roots of a nonzero integer polynomial."""
    p = fmpz_poly([int(x) for x in c])
    if p == 0 or p.degree() <= 0:
        return 0
    _, facs = p.factor()
    return sum(1 for f, _ in facs if f.degree() == 1 and int(f.coeffs()[1]) in (1, -1))


# --------------------------------------------------------------------------- lattices
def lll_rows(rows: Sequence[Sequence[int]]) -> list[list[int]]:
    """Exact LLL (flint, delta = 0.99 default) of the row lattice; returns reduced rows."""
    M = fmpz_mat([[int(x) for x in r] for r in rows])
    R = M.lll()
    return [[int(R[i, j]) for j in range(R.ncols())] for i in range(R.nrows())]


def kernel_lattice_mod_N(A: Sequence[Sequence[int]], N: int) -> list[list[int]]:
    """LLL-reduced basis of L = {z in Z^n : z A == 0 (mod N)}, A an n x d integer matrix.

    Embedding: rows [e_i | W A_i] (i < n) and [0 | W N e_j] (j < d) with W large; after LLL the rows with
    zero tail span L (W exceeds every needed norm).  Exact."""
    n, d = len(A), len(A[0])
    W = N * N * (n + 1)
    rows = []
    for i in range(n):
        rows.append([1 if k == i else 0 for k in range(n)] + [W * (int(A[i][j]) % N) for j in range(d)])
    for j in range(d):
        rows.append([0] * n + [W * N if k == j else 0 for k in range(d)])
    R = lll_rows(rows)
    basis = [r[:n] for r in R if all(x == 0 for x in r[n:])]
    assert len(basis) == n, "kernel extraction failed"
    return basis


def det_abs(rows: Sequence[Sequence[int]]) -> int:
    return abs(int(fmpz_mat([[int(x) for x in r] for r in rows]).det()))


# --------------------------------------------------------------------------- exact Gram-Schmidt + enumeration
def gram_schmidt(rows: Sequence[Sequence[int]]):
    """Exact (Fraction) Gram-Schmidt: returns mu (lower-triangular) and Bsq (squared norms of b_i*)."""
    n = len(rows)
    bstar: list[list[Fraction]] = []
    mu = [[Fraction(0)] * n for _ in range(n)]
    Bsq: list[Fraction] = []
    for i in range(n):
        v = [Fraction(x) for x in rows[i]]
        for j in range(i):
            num = sum(Fraction(rows[i][k]) * bstar[j][k] for k in range(len(v)))
            mu[i][j] = num / Bsq[j]
            v = [v[k] - mu[i][j] * bstar[j][k] for k in range(len(v))]
        bstar.append(v)
        Bsq.append(sum(x * x for x in v))
    return mu, Bsq


def _int_range(center: Fraction, rad_sq: Fraction) -> tuple[int, int]:
    """All integers x with (x - center)^2 <= rad_sq, as [lo, hi] (possibly empty: lo > hi). Exact."""
    if rad_sq < 0:
        return 1, 0
    # integer sqrt bound on the rational rad_sq = p/q
    p, q = rad_sq.numerator, rad_sq.denominator
    r = isqrt(p * q) // q + 1  # r >= sqrt(rad_sq)
    import math
    lo = math.floor(center - r)
    hi = math.ceil(center + r)
    while lo <= hi and (lo - center) ** 2 > rad_sq:
        lo += 1
    while hi >= lo and (hi - center) ** 2 > rad_sq:
        hi -= 1
    return lo, hi


def enumerate_ball(rows: Sequence[Sequence[int]], target: Sequence[Fraction], R2: Fraction, counter: dict):
    """All lattice vectors y = sum x_i b_i with ||y - target||^2 <= R2 (exact Fincke-Pohst / Schnorr-Euchner
    without pruning).  counter['nodes'] counts enumeration-tree nodes (the cost unit of Prop LIN)."""
    n = len(rows)
    mu, Bsq = gram_schmidt(rows)
    # coordinates of target in the GS basis: t = sum tau_i b_i*  (only the projection matters)
    m = len(rows[0])
    bstar = []
    for i in range(n):
        v = [Fraction(x) for x in rows[i]]
        for j in range(i):
            v = [v[k] - mu[i][j] * bstar[j][k] for k in range(m)]
        bstar.append(v)
    tau = [sum(Fraction(target[k]) * bstar[i][k] for k in range(m)) / Bsq[i] for i in range(n)]
    out = []
    x = [0] * n

    def rec(i: int, partial: Fraction):
        counter["nodes"] = counter.get("nodes", 0) + 1
        c = tau[i] - sum(mu[k][i] * x[k] for k in range(i + 1, n))
        rem = (R2 - partial) / Bsq[i]
        lo, hi = _int_range(c, rem)
        for xi in range(lo, hi + 1):
            x[i] = xi
            add = (xi - c) ** 2 * Bsq[i]
            if i == 0:
                out.append(list(x))
            else:
                rec(i - 1, partial + add)
        x[i] = 0

    rec(n - 1, Fraction(0))
    vecs = []
    for xs in out:
        vecs.append([sum(xs[i] * rows[i][k] for i in range(n)) for k in range(m)])
    return vecs


# --------------------------------------------------------------------------- Prop LIN: linearised lattice
def s_multisets_e(B: int, s: int):
    """All s-multisets of [-B, B] with their elementary symmetric vectors (e_1..e_s)."""
    out = []
    for S in combinations_with_replacement(range(-B, B + 1), s):
        c = coeffs(from_roots(S), s + 1)  # c[i] = coefficient of X^i, monic: c[s] = 1
        # P_S = X^s - e1 X^{s-1} + e2 X^{s-2} - ...  => e_i = (-1)^i c[s-i]
        e = [(-1) ** i * c[s - i] for i in range(1, s + 1)]
        out.append((S, e))
    return out


def lin_weights(c: int, mu: int, N: int, s: int) -> list[int]:
    """w_0 = c mu^s, w_i = (-1)^i c mu^{s-i} (mod N), so c*P_S(mu) == w_0 + sum_i e_i w_i (mod N)."""
    return [(pow(-1, i) * c * pow(mu, s - i, N)) % N for i in range(0, s + 1)]


def e_bounds(B: int, s: int) -> list[int]:
    return [comb(s, i) * B ** i for i in range(1, s + 1)]


def lin_box_volume(B: int, s: int, V: int) -> int:
    v = 2 * V + 1
    for b in e_bounds(B, s):
        v *= 2 * b + 1
    return v


def lin_split_hits(w: Sequence[int], N: int, V: int, E: np.ndarray) -> int:
    """Number of s-multisets (rows of E = their e-vectors, int64) with |centred(w0 + E.w)| <= V.
    Requires N < 2^40 and |e| * N < 2^62 (asserted by the caller)."""
    assert N < 2 ** 40
    acc = np.full(E.shape[0], w[0] % N, dtype=np.int64)
    for i in range(E.shape[1]):
        acc = np.mod(acc + _mulmod_vec(np.mod(E[:, i], N), w[i + 1], N), N)
    cen = np.where(acc > N // 2, acc - N, acc)
    return int(np.count_nonzero(np.abs(cen) <= V))


def lin_box_points(w: Sequence[int], N: int, B: int, s: int, V: int) -> int:
    """Exact number of (e_1..e_s, v) in the box |e_i| <= C(s,i)B^i, |v| <= V with
    w0 + sum e_i w_i == v (mod N).  Vectorised over (e_1..e_{s-1}); solved for e_s."""
    bnd = e_bounds(B, s)
    ws = w[s] % N
    if gcd(ws, N) != 1:
        raise ValueError("non-unit weight (reveals a factor)")
    inv = pow(ws, -1, N)
    # base over e_1..e_{s-1}
    base = np.array([w[0] % N], dtype=np.int64)
    for i in range(1, s):
        rng = np.arange(-bnd[i - 1], bnd[i - 1] + 1, dtype=np.int64)
        term = _mulmod_vec(np.mod(rng, N), w[i], N)
        base = np.mod(base[:, None] + term[None, :], N).ravel()
    bs = bnd[s - 1]
    total = 0
    if 2 * V + 1 <= 2 * bs + 1:
        # loop over v: e_s == inv*(v - base) mod N, need |centred(e_s)| <= bs
        for v in range(-V, V + 1):
            es = _mulmod_vec(np.mod(v - base, N), inv, N)
            cen = np.where(es > N // 2, es - N, es)
            total += int(np.count_nonzero(np.abs(cen) <= bs))
    else:
        for e in range(-bs, bs + 1):
            val = np.mod(base + (e * ws) % N, N)  # python-int scalar < N
            cen = np.where(val > N // 2, val - N, val)
            total += int(np.count_nonzero(np.abs(cen) <= V))
    return total


def _mulmod_vec(a: np.ndarray, c: int, N: int) -> np.ndarray:
    """(a * c) mod N for int64 arrays with entries in [0, N), N < 2^40, without overflow
    (split c into 20-bit limbs)."""
    c %= N
    lo = c & ((1 << 20) - 1)
    hi = c >> 20
    r = np.mod(a * lo, N)
    t = np.mod(a * hi, N)
    t = np.mod(t * (1 << 20), N)
    return np.mod(r + t, N)


def lin_cost_model(N: int, B: int, s: int, V: int, n_ms: int | None = None) -> dict:
    """Prop LIN model: expected box points p, expected split hits h per lattice, cost per hit
    max(1, p)/h (setup >= 1 per lattice, 1 per enumerated point), and brute-force cost per hit N/(2V+1)."""
    if n_ms is None:
        n_ms = comb(2 * B + s, s)
    p = Fraction(lin_box_volume(B, s, V), N)
    h = Fraction(n_ms * (2 * V + 1), N)
    return {"p": p, "h": h, "cost_per_hit": max(Fraction(1), p) / h, "brute_cost_per_hit": Fraction(N, 2 * V + 1)}


def lin_lattice_enumerate(w: Sequence[int], N: int, B: int, s: int, V: int, counter: dict) -> list[list[int]]:
    """OPERATIONAL-style (uses only N, mu-derived weights and public bounds): enumerate all box points
    (e_1..e_s, v) of the coset {w0 + sum e_i w_i - v == 0 (mod N)} with an LLL-reduced basis of the scaled
    lattice and exact Fincke-Pohst on the circumscribed ball.  Returns the box points."""
    bnd = e_bounds(B, s) + [V]
    n = s + 1
    L = 1
    for b in bnd:
        L = L * b // gcd(L, b)
    scale = [L // b for b in bnd]  # coordinate i scaled by L/b_i: box -> cube [-L, L]^n
    # homogeneous lattice {(e, v) : sum e_i w_i - v == 0 mod N}: basis e_i -> (unit_i, w_i), plus (0.., N)
    rows = []
    for i in range(s):
        r = [0] * n
        r[i] = 1
        r[s] = w[i + 1] % N
        rows.append(r)
    rN = [0] * n
    rN[s] = N
    rows.append(rN)
    rows = [[r[k] * scale[k] for k in range(n)] for r in rows]
    red = lll_rows(rows)
    red = [r for r in red if any(r)]
    # coset representative: e = 0, v = w0  => point z0 = (0,..,0, w0) ; we need y with y == z0 mod lattice,
    # i.e. enumerate lattice vectors y with ||(z0 - y) scaled|| small: target = z0 scaled.
    t = [Fraction(0)] * s + [Fraction((w[0] % N) * scale[s])]
    R2 = Fraction(n * L * L)
    vecs = enumerate_ball(red, t, R2, counter)
    pts = []
    for y in vecs:
        z = [t[k] - y[k] for k in range(n)]  # scaled coset point
        if all(abs(z[k]) <= L for k in range(n)):
            pts.append([int(z[k] / scale[k]) for k in range(n)])
    return pts


# --------------------------------------------------------------------------- Prop IH
def ih_relations(A: Sequence[int], A2: Sequence[int], r: Sequence[int], r2: Sequence[int], N: int, k: int):
    """Prop IH.  P = P_A, P' = P_A2 in one fibre (deg(P - P') < k), with P == r and P' == r2 mod (N, F),
    deg F >= k.  Then D := P - P' satisfies D == delta := r - r2 coefficientwise mod N; E := (D - delta)/N.
    Returns the list of integer relations (b, u_b, v_b, nontrivial) at roots b of P and P' that are not common:
      b root of P' : u_b = P(b) = D(b),  v_b = delta(b);
      b root of P  : u_b = P'(b) = -D(b), v_b = -delta(b).
    Each satisfies u_b == v_b (mod N); nontrivial iff E(b) != 0 (iff u_b != v_b)."""
    P, P2 = from_roots(A), from_roots(A2)
    D = P - P2
    assert D.degree() < k
    delta = fmpz_poly([int(a) - int(b) for a, b in zip(list(r) + [0] * k, list(r2) + [0] * k)])
    Ecoef = coeffs(D - delta, k)
    assert all(c % N == 0 for c in Ecoef), "D is not == delta mod N"
    E = fmpz_poly([c // N for c in Ecoef])
    rels = []
    common = set(A) & set(A2)
    for b in sorted(set(A2) - common):
        u, v = int(P(b)), int(delta(b))
        rels.append((b, u, v, int(E(b)) != 0))
    for b in sorted(set(A) - common):
        u, v = int(P2(b)), -int(delta(b))
        rels.append((b, u, v, int(E(b)) != 0))
    return rels, E


def ih_distinct_lower_bound(M: int, k: int, degD: int, n_supp: int | None = None) -> int:
    """REGISTERED (historical, used for the frozen P9(e) scoring; REFUTED as a theorem: it ignores the sign
    identification (u, v) ~ (-u, -v)).  Kept unchanged so the registered scoring is reproducible."""
    if n_supp is None:
        n_supp = M - k + 2
    good = n_supp - (k - 1)
    if degD <= 0:
        return 1 if good >= 1 else 0
    return max(0, -(-good // degD))


def ih_hypothesis(delta: Sequence[int], B: int, N: int) -> bool:
    """Prop IH's value-smallness hypothesis: ||delta||_inf * sum_{l<k} B^l < N (k = len(delta))."""
    k = len(delta)
    return height(delta) * sum(B ** l for l in range(k)) < N


def ih_good_points_bound(M: int, k: int) -> int:
    """Prop IH(b), under ih_hypothesis: at least M - 2k + 3 good points (non-common roots b with E(b) != 0)."""
    return max(0, M - 2 * k + 3)


def ih_distinct_bound(M: int, k: int, degD: int) -> int:
    """Prop IH(c), corrected (post hoc), under ih_hypothesis: distinct nontrivial relations up to sign
    >= ceil((M - 2k + 3)/(2 deg D)) for deg D >= 1; exactly 1 if deg D = 0 (then delta is constant)."""
    good = ih_good_points_bound(M, k)
    if good <= 0:
        return 0
    if degD <= 0:
        return 1
    return -(-good // (2 * degD))


# --------------------------------------------------------------------------- Prop LAT: multiplier lattices
def multiplier_lattice(P: fmpz_poly, F: fmpz_poly, N: int, j: int, k: int) -> list[list[int]]:
    """LLL-reduced basis of L_P = {(Q, r) : deg Q <= j, deg r < k, Q P - r in (N, F)} (F monic, k <= deg F).
    Coordinates: (q_0..q_j, r_0..r_{k-1})."""
    d = F.degree()
    assert k <= d
    A = []
    XiP = P
    for i in range(j + 1):
        A.append([c % N for c in coeffs(red_mod(XiP, F), d)])
        XiP = XiP * fmpz_poly([0, 1])
    for i in range(k):
        row = [0] * d
        row[i] = (-1) % N
        A.append(row)
    return kernel_lattice_mod_N(A, N)


def quadratic_presentation(mu: int, mu2: int, N: int) -> fmpz_poly:
    """Monic F = (X - mu)(X - mu2) reduced into [0, N) coefficientwise: F(mu) == F(mu2) == 0 (mod N)."""
    f0 = (mu * mu2) % N
    f1 = (-(mu + mu2)) % N
    return fmpz_poly([f0, f1, 1])


def classify_vector(vec: Sequence[int], j: int, k: int, B: int) -> dict:
    """For a basis vector (Q, r): Q_const (deg 0), Q_split (const * split with roots in [-B, B]),
    r_split likewise, useful = Q_split and r_split and Q != 0 and r != 0."""
    q = list(vec[: j + 1])
    r = list(vec[j + 1:])
    qn = any(q)
    rn = any(r)
    qs = qn and split_over_Z_small(q, B)
    rs = rn and split_over_Z_small(r, B)
    qdeg = max((i for i, c in enumerate(q) if c), default=-1)
    return {"Q_nonzero": qn, "r_nonzero": rn, "Q_deg": qdeg, "Q_split": qs, "r_split": rs,
            "useful": bool(qn and rn and qs and rs), "Q_int_roots": int_roots_count(q) if qn else 0}
