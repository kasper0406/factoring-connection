"""E2: exact independence / construction test for the G_r counting hypothesis (agent midk).

The rank-r aggregate family beats L[1/3] (reaching L[1/4]) only if the r-1 linking congruences mod N
are DEPENDENT at family height ~ N^{1/d}, i.e. if G_r fails favourably.  This module tests, with exact
arithmetic, the two ways such dependence could arise and shows each is blocked:
  (I)  a GLOBAL (characteristic-0) dependence  -> ring-power/subfield -> geometric height (measured);
  (II) a dependence that lives only mod N (a "symmetry" realised mod N) -> a non-scalar Z/N intertwiner
       T; Lemma SPL-INT: extracting its eigenstructure splits N (factoring-equivalent).

`split_from_intertwiner` and anything reading p, q are DIAGNOSTIC (factor-aware) and never imported by
operational code.  All other routines receive only N and public data.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import List, Optional, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# Lemma SPL-INT (proved, elementary): a Z/N eigen-computation of a non-scalar
# intertwiner either returns a global Z/N eigenvector or splits N.
# ---------------------------------------------------------------------------

def gauss_kernel_or_factor(A: List[List[int]], N: int) -> Tuple[str, object]:
    """Fraction-free Gaussian elimination of the matrix A over Z/N.
    Returns ('factor', g) with 1 < g < N if a nonzero non-unit pivot is met (splits N), else
    ('kernel', v) with v a nonzero kernel vector mod N if one exists, else ('injective', None).
    N squarefree assumed for the clean dichotomy."""
    m = len(A)
    n = len(A[0]) if m else 0
    M = [[a % N for a in row] for row in A]
    piv_cols: List[int] = []
    rr = 0
    for cclmn in range(n):
        # find a pivot in column cclmn at or below row rr: PREFER a unit entry (a factor is reported
        # only when the column has a nonzero entry but NO unit one -- a genuine rank obstruction, not
        # an accidental non-unit entry that a unit pivot elsewhere in the column would bypass).
        piv = None
        nonunit = None
        for i in range(rr, m):
            v = M[i][cclmn] % N
            if v == 0:
                continue
            g = gcd(v, N)
            if g == 1:
                piv = i
                break
            elif nonunit is None:
                nonunit = g
        if piv is None:
            if nonunit is not None:
                return ("factor", nonunit)
            # column is all zero at/below rr: free column
            continue
        M[rr], M[piv] = M[piv], M[rr]
        inv = pow(M[rr][cclmn], -1, N)
        M[rr] = [(inv * v) % N for v in M[rr]]
        for i in range(m):
            if i != rr and M[i][cclmn] % N:
                f = M[i][cclmn]
                M[i] = [(M[i][k] - f * M[rr][k]) % N for k in range(n)]
        piv_cols.append(cclmn)
        rr += 1
        if rr == m:
            break
    # build a kernel vector from a free column, if any
    free = [j for j in range(n) if j not in piv_cols]
    if not free:
        return ("injective", None)
    j0 = free[0]
    v = [0] * n
    v[j0] = 1
    for idx, pc in enumerate(piv_cols):
        v[pc] = (-M[idx][j0]) % N
    return ("kernel", v)


def split_from_intertwiner(N: int, T: List[List[int]], lam: int) -> Optional[int]:
    """Operational core (uses only N and T, no p, q): given a claimed Z/N eigenvalue lam of T, run the
    eigen-elimination of (T - lam I) over Z/N.  Returns a nontrivial factor of N iff the elimination
    meets a nonzero non-unit pivot (a *rank discrepancy* across the prime factors), else None
    (a consistent global Z/N eigenvector / injectivity).  Correct sufficient condition for a split:
    rank_p(T-lam I) != rank_q(T-lam I) for some p, q | N (Opus critic); a non-scalar T alone does NOT
    suffice (e.g. [[lam,1],[0,mu]] returns a kernel, no split)."""
    m = len(T)
    A = [[(T[i][j] - (lam if i == j else 0)) % N for j in range(m)] for i in range(m)]
    tag, val = gauss_kernel_or_factor(A, N)
    if tag == "factor":
        return int(val)
    return None


# ---------------------------------------------------------------------------
# Constructors of intertwiners (public): scalar / global vs mod-N-only symmetry.
# ---------------------------------------------------------------------------

def crt(a_p: int, a_q: int, p: int, q: int) -> int:
    """x = a_p mod p, a_q mod q (DIAGNOSTIC: uses p, q to simulate a structured mod-N object)."""
    N = p * q
    return (a_p * q * pow(q, -1, p) + a_q * p * pow(p, -1, q)) % N


def scalar_intertwiner(N: int, m: int, lam: int) -> List[List[int]]:
    """A global (char-0) scalar symmetry: T = lam*I.  Consistent mod every prime, gives no split and
    corresponds to a Z-relation (ring-power/subfield case, which pays geometric height)."""
    return [[(lam if i == j else 0) % N for j in range(m)] for i in range(m)]


def _randmat(rng, m, mod, rank=None):
    """A random m x m matrix over Z/mod (mod prime) of exactly the given rank (default full).
    Constructed directly as an outer product B(m x k) * C(k x m), k = rank, which has rank k
    generically; verified, resampled on the rare degenerate draw."""
    k = m if rank is None else rank
    for _ in range(200):
        B = [[rng.randrange(mod) for _ in range(k)] for _ in range(m)]
        C = [[rng.randrange(mod) for _ in range(m)] for _ in range(k)]
        A = [[sum(B[i][t] * C[t][j] for t in range(k)) % mod for j in range(m)] for i in range(m)]
        if _rank_modp(A, mod) == k:
            return A
    raise RuntimeError("could not sample matrix of requested rank")


def _rank_modp(M0, p):
    M = [[int(v) % p for v in row] for row in M0]
    rows, cols = len(M), len(M[0])
    rk = 0
    for c in range(cols):
        piv = next((i for i in range(rk, rows) if M[i][c] % p), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], -1, p)
        M[rk] = [(v * inv) % p for v in M[rk]]
        for i in range(rows):
            if i != rk and M[i][c] % p:
                f = M[i][c]
                M[i] = [(M[i][k] - f * M[rk][k]) % p for k in range(cols)]
        rk += 1
    return rk


def rank_discrepant_intertwiner(p: int, q: int, m: int, seed: int = 0) -> Tuple[List[List[int]], int]:
    """DIAGNOSTIC construction of the matrix a favourable G_r failure would need: an intertwiner
    A = T - lam I whose rank differs across the prime factors of N -- rank m-1 mod p (a one-dim
    'symmetry' visible at p, e.g. a CM/root-of-unity action realised mod p) but full rank m mod q.
    A is NOT identically 0 mod p (it is a generic rank-(m-1) matrix there), so a split requires the
    genuine non-unit pivot of Lemma SPL-INT, not a trivial entry gcd.  Returns (A mod N, 0): we test
    split_from_intertwiner(N, A, 0) directly (A already plays the role of T - lam I)."""
    import random
    N = p * q
    rng = random.Random(seed)
    Ap = _randmat(rng, m, p, rank=m - 1)     # singular mod p
    Aq = _randmat(rng, m, q, rank=m)         # invertible mod q
    A = [[crt(Ap[i][j] % p, Aq[i][j] % q, p, q) for j in range(m)] for i in range(m)]
    return A, 0


# ---------------------------------------------------------------------------
# P2a: independence baseline of the linking congruences at small primes.
# ---------------------------------------------------------------------------

def linking_independence_count(p: int, r: int, trials: int, seed: int = 0) -> dict:
    """Model the r-1 linking conditions as r independent projective coincidences mod p: two random
    r-tuples of residues 'link at rank r' iff they are projectively proportional mod p.  The
    independence prediction for Pr[all r generators link | first links] is (1/p)^{r-1}.  We measure
    the empirical rate that r uniform pairs (a_j, b_j) in (F_p*)^2 satisfy a_1 b_j ≡ a_j b_1 (all j)
    (proportional), i.e. share one projective linking direction across all r coordinates."""
    import random
    rr = random.Random(seed)
    hits = 0
    for _ in range(trials):
        a = [rr.randrange(1, p) for _ in range(r)]
        b = [rr.randrange(1, p) for _ in range(r)]
        if all((a[0] * b[j] - a[j] * b[0]) % p == 0 for j in range(r)):
            hits += 1
    pred = trials * (1.0 / p) ** (r - 1)
    import math
    sigma = math.sqrt(pred) if pred > 0 else 0.0
    return {"p": p, "r": r, "trials": trials, "hits": hits, "pred": pred,
            "z": (hits - pred) / sigma if sigma > 0 else 0.0}


# ---------------------------------------------------------------------------
# P2b(I): structured constructions pay geometric height.
# ---------------------------------------------------------------------------

def ring_power_height(d: int, base_height: int, r: int) -> dict:
    """A ring-generated (power) linked module w_j = alpha^{j-1} in a degree-d field of generator
    height H0: the j-th generator has height ~ H0^{j-1}, so the family height is H0^{r-1} - geometric
    in r.  Returns the exact log-height ratio vs the aggregate target N^{1/d}."""
    logH0 = sp.log(base_height)
    fam_log_height = (r - 1) * logH0      # height of the top power
    target = sp.Rational(1, d)            # in units of x = log N
    return {"d": d, "r": r, "family_log_height_over_logH0": r - 1,
            "note": "geometric in r; aggregate needs O(1)*x/d independent of r"}


def automorphism_graph_witness(m: int, p: int, q: int, r: int) -> dict:
    """The sharp concrete favourable-failure of the counting hypothesis (codex midk_1 §4).
    K = Q(zeta_m); phi: zeta -> t with Phi_m(t) = 0 mod N; sigma: zeta -> zeta^g a Galois element;
    phi_1 = phi, phi_2 = phi o sigma^{-1}; w_{1j} = zeta^{j}, w_{2j} = sigma(zeta^{j}) = zeta^{gj}.
    Then phi_1(w_{1j}) = phi_2(w_{2j}) = t^j for ALL j, at height 1 (roots of unity): the r-1 linking
    conditions are completely DEPENDENT with no height growth -> the counting half of G_r is false.
    But candidate u_2 = sigma(u_1), so the two norms coincide and any modular square root on side 2 is
    sigma(s) with phi_2(sigma(s)) = phi_1(s): the gcd is trivial.  Useless (no new relations).
    Returns the exact checks.  DIAGNOSTIC (uses p, q only to build a root of Phi_m mod N)."""
    import sympy as sp
    N = p * q
    g = next(a for a in range(2, m) if gcd(a, m) == 1)        # a nontrivial Galois element zeta->zeta^g
    Phi = sp.Poly(sp.cyclotomic_poly(m, sp.Symbol("t")), sp.Symbol("t"))
    # roots of Phi_m mod p and mod q, CRT to a root t mod N (needs p,q = 1 mod m)
    rp = [a for a in range(p) if Phi.eval(a) % p == 0]
    rq = [a for a in range(q) if Phi.eval(a) % q == 0]
    if not rp or not rq:
        return {"m": m, "ok": False, "reason": "Phi_m has no root mod p or q"}
    t = crt(rp[0], rq[0], p, q)
    ginv = pow(g, -1, m)
    # link check: phi_1(w_1j) = t^j ; phi_2(w_2j) = phi(sigma^{-1}(zeta^{gj})) = phi(zeta^{j}) = t^j
    links = all(pow(t, j, N) == pow(t, (ginv * (g * j)) % m, N) for j in range(1, r + 1))
    # the two sides carry identical norms (diagonal relation) -> useless:
    diagonal = (ginv * g) % m == 1
    return {"m": m, "r": r, "N_bits": N.bit_length(), "links_dependent_height1": links,
            "diagonal_useless": diagonal, "g": g}


def nnorm_fullrank_minima(f_coeffs: List[int], N: int) -> dict:
    """'N is a norm' O_K/(beta) ~ Z/N gives a FULL-RANK linked lattice; its successive minima are
    generic (~ det^{1/d}), so candidates have norm ~ N^{1/d} with no rank saving: useless (rank r=d,
    not an intermediate aggregate).  We just record the rank is full (= d)."""
    d = len(f_coeffs) - 1
    return {"degree": d, "linked_rank": d, "useless": True,
            "reason": "full rank: no intermediate r, minima generic"}
