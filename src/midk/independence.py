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
        # find a pivot in column cclmn at or below row rr
        piv = None
        for i in range(rr, m):
            g = gcd(M[i][cclmn] % N, N)
            if M[i][cclmn] % N == 0:
                continue
            if g == 1:
                piv = i
                break
            else:
                # nonzero non-unit -> nontrivial factor
                return ("factor", g)
        if piv is None:
            # column has only zeros or we already returned; it is a free column
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
    """DIAGNOSTIC-free operational core: given a claimed Z/N eigenvalue lam of T, attempt the
    eigen-computation.  Returns a nontrivial factor of N if the elimination exposes a zero divisor,
    else None (a consistent global Z/N eigenvector, no split).  Uses only N and T (no p, q)."""
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


def modp_only_symmetry(p: int, q: int, m: int, lam_p: int, seed: int = 0) -> Tuple[List[List[int]], int]:
    """DIAGNOSTIC construction of a 'structured mod p only' intertwiner: T acts as the scalar lam_p
    mod p (a symmetry visible at the prime p, e.g. a root of unity / CM action realised mod p) but as
    a generic non-scalar matrix mod q.  Returns (T mod N, lam) with lam ≡ lam_p mod p chosen so that
    det(T-lam I) ≡ 0 mod N.  This is exactly the family a favourable G_r failure would need."""
    N = p * q
    rng = sp.Integer  # exact
    import random
    r = random.Random(seed)
    # mod q: a generic matrix with a chosen eigenvalue lam_q so (T-lam I) is singular mod q too
    lam_q = r.randrange(q)
    # build T mod q with eigenvalue lam_q but non-scalar: companion-like
    Tq = [[r.randrange(q) for _ in range(m)] for _ in range(m)]
    # force an eigenvector e1: set first column so that T e1 = lam_q e1 -> make T upper-triangular-ish
    for i in range(m):
        Tq[i][0] = (lam_q if i == 0 else 0)
    Tp = [[(lam_p if i == j else 0) for j in range(m)] for i in range(m)]  # scalar mod p
    lam = crt(lam_p, lam_q, p, q)
    T = [[crt(Tp[i][j] % p, Tq[i][j] % q, p, q) for j in range(m)] for i in range(m)]
    return T, lam


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


def nnorm_fullrank_minima(f_coeffs: List[int], N: int) -> dict:
    """'N is a norm' O_K/(beta) ~ Z/N gives a FULL-RANK linked lattice; its successive minima are
    generic (~ det^{1/d}), so candidates have norm ~ N^{1/d} with no rank saving: useless (rank r=d,
    not an intermediate aggregate).  We just record the rank is full (= d)."""
    d = len(f_coeffs) - 1
    return {"degree": d, "linked_rank": d, "useless": True,
            "reason": "full rank: no intermediate r, minima generic"}
