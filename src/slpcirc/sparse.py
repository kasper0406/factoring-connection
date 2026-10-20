"""Lemma S1 / Theorem SP helpers: sparse (t-nomial) root bounds on cyclic groups and the
per-prime dichotomy check. Exact arithmetic only."""
from __future__ import annotations
import math
from sympy import factorint, primerange


def iroot_floor(n: int, k: int) -> int:
    """floor(n^(1/k)) for n >= 0, k >= 1 (exact)."""
    if k == 1:
        return n
    r = int(round(n ** (1.0 / k)))
    while r ** k > n:
        r -= 1
    while (r + 1) ** k <= n:
        r += 1
    return r


def s1_K(ell: int, t: int) -> int:
    return iroot_floor(ell - 1, t - 1)


def s1_coset_bound(ell: int, t: int) -> int:
    """Max zeros of a t-term polynomial (exponents distinct mod ell) on one coset of mu_ell."""
    K = s1_K(ell, t)
    if K < 1:
        return ell
    return min(ell, 2 * (-(-ell // K) - 1))


def ell0(t: int, lam_num: int, lam_den: int) -> int:
    """Threshold of Theorem SP: ell0 = (floor(2/lam) + 1)^(t-1) + 1 (lam = lam_num/lam_den).
    For every ell >= ell0, K_ell(t) = floor((ell-1)^(1/(t-1))) >= floor(2/lam) + 1 > 2/lam."""
    return (2 * lam_den // lam_num + 1) ** (t - 1) + 1


def reduced_terms(coeffs_exps, n, p):
    """Merge exponents mod n, coefficients mod p; drop zeros. Returns dict e mod n -> c mod p."""
    d = {}
    for c, e in coeffs_exps:
        k = e % n
        d[k] = (d.get(k, 0) + c) % p
    return {k: v for k, v in d.items() if v}


def count_zeros_Fp_star(coeffs_exps, p):
    """Exact count of u in F_p^* with sum c u^e = 0 (brute force over a generator's powers)."""
    n = p - 1
    red = reduced_terms(coeffs_exps, n, p)
    if not red:
        return n
    g = primitive_root(p)
    # evaluate at g^j: sum c * g^(e j)
    terms = [(c, pow(g, e, p)) for e, c in red.items()]
    cnt = 0
    vals = [c for c, _ in terms]
    for _ in range(n):
        if sum(vals) % p == 0:
            cnt += 1
        vals = [(v * m) % p for v, (_, m) in zip(vals, terms)]
    return cnt


def primitive_root(p: int) -> int:
    if p == 2:
        return 1
    fs = list(factorint(p - 1))
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fs):
            return g
    raise ValueError


def Lambda(exps):
    """prod_{i<j} (e_i - e_j) over distinct exponents (nonzero integer)."""
    es = sorted(set(exps))
    L = 1
    for i in range(len(es)):
        for j in range(i + 1, len(es)):
            L *= es[j] - es[i]
    return L


def dichotomy_violations(coeffs_exps, p, lam_num, lam_den, roots_in_Fp):
    """Theorem SP per-prime check. roots_in_Fp = #roots of F in F_p. Returns list of primes ell >= ell0
    dividing p-1 but not Lambda, when density >= lam; [] otherwise (no violation)."""
    if roots_in_Fp * lam_den < lam_num * p:
        return []
    t = len(set(e for _, e in coeffs_exps))
    # density on F_p^* >= (lam p - 1)/(p-1); use lam' = lam - 1/p (rational)
    ln, ld = lam_num * p - lam_den, lam_den * p
    if ln <= 0:
        return []
    l0 = ell0(max(t, 2), ln, ld)
    L = Lambda([e for _, e in coeffs_exps])
    return [ell for ell in factorint(p - 1) if ell >= l0 and L % ell != 0]
