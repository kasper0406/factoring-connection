"""Integer analogues of the systematic-equation relation family.  Exact arithmetic.

OPERATIONAL-SAFE: nothing here factors N or uses its factors.  Smoothness is tested against
the fixed primorial of B (Bernstein's test), which is a test of the candidate values, not of N.

E2  Cunningham-type "F_{1^n}" V-family (the closest integer instance of the BGJT axioms):
    N = Phi_n(m), so m^n == 1 (mod N): mu_n = <m> is an explicit "constant field" inside Z/N
    and X^n == 1 plays the role of X^q == h0/h1.  For f = A x^{n+1} + B x^n + C x + D,
        f(m) == r_f := (A + C) m + (B + D)   (mod N).
    A relation is a pair (f(m), r_f) with both sides B-smooth.  The BGJT analogue would need
    f to split over Z into linear factors with small values; Descartes forbids that (E3).
E3  Descartes count: completely split members of V_n(Z) of small height.
E4  Cyclotomic calibration: m^n - 1 = prod_{d|n} Phi_d(m); smoothness versus the per-piece
    random baseline (cyclotomic splitting is a size reduction, not smoothness by construction).
"""
from __future__ import annotations

import random
from fractions import Fraction
from typing import Dict, List

import gmpy2
from flint import fmpz_poly
from sympy import cyclotomic_poly, divisors, primerange
from sympy.abc import x as _x


# ---------------------------------------------------------------------------------------
# exact smoothness
_PRIMORIAL: Dict[int, gmpy2.mpz] = {}


def primorial(B: int) -> gmpy2.mpz:
    if B not in _PRIMORIAL:
        P = gmpy2.mpz(1)
        for p in primerange(2, B + 1):
            P *= p
        _PRIMORIAL[B] = P
    return _PRIMORIAL[B]


def is_smooth(v: int, B: int) -> bool:
    """Exact: |v| is B-smooth (v != 0).  v | P^(2^k) with 2^k >= log2|v|  <=>  B-smooth."""
    v = abs(int(v))
    if v == 0:
        raise ValueError("0")
    if v <= B:
        return True
    v = gmpy2.mpz(v)
    r = primorial(B) % v
    k = 1
    while (1 << k) < v.bit_length():
        k += 1
    for _ in range(k):
        if r == 0:
            return True
        r = r * r % v
    return r == 0


def random_baseline(bits: List[int], B: int, samples: int, seed: int) -> Dict[int, Fraction]:
    """Exact empirical B-smooth rate of uniform integers in [2^(l-1), 2^l), per bit length l.
    For l <= floor(log2 B) the rate is 1 exactly (every such integer is <= B)."""
    out = {}
    for l in sorted(set(bits)):
        if (1 << l) - 1 <= B:
            out[l] = Fraction(1)
            continue
        rng = random.Random(f"systeq|baseline|{B}|{l}|{seed}")
        hits = sum(is_smooth(rng.randrange(1 << (l - 1), 1 << l), B) for _ in range(samples))
        out[l] = Fraction(hits, samples)
    return out


def phi_n(n: int) -> List[int]:
    return [int(c) for c in reversed(cyclotomic_poly(n, _x, polys=True).all_coeffs())]


def ev(coeffs: List[int], m: int) -> int:
    r = 0
    for c in reversed(coeffs):
        r = r * m + c
    return r


# ---------------------------------------------------------------------------------------
# E2
def vfamily_candidates(n: int, m: int, H: int, count: int, seed: int):
    """Seeded sample of (A,B,C,D) in [-H,H]^4, A != 0, r_f != 0; yields (v=f(m), r_f)."""
    N = ev(phi_n(n), m)
    rng = random.Random(f"systeq|vfam|{n}|{m}|{H}|{seed}")
    mn = m ** n
    out = []
    while len(out) < count:
        A, Bc, C, D = (rng.randint(-H, H) for _ in range(4))
        if A == 0:
            continue
        r = (A + C) * m + (Bc + D)
        v = A * mn * m + Bc * mn + C * m + D
        if r == 0 or v == 0:
            continue
        out.append((v, r))
    # exactness check of the congruence on every candidate
    assert all((v - r) % N == 0 for v, r in out)
    return N, out


def vfamily_yield(n: int, m: int, H: int, count: int, B: int, base: Dict[int, Fraction],
                  seed: int, exclude_degenerate: bool = False) -> Dict:
    """exclude_degenerate (POST-HOC, critic fix): drop candidates with A m + B = 0 (then
    f(m) = r_f exactly) or C m + D = 0 (then f(m) = m^n r_f): their 'relations' only restate
    m^n == 1.  The candidate list is regenerated with the same seed and filtered."""
    N, cands = vfamily_candidates(n, m, H, count, seed)
    if exclude_degenerate:
        mn = m ** n
        cands = [(v, r) for v, r in cands if v != r and v != mn * r]
    obs = 0
    pred = Fraction(0)
    for v, r in cands:
        lv, lr = abs(v).bit_length(), abs(r).bit_length()
        pred += base[lv] * base[lr]
        if is_smooth(v, B) and is_smooth(r, B):
            obs += 1
    return {"n": n, "m": m, "N_bits": N.bit_length(), "H": H, "B": B, "candidates": len(cands),
            "min_v_bits": min(abs(v).bit_length() for v, _ in cands),
            "max_r_abs": max(abs(r) for _, r in cands),
            "observed": obs, "predicted": pred, "G": (Fraction(obs) / pred) if pred else None,
            "v_bits": sorted({abs(v).bit_length() for v, _ in cands})}


def vfamily_bits(n: int, m: int, H: int, count: int, seed: int) -> List[int]:
    _, cands = vfamily_candidates(n, m, H, count, seed)
    return [abs(v).bit_length() for v, _ in cands] + [abs(r).bit_length() for _, r in cands]


# ---------------------------------------------------------------------------------------
# E3
def descartes_count(n: int, H: int) -> Dict:
    """Exhaustive over A x^{n+1} + B x^n + C x + D, coefficients in [-H, H], A != 0.
    Counts members that split completely over Q, separating the degenerate x^n (A x + B)."""
    split_nondeg, split_deg, total = 0, 0, 0
    examples = []
    for A in range(-H, H + 1):
        if A == 0:
            continue
        for Bc in range(-H, H + 1):
            for C in range(-H, H + 1):
                for D in range(-H, H + 1):
                    total += 1
                    coeffs = [D, C] + [0] * (n - 2) + [Bc, A]
                    f = fmpz_poly(coeffs)
                    _, facs = f.factor()
                    if all(g.degree() == 1 for g, _ in facs):
                        if C == 0 and D == 0:
                            split_deg += 1
                        else:
                            split_nondeg += 1
                            if len(examples) < 5:
                                examples.append((A, Bc, C, D))
    return {"n": n, "H": H, "total": total, "split_nondegenerate": split_nondeg,
            "split_degenerate": split_deg, "examples": examples}


# ---------------------------------------------------------------------------------------
# E4
def cyclo_calibration(n: int, m_lo: int, m_hi: int, B: int, base: Dict[int, Fraction]) -> Dict:
    ds = divisors(n)
    polys = {d: phi_n(d) for d in ds}
    obs = 0
    pred = Fraction(0)
    for m in range(m_lo, m_hi):
        pieces = [ev(polys[d], m) for d in ds]
        assert eval_prod(pieces) == m ** n - 1
        pr = Fraction(1)
        for v in pieces:
            pr *= base[abs(v).bit_length()]
        pred += pr
        if all(is_smooth(v, B) for v in pieces):
            obs += 1
    return {"n": n, "m_range": [m_lo, m_hi], "B": B, "observed": obs, "predicted": pred,
            "R_piece": (Fraction(obs) / pred) if pred else None}


def cyclo_bits(n: int, m_lo: int, m_hi: int) -> List[int]:
    ds = divisors(n)
    return [abs(ev(phi_n(d), m)).bit_length() for m in range(m_lo, m_hi) for d in ds]


def eval_prod(xs):
    r = 1
    for v in xs:
        r *= v
    return r
