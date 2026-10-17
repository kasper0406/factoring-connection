"""Exhaustive check of the GPT-6.1-sol round-4 claims (notes/codex_r4_sos_classgroup.md §2–§3).  DIAGNOSTIC.

Everything here knows p and q.  Exact integer arithmetic only (Python ints).

Claims checked, for every pair of distinct odd primes p < q < BOUND:

  R4-1  exact identity      s = eps(N+1) - 16*eps*u*v,  p = eps + 4u, q = eps + 4v   (needs p = q mod 4)
  R4-2  mod-16 table        s mod 16 = eps(N+1) mod 16, i.e. the printed 4x2 table
  R4-3  mod-32 formula      s = eps(N+1) - 16*eps*[(2/p) = -1]*[(2/q) = -1]  (mod 32)
        and, for N = 1 mod 8: (2/p) = (2/q), and s mod 32 is eps(N+1) - 16*eps*[2 is a non-square mod N]
  R4-4  4-rank criterion    for N = 1 mod 8:  2 is a square mod N  <=>  r4 Cl(-8N) > 0 (eps = +1),
                                                                         r4 Cl(-4N) > 0 (eps = -1)

The class groups are computed by brute force: enumerate the reduced primitive positive definite forms of
discriminant D, compose by Dirichlet/Shanks composition, and count elements with x^2 = 1 and x^4 = 1.
For a finite abelian group G, |G[2]| = 2^{r2} and |G[4]| = 2^{r2 + r4}.

Extras (also printed in the r4 note, checked here because they are cheap):
  X1  mod-64 formula for N = 5 mod 8
  X2  the public class T = [(2, 2, (N+1)/2)] is a square in Cl(-4N)  <=>  (2/p) = (2/q) = +1   (p = q mod 4)
  X3  r4 = 2 - rank_F2(M) for the printed genus matrices M_{-4N}, M_{-8N}            (N = 1 mod 8)
  X4  Milovic Prop. 1 as quoted: eps = +1:  r4 Cl(-4N) = 2  <=>  p = q = 1 mod 8 and (p/q) = 1   (N = 1 mod 8)
  X5  r2 Cl(D) = 2 for D in {-4N, -8N} (genus theory sanity check)
  X6  the "wrong" discriminant (Cl(-4N) when eps = +1, Cl(-8N) when eps = -1): how often it agrees

Run:  .venv/bin/python -m src.repair.verify_r4 --bound 300 --out experiments/raw/verify_r4.json
"""
from __future__ import annotations

import argparse
import json
import time
from math import gcd, isqrt
from typing import Dict, List, Tuple

Form = Tuple[int, int, int]


# ---------------------------------------------------------------------------
# elementary number theory
# ---------------------------------------------------------------------------
def odd_primes_below(bound: int) -> List[int]:
    sieve = bytearray([1]) * bound
    sieve[0:2] = b"\x00\x00"
    for i in range(2, isqrt(bound - 1) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(3, bound) if sieve[i]]


def legendre(a: int, p: int) -> int:
    """Legendre symbol (a/p) for an odd prime p, by Euler's criterion."""
    t = pow(a % p, (p - 1) // 2, p)
    return -1 if t == p - 1 else t


def kronecker(D: int, n: int) -> int:
    """Kronecker symbol (D/n) for n >= 1 (used only for the analytic class number check)."""
    if n == 1:
        return 1
    res = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0:
            return 0
        res *= 1 if D % 8 in (1, 7) else -1
    # Jacobi symbol (D/n), n odd
    a, m = D % n, n
    while a:
        while a % 2 == 0:
            a //= 2
            if m % 8 in (3, 5):
                res = -res
        a, m = m, a
        if a % 4 == 3 and m % 4 == 3:
            res = -res
        a %= m
    return res if m == 1 else 0


def class_number_analytic(D: int) -> int:
    """h(D) = -(1/|D|) sum_{a=1}^{|D|} (D/a) a, valid for fundamental D < -4.  Independent of the form code."""
    n = -D
    s = sum(kronecker(D, a) * a for a in range(1, n + 1))
    assert s % n == 0
    return -s // n


# ---------------------------------------------------------------------------
# binary quadratic forms (positive definite, discriminant D < 0)
# ---------------------------------------------------------------------------
def _egcd(a: int, b: int) -> Tuple[int, int, int]:
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        k = a // b
        a, b = b, a - k * b
        x0, x1 = x1, x0 - k * x1
        y0, y1 = y1, y0 - k * y1
    return a, x0, y0


def reduce_form(f: Form) -> Form:
    a, b, c = f
    D = b * b - 4 * a * c
    while True:
        # normalize: -a < b <= a
        if not (-a < b <= a):
            k = (a - b) // (2 * a)
            b += 2 * a * k
            c = (b * b - D) // (4 * a)
        if a > c:
            a, b, c = c, -b, a
            continue
        if a == c and b < 0:
            b = -b
        return a, b, c


def compose(f: Form, g: Form) -> Form:
    """Dirichlet composition (Shanks' formula), then reduction."""
    a1, b1, c1 = f
    a2, b2, c2 = g
    D = b1 * b1 - 4 * a1 * c1
    assert b2 * b2 - 4 * a2 * c2 == D
    s = (b1 + b2) // 2
    g1, x1, y1 = _egcd(a1, a2)
    e, x2, w = _egcd(g1, s)
    u, v = x2 * x1, x2 * y1
    a3 = a1 * a2 // (e * e)
    num = u * a1 * b2 + v * a2 * b1 + w * (b1 * b2 + D) // 2
    assert num % e == 0
    b3 = (num // e) % (2 * a3)
    assert (b3 * b3 - D) % (4 * a3) == 0
    c3 = (b3 * b3 - D) // (4 * a3)
    return reduce_form((a3, b3, c3))


def identity_form(D: int) -> Form:
    return (1, 0, -D // 4) if D % 4 == 0 else (1, 1, (1 - D) // 4)


def reduced_forms(D: int) -> List[Form]:
    """All reduced primitive forms (a, b, c): |b| <= a <= c, b >= 0 if |b| = a or a = c."""
    assert D < 0 and D % 4 in (0, 1)
    out = []
    amax = isqrt(-D // 3)
    for a in range(1, amax + 1):
        for b in range(D % 2, a + 1, 2):  # b = D (mod 2)
            t = b * b - D
            if t % (4 * a):
                continue
            c = t // (4 * a)
            if c < a:
                continue
            for bb in ((b,) if (b == 0 or b == a or a == c) else (b, -b)):
                if gcd(gcd(a, abs(bb)), c) == 1:
                    out.append((a, bb, c))
    return out


def two_four_torsion(D: int) -> Dict[str, object]:
    """Brute force: h, |G[2]|, |G[4]|, r2, r4, and the set of squares G^2."""
    forms = reduced_forms(D)
    one = identity_form(D)
    n2 = n4 = 0
    squares = set()
    for f in forms:
        f2 = compose(f, f)
        squares.add(f2)
        if f2 == one:
            n2 += 1
            n4 += 1
        elif compose(f2, f2) == one:
            n4 += 1
    r2 = n2.bit_length() - 1
    r24 = n4.bit_length() - 1
    assert n2 == 1 << r2 and n4 == 1 << r24, (D, n2, n4)
    return {"h": len(forms), "G2": n2, "G4": n4, "r2": r2, "r4": r24 - r2, "squares": squares}


def rank_f2(M: List[List[int]]) -> int:
    rows = [r[:] for r in M]
    rank, ncol = 0, len(rows[0])
    for col in range(ncol):
        piv = next((i for i in range(rank, len(rows)) if rows[i][col] % 2), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i][col] % 2:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


# ---------------------------------------------------------------------------
# the checks
# ---------------------------------------------------------------------------
TABLE16 = {1: (2, 14), 5: (6, 10), 9: (10, 6), 13: (14, 2)}  # N mod 16 -> (eps=+1, eps=-1), as printed


def check_pair(p: int, q: int, class_groups: bool = True) -> Dict[str, object]:
    """All checks for one pair p < q with p = q (mod 4).  Returns booleans (True = claim holds)."""
    N, s = p * q, p + q
    assert p % 4 == q % 4
    eps = 1 if p % 4 == 1 else -1
    u, v = (p - eps) // 4, (q - eps) // 4
    bp, bq = legendre(2, p), legendre(2, q)
    two_sq = bp == 1 and bq == 1  # 2 is a square mod N  <=>  square mod p and mod q (CRT)
    rec: Dict[str, object] = {"p": p, "q": q, "N": N, "eps": eps, "N_mod_8": N % 8}

    rec["R4_1_identity"] = s == eps * (N + 1) - 16 * eps * u * v
    rec["R4_2_table16"] = s % 16 == TABLE16[N % 16][0 if eps == 1 else 1] and s % 16 == eps * (N + 1) % 16
    ind = ((1 - bp) // 2) * ((1 - bq) // 2)
    rec["R4_3_mod32_formula"] = s % 32 == (eps * (N + 1) - 16 * eps * ind) % 32
    if N % 8 == 1:
        rec["R4_3_bp_eq_bq"] = bp == bq
        rec["R4_3_bit_is_2_square"] = s % 32 == (eps * (N + 1) - 16 * eps * (0 if two_sq else 1)) % 32
    if N % 8 == 5:
        rec["X1_bp_ne_bq"] = bp != bq
        rec["X1_mod32_no_new_bit"] = s % 32 == eps * (N + 1) % 32
        r = p if bp == 1 else q
        rec["X1_mod64_formula"] = s % 64 == (eps * (N + 1) - 32 * eps * (((r - eps) // 8) % 2)) % 64

    if not class_groups:
        return rec
    # Cl(-4N): always computed (needed for X2); Cl(-8N): computed for N = 1 mod 8.
    G4N = two_four_torsion(-4 * N)
    T = reduce_form((2, 2, (N + 1) // 2))
    rec["X2_T_square_iff_2_square"] = (T in G4N["squares"]) == two_sq
    rec["X5_r2_minus4N"] = G4N["r2"] == 2
    rec["h_minus4N"], rec["r4_minus4N"] = G4N["h"], G4N["r4"]
    if N % 8 == 1:
        G8N = two_four_torsion(-8 * N)
        rec["X5_r2_minus8N"] = G8N["r2"] == 2
        rec["h_minus8N"], rec["r4_minus8N"] = G8N["h"], G8N["r4"]
        sel, other = (G8N, G4N) if eps == 1 else (G4N, G8N)
        rec["R4_4_criterion"] = (sel["r4"] > 0) == two_sq
        rec["X6_wrong_disc_agrees"] = (other["r4"] > 0) == two_sq
        beta = (1 - bp) // 2
        eta = (1 - eps) // 2
        lam = (1 - legendre(p, q)) // 2
        M4 = [[beta, lam + eta], [beta, lam]]
        M8 = [[beta, beta + lam + eta], [beta, lam]]
        rec["X3_matrix_minus4N"] = G4N["r4"] == 2 - rank_f2(M4)
        rec["X3_matrix_minus8N"] = G8N["r4"] == 2 - rank_f2(M8)
        if eps == 1:
            rec["X4_milovic"] = (G4N["r4"] == 2) == (p % 8 == 1 and q % 8 == 1 and legendre(p, q) == 1)
    return rec


def run(bound: int = 300, class_groups: bool = True) -> Dict[str, object]:
    t0 = time.time()
    P = odd_primes_below(bound)
    pairs = [(p, q) for i, p in enumerate(P) for q in P[i + 1:]]
    in_scope = [(p, q) for p, q in pairs if p % 4 == q % 4]
    records = [check_pair(p, q, class_groups) for p, q in in_scope]
    keys = sorted({k for r in records for k, v in r.items() if isinstance(v, bool)})
    summary = {}
    for k in keys:
        vals = [r[k] for r in records if k in r]
        fails = [(r["p"], r["q"]) for r in records if k in r and not r[k]]
        summary[k] = {"checked": len(vals), "hold": sum(vals), "fail": len(vals) - sum(vals), "first_failures": fails[:10]}
    return {
        "bound": bound,
        "odd_primes": len(P),
        "pairs_total": len(pairs),
        "pairs_p_eq_q_mod4": len(in_scope),
        "pairs_N_1_mod_8": sum(1 for r in records if r["N_mod_8"] == 1),
        "pairs_N_5_mod_8": sum(1 for r in records if r["N_mod_8"] == 5),
        "summary": summary,
        "max_h_minus8N": max((r.get("h_minus8N", 0) for r in records), default=0),
        "seconds": round(time.time() - t0, 2),
        "records": records,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bound", type=int, default=300)
    ap.add_argument("--out", default=None, help="write the full result (with per-pair records) as JSON")
    args = ap.parse_args()
    res = run(args.bound)
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(res, fh, sort_keys=True)
            fh.write("\n")
    print(json.dumps({k: v for k, v in res.items() if k != "records"}, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
