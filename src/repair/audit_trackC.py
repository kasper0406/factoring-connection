"""Reproductions of the r3 audit's Track C items (agent `repair`).  DIAGNOSTIC.

Run:  .venv/bin/python -m src.repair.audit_trackC --out experiments/raw/repair_audit_trackC.jsonl
"""
from __future__ import annotations

import argparse
import json
from math import comb, gcd, isqrt, lcm
from typing import Dict, List

from src.repair.c2r_search import is_prime


def factor_small(n: int) -> Dict[int, int]:
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


# ---------------------------------------------------------------------------
# O1 polynomial-readout corollary: the bound is on the lifted dimension
# ---------------------------------------------------------------------------
def o1_poly_readout(primes=(7, 11, 13, 31, 101, 257)) -> Dict:
    rows = []
    for p in primes:
        nr = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)  # s^2 = nr
        mul = lambda u, v: ((u[0] * v[0] + nr * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p)

        def pw(u, e):
            r = (1, 0)
            while e:
                if e & 1:
                    r = mul(r, u)
                u = mul(u, u)
                e >>= 1
            return r
        order = p * p - 1
        primes_div = list(factor_small(order))
        alpha = next((a, b) for a in range(p) for b in range(1, p)
                     if all(pw((a, b), order // l) != (1, 0) for l in primes_div))
        norm = lambda u: (u[0] * u[0] - nr * u[1] * u[1]) % p  # quadratic form in coordinates
        t1, t2 = pw(alpha, p + 1), pw(alpha, p + 2)
        u = []
        z = (1, 0)
        for K in range(1, p + 3):
            z = mul(z, alpha)
            u.append(norm(((z[0] - t1[0]) % p, (z[1] - t1[1]) % p)) * norm(((z[0] - t2[0]) % p, (z[1] - t2[1]) % p)) % p)
        rows.append({"p": p, "state_dim": 2, "readout_degree": 4, "lifted_dim": comb(2 + 4, 4),
                     "nonzero_1..p": all(u[K - 1] for K in range(1, p + 1)),
                     "zero_at_p+1,p+2": u[p] == 0 and u[p + 1] == 0})
    # interval version: state (A, 1) (dimension 2), readout prod_{i<L}(A+i) of degree L
    p, L = 101, 7
    ZL = {(-i) % p for i in range(L)}
    ok = all((eval_poch(A, L, p) == 0) == (A % p in ZL) for A in range(1, 2 * p))
    return {"item": "C-O1-poly", "prefix_examples": rows,
            "interval_example": {"p": p, "L": L, "state_dim": 2, "readout_degree": L,
                                 "zero_iff_A_mod_p_in_Z_L_for_1<=A<2p": ok}}


def eval_poch(A: int, L: int, p: int) -> int:
    v = 1
    for i in range(L):
        v = v * (A + i) % p
    return v


# ---------------------------------------------------------------------------
# O2 strengthening: deg Q >= ceil(L/2)
# ---------------------------------------------------------------------------
def o2_strengthening(pmax: int = 100) -> Dict:
    viol = 0
    gviol = 0
    configs = 0
    attained = {}
    for p in range(5, pmax):
        if not is_prime(p):
            continue
        for e in range(1, p):
            pw = [pow(y, e, p) for y in range(p)]
            fib = {}
            for y in range(p):
                fib[pw[y]] = fib.get(pw[y], 0) + 1
            g = gcd(e, p - 1)
            for u in range(p):
                img = set()
                tot = 0
                for L in range(1, p - 1):
                    y = (u - (L - 1)) % p  # T = Z_L + u grows by -(L-1)+u
                    if pw[y] not in img:
                        img.add(pw[y])
                        tot += fib[pw[y]]
                    if L >= 2 and tot == L:  # T is a union of fibres => exact one-layer fold
                        configs += 1
                        viol += len(img) < -(-L // 2)
                        gviol += g > 2
                        key = (p, L)
                        attained[key] = min(attained.get(key, 10 ** 9), len(img))
    att = sum(1 for (p, L), d in attained.items() if d == -(-L // 2))
    return {"item": "C-O2", "pmax": pmax, "exact_configs": configs, "violations_ceil_L_over_2": viol,
            "violations_g<=2": gviol, "(p,L)_with_exact_fold": len(attained),
            "(p,L)_attaining_ceil_L_over_2": att}


def o2prime_consequences() -> Dict:
    p, B = 257, 10
    g = gcd(lcm(*range(1, B + 1)), p - 1)
    smooth = 1
    for l, a in factor_small(p - 1).items():
        if l <= B:
            smooth *= l ** a
    differ = 0
    n = 0
    for q in range(3, 20000):
        if not is_prime(q):
            continue
        n += 1
        f = factor_small(q - 1)
        sm = 1
        for l, a in f.items():
            if l <= B:
                sm *= l ** a
        differ += gcd(lcm(*range(1, B + 1)), q - 1) != sm
    L = 1000
    k = (L - 1).bit_length()
    return {"item": "C-O2prime", "p": p, "B": B, "gcd(lcm(1..B),p-1)": g, "B_smooth_part": smooth,
            "primes<20000_where_they_differ": differ, "of": n,
            "bounded_gcd_example": {"G": 2, "L": L, "k": k, "bound_L/G^k": L / 2 ** k}}


def o3_degree() -> Dict:
    p, L = 101, 10
    zero = {A for A in range(p) if eval_poch(A, L, p) == 0}
    return {"item": "C-O3", "p": p, "L": L, "degree": L, "degree<p": L < p,
            "zero_set_is_Z_L": zero == {(-i) % p for i in range(L)}}


# ---------------------------------------------------------------------------
# O4: annihilation needs the element's order, not the group order; and G(p) vs sqrt(p)
# ---------------------------------------------------------------------------
def ec_add(P, Q, a, p):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def ec_mul(k, P, a, p):
    R = None
    while k:
        if k & 1:
            R = ec_add(R, P, a, p)
        P = ec_add(P, P, a, p)
        k >>= 1
    return R


def o4_order_example(p: int = 101, a: int = 2, b: int = 3) -> Dict:
    pts = [None] + [(x, y) for x in range(p) for y in range(p) if (y * y - x ** 3 - a * x - b) % p == 0]
    nE = len(pts)
    best = None
    for P in pts[1:]:
        n = next(d for d in range(1, nE + 1) if nE % d == 0 and ec_mul(d, P, a, p) is None)
        if n < nE:
            best = (P, n)
            break
    P, n = best
    return {"item": "C-O4-order", "curve": f"y^2=x^3+{a}x+{b} over F_{p}", "#E": nE, "point": P,
            "order": n, "n*P=O": ec_mul(n, P, a, p) is None, "#E divides n": n % nE == 0}


def generation_bound(pmax: int = 20000) -> Dict:
    """G(p) = least y with <1..y> = F_p^*.  O4's lemma gives G(p) <= sqrt(p) + 1; known
    results (Burthe; Martin-Parvardi Cor. 3) give G(p) << p^(1/(4 sqrt e) + eps)."""
    worst_ratio = (0.0, None)
    worst = (0, None)
    for p in range(5, pmax):
        if not is_prime(p):
            continue
        f = list(factor_small(p - 1))
        ordl = 1
        y = 1
        while ordl != p - 1:
            y += 1
            if not is_prime(y):
                continue
            o = p - 1
            for l in f:
                while o % l == 0 and pow(y, o // l, p) == 1:
                    o //= l
            ordl = lcm(ordl, o)
        r = y / (isqrt(p) + 1)
        if r > worst_ratio[0]:
            worst_ratio = (r, p)
        if y > worst[0]:
            worst = (y, p)
    return {"item": "C-O4-G(p)", "pmax": pmax, "max_G(p)": worst, "max_G(p)/(sqrt p + 1)": worst_ratio}


def m2_degree(p: int = 7, e: int = 10, L: int = 3) -> Dict:
    direct = [sum(pow(A + i, e, p) for i in range(L)) % p for A in range(p)]
    expanded = [sum(comb(e, k) * pow(A, k, p) * sum(pow(i, e - k, p) if e - k else 1 for i in range(L))
                    for k in range(e + 1)) % p for A in range(p)]
    return {"item": "C-M2", "p": p, "e": e, "e>=p-1": e >= p - 1, "L": L,
            "direct_expansion_polynomial_of_degree<=e_matches": direct == expanded}


ALL = [o1_poly_readout, o2_strengthening, o2prime_consequences, o3_degree, o4_order_example,
       generation_bound, m2_degree]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/raw/repair_audit_trackC.jsonl")
    a = ap.parse_args(argv)
    with open(a.out, "w") as fh:
        for f in ALL:
            r = f()
            fh.write(json.dumps(r, default=str) + "\n")
            print(json.dumps(r, default=str))


if __name__ == "__main__":
    main()
