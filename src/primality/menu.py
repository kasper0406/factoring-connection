"""HARNESS / DIAGNOSTIC: the exponent x group menu ("are there more favourable checks than a^(N-1)?").

For a uniform element of a finite abelian group G(F_r) = C_(d_1) x ... x C_(d_k) (d_i | d_(i+1)),
    A_r(E) := P(x^E = 1) = |G[E]| / |G| = prod_i gcd(E, d_i) / |G|            (Lemma M1)
and, for N = pq with independent CRT components, the gcd(x^E - 1, N) leak probability is
    A_p + A_q - 2 A_p A_q.
Groups: (Z/N)^* ; the norm-one torus T (Selfridge D) ; the full unit group of Z/N[sqrt D]
(inert: cyclic p^2-1, split: C_(p-1)^2) ; the idealised cyclotomic tori T_k (order Phi_k(r)) and
F_(r^k)^* (as if a degree-k field were available at both primes; Bach-Shallit) ; random elliptic
curves (32-bit N only, cyclic approximation gcd(E,#E)/#E, exact formula checked by enumeration for r < 200).
Exponents: N-1, N+1, N, (N-1)/2, N^2-1, N^k-1 (k = 3, 4, 6), L_B = lcm(1..B), and products.

Run: .venv/bin/python -m src.primality.experiments menu
"""
from __future__ import annotations

import json
import math
import random
from fractions import Fraction as Fr
from math import gcd
from multiprocessing import Pool

import numpy as np
import sympy

from src.common.instances import balanced_semiprime
from src.common.metrics import append_jsonl
from src.primality import ops

RAW = "experiments/raw/"


def lcm_upto(B: int) -> int:
    L = 1
    for p in sympy.primerange(2, B + 1):
        pk = p
        while pk * p <= B:
            pk *= p
        L *= pk
    return L


LB = {B: lcm_upto(B) for B in (100, 1000, 10000, 100000)}


def exponents(N: int):
    E = {"N-1": N - 1, "N+1": N + 1, "N": N, "(N-1)/2": (N - 1) // 2, "N^2-1": N * N - 1,
         "N^3-1": N ** 3 - 1, "N^4-1": N ** 4 - 1, "N^6-1": N ** 6 - 1}
    for B, L in LB.items():
        E[f"L{B}"] = L
    E["(N-1)L1000"] = (N - 1) * LB[1000]
    E["(N^2-1)L1000"] = (N * N - 1) * LB[1000]
    E["(N^6-1)L1000"] = (N ** 6 - 1) * LB[1000]
    E["(N^2-1)(N^6-1)L10000"] = (N * N - 1) * (N ** 6 - 1) * LB[10000]
    return E


def groups(r: int, D: int):
    """Invariant factors of G(F_r) for each group family (D = public Selfridge discriminant)."""
    chi = int(sympy.legendre_symbol(D % r, r))
    G = {"Zn*": [r - 1], "torus_D": [r - chi],
         "unit_Z[sqrtD]": [r * r - 1] if chi == -1 else [r - 1, r - 1]}
    for k in (2, 3, 4, 6):
        G[f"T{k}_ideal"] = [int(sympy.cyclotomic_poly(k, r))]
        G[f"F{k}_ideal"] = [r ** k - 1]
    return G


def A(E: int, inv) -> Fr:
    num = 1
    for d in inv:
        num *= gcd(E, d)
    den = 1
    for d in inv:
        den *= d
    return Fr(num, den)


def ec_order(r: int, a: int, b: int, qr: np.ndarray) -> int:
    x = np.arange(r, dtype=np.int64)
    f = (x * x % r * x + a * x + b) % r
    # number of y with y^2 = f: 1 + legendre(f)
    leg = np.where(f == 0, 0, np.where(qr[f], 1, -1))
    return int(r + 1 + leg.sum())


def menu_instance(args):
    bits, seed, with_ec = args
    inst = balanced_semiprime(bits, seed)
    N = inst.N
    p, q = inst.labels["factors"]
    sel = ops.selfridge(N)
    D = sel.get("D", 5)
    Es = exponents(N)
    Gp, Gq = groups(p, D), groups(q, D)
    rec = {"bits": bits, "seed": seed, "D": D, "leak": {}, "Ap": {}, "Aq": {}}
    for gname in Gp:
        for ename, E in Es.items():
            ap, aq = A(E, Gp[gname]), A(E, Gq[gname])
            lk = ap + aq - 2 * ap * aq
            key = f"{gname}|{ename}"
            rec["leak"][key] = float(lk) if lk > Fr(1, 10 ** 300) else 0.0
            rec["Ap"][key] = math.log2(ap.numerator) - math.log2(ap.denominator)
    if with_ec:
        rng = random.Random(f"ec|{bits}|{seed}")
        qrs = {}
        for r in (p, q):
            t = np.zeros(r, dtype=bool)
            t[(np.arange(1, r, dtype=np.int64) ** 2) % r] = True
            qrs[r] = t
        ec = []
        for _ in range(3):
            a_, b_ = rng.randrange(N), rng.randrange(N)
            np_, nq_ = ec_order(p, a_ % p, b_ % p, qrs[p]), ec_order(q, a_ % q, b_ % q, qrs[q])
            row = {}
            for ename in ("N-1", "N+1", "N", "L100", "L1000", "L10000", "(N+1)L1000"):
                E = Es[ename] if ename in Es else (N + 1) * LB[1000]
                ap, aq = Fr(gcd(E, np_), np_), Fr(gcd(E, nq_), nq_)
                row[ename] = float(ap + aq - 2 * ap * aq)
            ec.append(row)
        rec["ec"] = ec
    return rec


def check_ec_formula(seed=0):
    """Exact check of |E[e]| = gcd(e,d1) gcd(e,d2) by enumerating all points, r < 200."""
    rng = random.Random(seed)
    bad = checked = 0
    for r in sympy.primerange(5, 200):
        for _ in range(4):
            a, b = rng.randrange(r), rng.randrange(r)
            if (4 * a ** 3 + 27 * b * b) % r == 0:
                continue
            pts = [None] + [(x, y) for x in range(r) for y in range(r) if (y * y - x ** 3 - a * x - b) % r == 0]

            def add(P, Q):
                if P is None:
                    return Q
                if Q is None:
                    return P
                if P[0] == Q[0] and (P[1] + Q[1]) % r == 0:
                    return None
                if P == Q:
                    l = (3 * P[0] ** 2 + a) * pow(2 * P[1], -1, r) % r
                else:
                    l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, r) % r
                x3 = (l * l - P[0] - Q[0]) % r
                return (x3, (l * (P[0] - x3) - P[1]) % r)

            def mul(e, P):
                R = None
                while e:
                    if e & 1:
                        R = add(R, P)
                    P = add(P, P)
                    e >>= 1
                return R

            n = len(pts)
            orders = []
            for P in pts:
                o = min(d for d in sympy.divisors(n) if mul(d, P) is None)
                orders.append(o)
            d2 = max(orders)
            d1 = n // d2
            for e in rng.sample(range(1, 4 * n), 6) + [n, d2, 2, 3, 12]:
                cnt = sum(1 for P in pts if mul(e, P) is None)
                checked += 1
                bad += cnt != gcd(e, d1) * gcd(e, d2)
    return {"checked": checked, "bad": bad}


def check_identities(bound=2000):
    """gcd(E(N), Psi(p)) = gcd(E reduced via p's relation, Psi(p)) for N = pq: exhaustive over pairs."""
    bad = cnt = 0
    for p in sympy.primerange(3, 120):
        for q in sympy.primerange(p + 2, bound // p + 1):
            if q == p:
                continue
            N = p * q
            for k in (1, 2, 3, 4, 6):
                cnt += 1
                bad += gcd(N ** k - 1, p ** k - 1) != gcd(q ** k - 1, p ** k - 1)
                ph = int(sympy.cyclotomic_poly(k, p))
                bad += gcd(N ** k - 1, ph) != gcd(q ** k - 1, ph)
            bad += gcd(N + 1, p + 1) != gcd(q - 1, p + 1)
            bad += gcd(N - 1, p + 1) != gcd(q + 1, p + 1)
            bad += gcd(N, p - 1) != gcd(q, p - 1)
            bad += gcd(N * N - 1, p - 1) != gcd(q * q - 1, p - 1)
    return {"pairs_x_k": cnt, "bad": bad}


def run(a):
    cfg = {"bits": {32: 400, 48: 300, 64: 300}, "ec_bits": 32}
    json.dump(cfg, open("experiments/configs/primality_menu.json", "w"), indent=1)
    out = RAW + "primality_menu.jsonl"
    open(out, "w").close()
    jobs = [(b, s, b == cfg["ec_bits"]) for b, n in cfg["bits"].items() for s in range(n)]
    with Pool(a.workers) as pool:
        for rec in pool.imap_unordered(menu_instance, jobs, chunksize=4):
            append_jsonl(out, rec)
    chk = {"ec_formula": check_ec_formula(), "identities": check_identities()}
    print(chk)
    json.dump(chk, open(RAW + "primality_menu_checks.json", "w"), indent=1)
