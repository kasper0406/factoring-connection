"""E1 (Lemma S1 exhaustive/random check) and E3 (Theorem SP per-prime dichotomy check).

Usage: python -m src.slpcirc.e1_e3 E1|E3 OUT.jsonl
"""
from __future__ import annotations
import json, random, sys, time
import numpy as np
import flint
from sympy import primerange, factorint
from .sparse import primitive_root, s1_coset_bound, Lambda, dichotomy_violations


def zeros_table(p, exps, coeff_rows):
    """For fixed exponents, count zeros on F_p^* for each coefficient row (vectorized)."""
    n = p - 1
    g = primitive_root(p)
    pw = np.array([pow(g, j, p) for j in range(n)], dtype=np.int64)  # g^j
    # value of u^e at u = g^j is g^(e j mod n)
    idx = [((e * np.arange(n)) % n) for e in exps]
    mons = np.stack([pw[ix] for ix in idx])            # t x n
    C = np.array(coeff_rows, dtype=np.int64)            # R x t
    vals = (C @ mons) % p                               # R x n  (entries < t*p^2, safe)
    return (vals == 0).sum(axis=1)


def check_bound(p, exps, z):
    n = p - 1
    viol = []
    for ell in factorint(n):
        if len({e % ell for e in exps}) == len(exps):
            if z > (n // ell) * s1_coset_bound(ell, len(exps)):
                viol.append(ell)
    return viol


def run_E1(out):
    t0 = time.time()
    checked = 0; violations = []
    maxratio = 0.0
    for p in primerange(3, 44):
        n = p - 1
        # t = 2: exponents {0, a}
        for a in range(1, n):
            rows = [[1, c] for c in range(1, p)]
            zs = zeros_table(p, [0, a], rows)
            for c, z in zip(range(1, p), zs):
                checked += 1
                v = check_bound(p, [0, a], int(z))
                if v:
                    violations.append([p, [0, a], [1, c], int(z), v])
        # t = 3
        rows = [[1, c2, c3] for c2 in range(1, p) for c3 in range(1, p)]
        for a in range(1, n):
            for b in range(a + 1, n):
                zs = zeros_table(p, [0, a, b], rows)
                for row, z in zip(rows, zs):
                    checked += 1
                    v = check_bound(p, [0, a, b], int(z))
                    if v:
                        violations.append([p, [0, a, b], row, int(z), v])
    rng = random.Random(20261020)
    primes = list(primerange(3, 212))
    for _ in range(200000):
        p = rng.choice(primes); n = p - 1
        t = rng.choice([4, 5])
        if n < t:
            continue
        exps = [0] + rng.sample(range(1, n), t - 1)
        row = [1] + [rng.randrange(1, p) for _ in range(t - 1)]
        z = int(zeros_table(p, exps, [row])[0])
        checked += 1
        v = check_bound(p, exps, z)
        if v:
            violations.append([p, exps, row, z, v])
    res = {"exp": "E1", "checked": checked, "violations": len(violations), "examples": violations[:20],
           "seconds": round(time.time() - t0, 1)}
    with open(out, "a") as fh:
        fh.write(json.dumps(res) + "\n")
    print(json.dumps(res))


def roots_count_sparse(coeffs_exps, p):
    """#roots in F_p of sum c x^e (exponents reduced mod p-1 for x != 0; x = 0 checked directly)."""
    n = p - 1
    red = {}
    for c, e in coeffs_exps:
        k = e % n if e > 0 else 0
        # x^e for x != 0 equals x^(e mod n); e = 0 gives 1
        red[k] = (red.get(k, 0) + c) % p
    poly = [0] * (max(red) + 1 if red else 1)
    for k, c in red.items():
        poly[k] = c
    f = flint.nmod_poly(poly, p)
    zero_at_0 = sum(c for c, e in coeffs_exps if e == 0) % p == 0
    if f.is_zero():
        cnt_star = n
    else:
        xp = flint.nmod_poly([0, 1], p).pow_mod(n, f) if f.degree() > 0 else None
        if f.degree() <= 0:
            cnt_star = 0
        else:
            h = f.gcd(xp - 1)
            cnt_star = h.degree()
    return cnt_star + (1 if zero_at_0 else 0)


def lcm_upto(B):
    from math import lcm
    L = 1
    for k in range(2, B + 1):
        L = lcm(L, k)
    return L


_PW = {}


def roots_count_sparse_np(coeffs_exps, p):
    """Same as roots_count_sparse, by exact vectorized evaluation at all g^j (fast for p < 2^20)."""
    n = p - 1
    if p not in _PW:
        g = primitive_root(p)
        pw = np.empty(n, dtype=np.int64); v = 1
        for j in range(n):
            pw[j] = v; v = v * g % p
        _PW.clear(); _PW[p] = pw
    pw = _PW[p]
    j = np.arange(n, dtype=np.int64)
    acc = np.zeros(n, dtype=np.int64)
    for c, e in coeffs_exps:
        k = e % n if e > 0 else 0
        acc += (c % p) * pw[(k * j) % n]
    cnt_star = int(((acc % p) == 0).sum())
    zero_at_0 = sum(c for c, e in coeffs_exps if e == 0) % p == 0
    return cnt_star + (1 if zero_at_0 else 0)


def run_E3(out):
    t0 = time.time()
    X = 2 ** 14
    primes = list(primerange(X + 1, 2 * X + 1))
    fams = []
    for B in (50, 200, 1000):
        M = lcm_upto(B)
        fams.append((f"x^lcm({B})-1", [(1, M), (-1, 0)]))
    rng = random.Random(1729)
    for i in range(200):
        t = 3 if i < 100 else 4
        exps = rng.sample(range(1, 2 ** 40), t - 1) + [0]
        cs = [rng.choice([c for c in range(-5, 6) if c]) for _ in range(t)]
        fams.append((f"rand{t}_{i}", list(zip(cs, exps))))
    M1, M2 = lcm_upto(300), lcm_upto(120) * 7 ** 3
    fams.append(("x^M+x^M'-2", [(1, M1), (1, M2), (-2, 0)]))
    lam_num, lam_den = 1, 20
    viol = []; summary = []
    hits_by = {name: 0 for name, _ in fams}
    for p in primes:
        for name, ce in fams:
            r = roots_count_sparse_np(ce, p)
            if 20 * r >= p:
                hits_by[name] += 1
                v = dichotomy_violations(ce, p, lam_num, lam_den, r)
                if v:
                    viol.append([name, p, r, v])
    summary = [[name, hits_by[name]] for name, _ in fams]
    res = {"exp": "E3", "X": X, "nprimes": len(primes), "lambda": "1/20", "violations": len(viol),
           "viol_examples": viol[:20], "hits_per_family": summary, "seconds": round(time.time() - t0, 1)}
    with open(out, "a") as fh:
        fh.write(json.dumps(res) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "hits_per_family"}))
    print([s for s in summary if s[1]])


if __name__ == "__main__":
    {"E1": run_E1, "E3": run_E3}[sys.argv[1]](sys.argv[2])
