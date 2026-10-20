"""E2 (notes/dlog25_prereg.md): Lemma XK4-U (unit elimination) at degrees d = 3, 5, 7 with the EXPLICIT
precision kappa and relation bound R of notes/dlog25_results.md §4 (no ad-hoc precision).

Setup (toy, disclosed): monic irreducible f of odd degree d with small coefficients; pair elements a - b*alpha
(gcd(a, b) = 1, b >= 1) whose norm F(a, b) is B-smooth with no prime factor dividing disc(f) (good primes only:
then the prime ideals are (r, a/b mod r) with valuation v_r(F)).  Exact integer kernel of the valuation matrix
=> units gamma_x (Lemma XK1).  Take s = r_U + 1 kernel vectors, compute certified Log vectors (arb), run LLL on
b_j = (e_j, round(2^kappa v_j)) with kappa from the explicit bound, certify the first vector by the
Blanksby-Montgomery gap, and check torsion (gamma = +-1, d odd) at 5 auxiliary degree-one primes (exact).

DIAGNOSTIC code: no secrets exist in this experiment.
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from math import gcd

import flint
from sympy import primerange, isprime


def hom_eval(f, a, b):
    d = len(f) - 1
    return sum(c * a ** i * b ** (d - i) for i, c in enumerate(f))


def random_field(d, rng, coef=3):
    while True:
        f = [rng.randint(-coef, coef) for _ in range(d)] + [1]
        if f[0] == 0:
            continue
        P = flint.fmpz_poly(f)
        fac = P.factor()[1]
        if len(fac) == 1 and fac[0][1] == 1:
            disc = int(flint.fmpz_poly(f).discriminant()) if hasattr(P, "discriminant") else None
            return f, disc


def disc_of(f):
    import sympy as sp
    x = sp.symbols("x")
    return int(sp.discriminant(sum(c * x ** i for i, c in enumerate(f)), x))


def smooth_factor(n, B):
    n = abs(n)
    fac = flint.fmpz(n).factor()
    out = {}
    for q, e in fac:
        q = int(q)
        if q > B:
            return None
        out[q] = int(e)
    return out


def bm_gap(d):
    x = flint.arb(1) / (52 * d * flint.arb(6 * d).log())
    return (1 + x).log() / d


def collect(f, disc, B, box_a, box_b, cap):
    rels = []
    for b in range(1, box_b + 1):
        for a in range(-box_a, box_a + 1):
            if gcd(a, b) != 1:
                continue
            F = hom_eval(f, a, b)
            if F == 0:
                continue
            fac = smooth_factor(F, B)
            if fac is None or any(disc % r == 0 for r in fac):
                continue
            cols = {(r, a * pow(b, -1, r) % r): e for r, e in fac.items()}
            rels.append((a, b, cols))
            if len(rels) >= cap:
                return rels
    return rels


def exact_kernel(rels):
    keys = sorted({k for _, _, c in rels for k in c})
    idx = {k: j for j, k in enumerate(keys)}
    M = [[0] * len(rels) for _ in keys]
    for i, (_, _, c) in enumerate(rels):
        for k, v in c.items():
            M[idx[k]][i] = v
    if not keys:
        return len(keys), [[1 if i == j else 0 for i in range(len(rels))] for j in range(len(rels))]
    K, nul = flint.fmpz_mat(M).nullspace()
    basis = [[int(K[i, j]) for i in range(len(rels))] for j in range(nul)]
    if basis:
        Lb = flint.fmpz_mat(basis).lll()          # small entries (harmless: any basis of Lambda_1 works)
        basis = [[int(Lb[i, j]) for j in range(len(rels))] for i in range(Lb.nrows())]
    return len(keys), basis


def embeddings(f, prec):
    flint.ctx.prec = prec
    roots = flint.fmpz_poly(f).complex_roots()
    real = [r for r, m in roots if r.imag.is_zero()]
    cplx = [r for r, m in roots if not r.imag.is_zero() and r.imag > 0]
    assert len(real) + 2 * len(cplx) == len(f) - 1
    return real, cplx


def log_vec(rels, x, real, cplx):
    out = [flint.arb(0)] * (len(real) + len(cplx))
    for xi, (a, b, _) in zip(x, rels):
        if xi == 0:
            continue
        lv = [abs(flint.acb(a) - b * r).log() for r in real] + [2 * abs(flint.acb(a) - b * r).log() for r in cplx]
        out = [o + xi * l for o, l in zip(out, lv)]
    return out


def explicit_kappa(d, s, r1r2, rU, vmax2, beta_d):
    """kappa from §4 (Lemma XK4-U): C1 = (3/2)^rU 2^{rU/2} vmax/beta_d; R = max(1, rU! C1^rU);
    R' = 2^{(s-1)/2}(1 + s sqrt(r1r2)) R; need 2^kappa beta_d > (1 + s sqrt(r1r2)) R'."""
    l2 = math.log2
    C1 = rU * l2(1.5) + rU / 2 + l2(max(vmax2, 1e-300)) - l2(beta_d)
    logR = max(0.0, l2(math.factorial(rU)) + rU * C1) if rU > 0 else 0.0
    logRp = (s - 1) / 2 + l2(1 + s * math.sqrt(r1r2)) + logR
    kappa = math.ceil(l2(1 + s * math.sqrt(r1r2)) + logRp - l2(beta_d)) + 2
    return kappa, logR


def aux_primes(f, disc, rng, count=5, bits=40):
    out = []
    while len(out) < count:
        q = rng.randrange(2 ** (bits - 1), 2 ** bits)
        if not isprime(q) or disc % q == 0:
            continue
        roots = flint.nmod_poly(f, q).roots()
        if roots:
            out.append((q, int(roots[0][0])))
    return out


def psi(rels, x, q, s_):
    v = 1
    for xi, (a, b, _) in zip(x, rels):
        if xi == 0:
            continue
        z = (a - b * s_) % q
        v = v * pow(z, xi, q) % q
    return v


def run_one(d, seed, B, box_a, box_b, cap, coef=3):
    rng = random.Random(seed)
    f, _ = random_field(d, rng, coef)
    disc = disc_of(f)
    t0 = time.time()
    rels = collect(f, disc, B, box_a, box_b, cap)
    ncols, basis = exact_kernel(rels)
    real, cplx = embeddings(f, 64)
    r1, r2 = len(real), len(cplx)
    rU = r1 + r2 - 1
    s = rU + 1
    rec = {"d": d, "seed": seed, "f": f, "r1": r1, "r2": r2, "rU": rU, "rels": len(rels), "cols": ncols,
           "kernel_dim": len(basis)}
    # kernel vectors with nontrivial gamma: all nonzero x give gamma_x != 1 as an element unless x = 0 (distinct
    # coprime pairs are multiplicatively independent only generically) -- we just use the first s vectors.
    if len(basis) < s:
        rec["status"] = "too_few_kernel_vectors"
        rec["secs"] = round(time.time() - t0, 2)
        return rec
    X = basis[:s]
    bx = max(abs(v).bit_length() for x in X for v in x)
    # low-precision pass for vmax (upper bound used only in kappa)
    lv0 = [log_vec(rels, x, *embeddings(f, 128)) for x in X]
    vmax2 = max(math.sqrt(sum(float(c.mid()) ** 2 for c in v)) for v in lv0) * 1.01 + 1e-9
    beta_d = float(bm_gap(d))
    kappa, logR = explicit_kappa(d, s, r1 + r2, rU, vmax2, beta_d)
    P = kappa + bx + max(1, len(rels)).bit_length() + 40 + d * 8   # precision ledger (+ slack for root error)
    real, cplx = embeddings(f, P + 64)
    V = [log_vec(rels, x, real, cplx) for x in X]
    nonunit = sum(1 for v in V if not sum(v, flint.arb(0)).contains(0))
    scale = flint.arb(2) ** kappa
    rows = []
    for j, v in enumerate(V):
        r = [0] * s
        r[j] = 1
        r += [int((c * scale).mid().floor().unique_fmpz()) for c in v]
        rows.append(r)
    L = flint.fmpz_mat(rows).lll()
    n = [int(L[0, j]) for j in range(s)]
    comb = [flint.arb(0)] * (r1 + r2)
    for nj, v in zip(n, V):
        comb = [c + nj * vv for c, vv in zip(comb, v)]
    gap = flint.arb(beta_d) / 2
    certified = any(n) and all(abs(c) < gap for c in comb)
    xstar = [sum(nj * x[i] for nj, x in zip(n, X)) for i in range(len(rels))]
    g = 0
    for v in xstar:
        g = gcd(g, v)
    if g:
        xstar = [v // g for v in xstar]
    aux = aux_primes(f, disc, rng)
    vals = [psi(rels, xstar, q, s_) for q, s_ in aux]
    pm1 = [v in (1, q - 1) for v, (q, _) in zip(vals, aux)]
    signs = {1 if v == 1 else -1 for v, (q, _) in zip(vals, aux) if v in (1, q - 1)}
    # all LLL rows: count certified relations (implied unit rank) and check non-relations fail the torsion test
    nrel, nonrel_pass = 0, 0
    for i in range(L.nrows()):
        ni = [int(L[i, j]) for j in range(s)]
        if not any(ni):
            continue
        ci = [flint.arb(0)] * (r1 + r2)
        for nj, v in zip(ni, V):
            ci = [c + nj * vv for c, vv in zip(ci, v)]
        xi = [sum(nj * x[t] for nj, x in zip(ni, X)) for t in range(len(rels))]
        if all(abs(c) < gap for c in ci):
            nrel += 1
        elif any(abs(c) > gap for c in ci):
            q, s_ = aux[0]
            if psi(rels, xi, q, s_) in (1, q - 1):
                nonrel_pass += 1
    rec.update({"status": "ok", "beta_x_bits": bx, "beta_d": beta_d, "kappa": kappa, "log2R": round(logR, 2),
                "log2_n": round(math.log2(max(1, math.sqrt(sum(v * v for v in n)))), 2), "n_nonzero": any(n),
                "certified": certified, "aux_pm1": all(pm1), "aux_sign_consistent": len(signs) == 1,
                "nonunit_kernel_vectors": nonunit, "lll_relations": nrel, "unit_rank_implied": s - nrel,
                "nonrelation_rows_passing_aux": nonrel_pass, "secs": round(time.time() - t0, 2)})
    return rec


BLOCKS = [
    # d, seeds, B, box_a, box_b, cap, coef
    (3, range(0, 8), 600, 150, 150, 400, 3),
    (5, range(10, 18), 1500, 150, 150, 600, 3),
    (7, range(20, 24), 5000, 300, 300, 1500, 1),     # descriptive only (calibration: kernel too small)
]


def main(out_path):
    with open(out_path, "w") as fh:
        for d, seeds, B, ba, bb, cap, coef in BLOCKS:
            for sd in seeds:
                rec = run_one(d, sd, B, ba, bb, cap, coef)
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                print(rec, flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
