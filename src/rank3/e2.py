"""Experiment E2 (pre-registered, notes/rank3_prereg.md): toy-size construction comparison, rank 2 vs rank 3.

Measured only; NOT a test of asymptotics. All norms and smoothness tests are exact integers. The LLL weight search
for design R3-lll uses rounded real embeddings (a heuristic search step, disclosed); the counts are exact.
Harness code: N = pq is generated here (primality test used), and the design code receives only N and m-free data
derived from N (base-m / lattice polynomials); p, q are never used except to label output.

Run: .venv/bin/python -m src.rank3.e2      (writes experiments/raw/rank3_E2.jsonl)
"""
from __future__ import annotations

import itertools
import json
import math
import random
import statistics
import sys
import time
from typing import Dict, List, Sequence

import flint
import gmpy2
import sympy
from flint import acb, arb, ctx, fmpz_mat, fmpz_poly

from src.rank3.forms import fpoly, poly_comb

SEED = 20261020
OUT = "experiments/raw/rank3_E2.jsonl"
BSM = 1 << 12
PRIMORIAL = gmpy2.mpz(math.prod(int(p) for p in sympy.primerange(2, BSM + 1)))
E2D, E3D = 256, 32  # rank-2 box |a| <= 256, 1 <= b <= 256; rank-3 box |a|,|b| <= 32, 1 <= c <= 32


def iroot(n: int, k: int) -> int:
    x = int(gmpy2.iroot(gmpy2.mpz(n), k)[0])
    return x


def smooth(n: int) -> bool:
    n = gmpy2.mpz(abs(n))
    if n == 0:
        return False
    while n > 1:
        g = gmpy2.gcd(n, PRIMORIAL)
        if g == 1:
            return False
        n //= g
    return True


def irreducible(c) -> bool:
    fac = fmpz_poly([int(x) for x in c]).factor()[1]
    return len(fac) == 1 and fac[0][1] == 1


def make_N(rnd: random.Random, n: int):
    while True:
        p = int(sympy.nextprime(rnd.randrange(1 << (n // 2 - 1), 1 << (n // 2))))
        q = int(sympy.nextprime(rnd.randrange(1 << (n - n // 2 - 1), 1 << (n - n // 2))))
        if p != q and (p * q).bit_length() == n:
            return p * q, p, q


def base_m(N: int, d: int):
    m = iroot(N, d + 1)
    c, r = [], N
    for _ in range(d):
        c.append(r % m)
        r //= m
    c.append(r)
    return c, [-m, 1], m


def lattice_pair(N: int, m: int, d: int = 3):
    """Two short irreducible coprime polynomials of degree >= 2 in {h in Z[x]_{<=d} : h(m) = 0 mod N} (LLL)."""
    rows = [[N] + [0] * d] + [[-pow(m, k, N)] + [0] * (k - 1) + [1] + [0] * (d - k) for k in range(1, d + 1)]
    L = fmpz_mat(rows).lll()
    cands = []
    red = [[int(L[i, j]) for j in range(d + 1)] for i in range(L.nrows())]
    combos = sorted(itertools.product((-1, 0, 1), repeat=L.nrows()), key=lambda c: sum(map(abs, c)))
    for c in combos:
        if not any(c):
            continue
        h = [sum(ci * row[j] for ci, row in zip(c, red)) for j in range(d + 1)]
        while len(h) > 1 and h[-1] == 0:
            h.pop()
        if len(h) - 1 >= 2 and irreducible(h):
            cands.append(h)
    for h1, h2 in itertools.combinations(cands, 2):
        if fpoly(h1).resultant(fpoly(h2)) != 0:
            return h1, h2
    raise RuntimeError("no lattice pair")


def evaluate(f: Sequence[int], g: Sequence[int], gens: List[List[int]], box, gens2: List[List[int]] = None,
             integral: bool = False) -> Dict:
    """Exact evaluation over the box: candidate a -> A = sum a_j gens_j (K side), A' = sum a_j gens2_j (K' side,
    default = gens); W = Res(f, A), W' = Res(g, A'). With integral=True, W = Res(f, A)/lc(f)^deg A (the norm of the
    algebraic integer A(alpha)), likewise W'."""
    ff, gg = fpoly(f), fpoly(g)
    lf, lg = int(f[-1]), int(g[-1])
    gens2 = gens if gens2 is None else gens2
    logs, n_smooth, n_cand = [], 0, 0
    for a in box:
        if math.gcd(*a) != 1:
            continue
        A = poly_comb(gens, a)
        A2 = poly_comb(gens2, a)
        if not any(A) or not any(A2):
            continue
        w = int(ff.resultant(fpoly(A)))
        w2 = int(gg.resultant(fpoly(A2)))
        if integral:
            qa, ra = divmod(w, lf ** (len(A) - 1))
            qb, rb = divmod(w2, lg ** (len(A2) - 1))
            assert ra == 0 and rb == 0
            w, w2 = qa, qb
        if w == 0 or w2 == 0:
            continue
        n_cand += 1
        logs.append((abs(w) * abs(w2)).bit_length())
        if smooth(w2) and smooth(w):
            n_smooth += 1
    return {"n_cand": n_cand, "median_log2F": statistics.median(logs) if logs else None, "yield": n_smooth}


def box2():
    return ((a, b) for b in range(1, E2D + 1) for a in range(-E2D, E2D + 1))


def box3():
    return ((a, b, c) for c in range(1, E3D + 1) for b in range(-E3D, E3D + 1) for a in range(-E3D, E3D + 1))


def omega_basis(f: Sequence[int]) -> List[List[int]]:
    """Integral elements omega_0 = 1, omega_k = c_d x^k + c_{d-1} x^{k-1} + ... + c_{d-k+1} x (k < d)."""
    d = len(f) - 1
    out = [[1]]
    for k in range(1, d):
        h = [0] * (k + 1)
        for t in range(1, k + 1):
            h[t] = int(f[d - k + t])
        out.append(h)
    return out


def lll_rank3(f, g, N: int, m: int, rnd: random.Random, n_sample: int = 1500):
    """Design R3-lll. Linked lattice of integral elements Lambda = {(y, y') in Z[omega] x Z[omega'] :
    phi(y) = phi'(y') mod N} (rank d + d'), embedded with K'-side weight 2^w; integer LLL (rounded embedding) and
    the 3 shortest reduced vectors as generators (y_j, y'_j) -- different polynomials on the two sides, linked mod N.
    Keep the weight with the smallest sampled median log2|F| (heuristic search; counts are exact)."""
    om, om2 = omega_basis(f), omega_basis(g)
    ev = lambda h: sum(c * pow(m, k, N) for k, c in enumerate(h)) % N
    v, v2 = [ev(h) for h in om], [ev(h) for h in om2]
    d, d2 = len(om), len(om2)
    # basis of Lambda in coordinates (b_0..b_{d-1}, b'_0..b'_{d2-1}); b'_0 absorbs the congruence (v'_0 = 1)
    basis = []
    for k in range(d):
        e = [0] * (d + d2)
        e[k] = 1
        e[d] = v[k] % N
        basis.append(e)
    for k in range(1, d2):
        e = [0] * (d + d2)
        e[d + k] = 1
        e[d] = (-v2[k]) % N
        basis.append(e)
    e = [0] * (d + d2)
    e[d] = N
    basis.append(e)
    # this spans Lambda (d + d2 + 1 generators, rank d + d2); LLL removes the dependency
    ctx.prec = 200
    rf = [z for z, _ in fpoly(f).complex_roots()]
    rg = [z for z, _ in fpoly(g).complex_roots()]

    def emb(poly, roots, scale):
        out = []
        for z in roots:
            s = acb(0)
            for c in reversed(poly):
                s = s * z + int(c)
            if z.imag.contains(0):
                out.append(int(float(s.real.mid()) * scale))
            elif float(z.imag.mid()) > 0:
                out.append(int(float(s.real.mid()) * scale))
                out.append(int(float(s.imag.mid()) * scale))
        return out

    def polys_of(b):
        y = poly_comb(om, b[:d])
        y2 = poly_comb(om2, b[d:])
        return y, y2
    sample = [(rnd.randrange(-E3D, E3D + 1), rnd.randrange(-E3D, E3D + 1), rnd.randrange(1, E3D + 1))
              for _ in range(n_sample)]
    best = None
    for w in range(-16, 17, 2):
        sK, sK2 = 2.0 ** 30, 2.0 ** (30 + w)
        rows = []
        for b in basis:
            y, y2 = polys_of(b)
            rows.append(emb(y, rf, sK) + emb(y2, rg, sK2))
        L, T = fmpz_mat(rows).lll(transform=True)
        idx = [i for i in range(L.nrows()) if any(int(L[i, j]) for j in range(L.ncols()))]
        idx.sort(key=lambda i: sum(int(L[i, j]) ** 2 for j in range(L.ncols())))
        gensK, gensK2 = [], []
        for i in idx[:3]:
            b = [sum(int(T[i, t]) * basis[t][c] for t in range(len(basis))) for c in range(d + d2)]
            y, y2 = polys_of(b)
            gensK.append(gtrim(y))
            gensK2.append(gtrim(y2))
        if rank_of(gensK) < 3 or rank_of(gensK2) < 3:  # both projections must be injective (rank-3 module)
            continue
        evs = evaluate(f, g, gensK, sample, gens2=gensK2, integral=True)
        if evs["n_cand"] < n_sample // 2:
            continue
        if best is None or evs["median_log2F"] < best[0]:
            best = (evs["median_log2F"], w, gensK, gensK2)
    return best


def rank_of(polys) -> int:
    n = max(len(p) for p in polys)
    return sympy.Matrix([list(p) + [0] * (n - len(p)) for p in polys]).rank()


def gtrim(h):
    h = list(h)
    while len(h) > 1 and h[-1] == 0:
        h.pop()
    return h


def main() -> None:
    rows = []
    with open(OUT, "w") as fh:
        for n in (40, 50, 60):
            for s in range(3):
                rnd = random.Random(f"{SEED}-E2-{n}-{s}")
                N, p, q = make_N(rnd, n)
                f_rat, g_rat, m = base_m(N, 3)
                h1, h2 = lattice_pair(N, m, 3)
                res = {"n": n, "seed": s, "N": N, "f_rat": f_rat, "m": m, "h1": h1, "h2": h2}
                t0 = time.time()
                res["R2-rat"] = evaluate(f_rat, g_rat, [[1], [0, -1]], box2())
                res["R2-nn"] = evaluate(h1, h2, [[1], [0, -1]], box2())
                P3 = [[1], [0, 1], [0, 0, 1]]
                res["R3-rat"] = evaluate(f_rat, g_rat, P3, box3())
                res["R3-nn"] = evaluate(h1, h2, P3, box3())
                best = lll_rank3(h1, h2, N, m, rnd)
                res["R3-lll_w"] = best[1]
                res["R3-lll_gens"] = [best[2], best[3]]
                res["R3-lll"] = evaluate(h1, h2, best[2], box3(), gens2=best[3], integral=True)
                res["secs"] = round(time.time() - t0, 1)
                rows.append(res)
                fh.write(json.dumps(res) + "\n")
                fh.flush()
                print(n, s, {k: (v["yield"], v["median_log2F"], v["n_cand"]) for k, v in res.items()
                             if isinstance(v, dict)}, "w", best[1], res["secs"], "s", flush=True)
    viol6 = viol7 = 0
    for res in rows:
        y2 = max(res["R2-rat"]["yield"], res["R2-nn"]["yield"])
        y3 = max(res[k]["yield"] for k in ("R3-rat", "R3-nn", "R3-lll"))
        b2 = min(res["R2-rat"]["median_log2F"], res["R2-nn"]["median_log2F"]) if True else None
        best2 = "R2-rat" if res["R2-rat"]["yield"] >= res["R2-nn"]["yield"] else "R2-nn"
        m2 = res[best2]["median_log2F"]
        viol6 += y3 > 2 * y2
        viol7 += any(res[k]["median_log2F"] < m2 - 2 for k in ("R3-rat", "R3-nn", "R3-lll"))
        print("n", res["n"], "seed", res["seed"], "best2", best2, "y2", y2, "y3", y3, "ratio",
              round(y3 / max(y2, 1), 3))
    print(f"P6 violations {viol6}/{len(rows)}; P7 violations {viol7}/{len(rows)}")


if __name__ == "__main__":
    main()
