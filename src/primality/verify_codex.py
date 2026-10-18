"""HARNESS: exhaustive check of the GPT-6.1-sol r6 claims (notes/codex_r6_primality_lit.md).

Checked, for every pair of odd primes p < q with pq < bound and EVERY unit base:
  (1)  P(gcd(a^e - 1, N) proper) = A_p + A_q - 2 A_p A_q   (e = N-1, (N-1)/2, N+1, m)
  (2)  g_p = g_q = gcd(p-1, q-1); t = k iff u != v, t > k iff u = v
  §2 table: Fermat pass F; Fermat gcd A+B-2F; MR pass F C_k; nontrivial sqrt F(1-C_k);
       augmented MR A+B-F(1+C_k); SS pass F/2; SS residual (b vs J) gcd (A+B)/2 - F;
       SS both gcds A+B-2F (u=v) / A+B-3F/2 (u!=v)
  (8)  P(x_i = 1) = F 4^(min(i,k)-k);  P(x_i = -1) = F 4^(i-k) (i<k), 0 (k<=i<t)
  (10) P(x_i = -1 for some i >= 1) = F (4^k - 4)/(3 4^k)
  (11) P(b = -1) = 0 if u = v, F/4 if u != v       (b = a^((N-1)/2))
  Carmichael: P(MR gives no factor) <= 2^(1-r)  (exact engine, verified exhaustively in verify.py)
  (5)  strong-Lucas P-sampled formula; (7) torus table incl. Q^K normalisation; z^m = (V+sqrt(D)U)/(V-sqrt(D)U)
  weighted identity a^(Aq+Bp) = a^(AN+B) mod p (and one-sidedness mod q)
Run: .venv/bin/python -m src.primality.verify_codex [--bound 3000] [--workers 32]
"""
from __future__ import annotations

import argparse
import json
import random
import time
from fractions import Fraction as Fr
from math import gcd
from multiprocessing import Pool

import sympy

from src.common.instances import CARMICHAEL_SMALL, carmichael_chernick
from src.primality import formulas as F
from src.primality import ops
from src.primality.verify import semiprimes, qmul, qpow, conj


def check_pair(pq):
    p, q = pq
    N = p * q
    t, m = ops.split2(N - 1)
    u, v = F.v2(p - 1), F.v2(q - 1)
    k = min(u, v)
    g = gcd(p - 1, q - 1)
    out = {"p": p, "q": q, "bad": []}
    if not (gcd(N - 1, p - 1) == g == gcd(N - 1, q - 1)):
        out["bad"].append("(2) g")
    if (u != v and t != k) or (u == v and t <= k):
        out["bad"].append("(2) t")
    A, B = Fr(g, p - 1), Fr(g, q - 1)
    Fq = A * B
    Ck = Fr(4 ** k + 2, 3 * 4 ** k)
    phi = (p - 1) * (q - 1)
    exps = {"N-1": N - 1, "(N-1)/2": (N - 1) // 2, "N+1": N + 1, "m": m}
    c = {key: 0 for key in ["fpass", "fgcd", "mrpass", "sqrt", "aug", "sspass", "ssres", "ssboth",
                            "sqrtm1", "bm1"] + [f"e:{x}" for x in exps]}
    x1 = [0] * (t + 1)
    xm1 = [0] * (t + 1)
    for a in range(1, N):
        if gcd(a, N) != 1:
            continue
        mr = ops.miller_rabin(N, a)
        seq = mr["seq"]
        c["fpass"] += seq[t] == 1
        c["fgcd"] += 1 < gcd(seq[t] - 1, N) < N
        c["mrpass"] += mr["verdict"]
        c["sqrt"] += mr["sqrt_split"] is not None and 1 < mr["sqrt_split"] < N
        aug = any(1 < gcd(seq[i] - 1, N) < N or 1 < gcd(seq[i] + 1, N) < N for i in range(t)) or \
            1 < gcd(seq[t] - 1, N) < N
        c["aug"] += aug
        c["sqrtm1"] += any(seq[i] == N - 1 for i in range(1, t))
        for i in range(t + 1):
            x1[i] += seq[i] == 1
            xm1[i] += seq[i] == N - 1
        ss = ops.solovay_strassen(N, a)
        c["sspass"] += ss["verdict"]
        c["ssres"] += ss["residual_split"] is not None
        c["ssboth"] += any(l["class"] == "structural" for l in ss["leaks"])
        c["bm1"] += ss["euler"] == N - 1
        for key, e in exps.items():
            c[f"e:{key}"] += 1 < gcd(pow(a, e, N) - 1, N) < N
    pred = {"fpass": Fq, "fgcd": A + B - 2 * Fq, "mrpass": Fq * Ck, "sqrt": Fq * (1 - Ck),
            "aug": A + B - Fq * (1 + Ck), "sspass": Fq / 2, "ssres": (A + B) / 2 - Fq,
            "ssboth": A + B - 2 * Fq if u == v else A + B - Fr(3, 2) * Fq,
            "sqrtm1": Fq * Fr(4 ** k - 4, 3 * 4 ** k), "bm1": Fr(0) if u == v else Fq / 4}
    for key, e in exps.items():
        Ap, Aq = Fr(gcd(e, p - 1), p - 1), Fr(gcd(e, q - 1), q - 1)
        pred[f"e:{key}"] = Ap + Aq - 2 * Ap * Aq
    for key, pr in pred.items():
        if pr * phi != c[key]:
            out["bad"].append(f"{key}: count {c[key]} pred {pr * phi}")
    for i in range(t + 1):
        p1 = Fq * Fr(4) ** (min(i, k) - k)
        pm = Fq * Fr(4) ** (i - k) if i < k else Fr(0)
        if i < t and p1 * phi != x1[i]:
            out["bad"].append(f"(8)+ i={i}")
        if i == t and p1 * phi != x1[i]:
            out["bad"].append(f"(8)+ i=t")
        if i < t and pm * phi != xm1[i]:
            out["bad"].append(f"(8)- i={i}: {xm1[i]} vs {pm * phi}")
    out["counts"] = c
    return out


def check_lucas_codex(pq):
    """(5): P-sampled strong Lucas with U_(2^i d), V_(2^i d) (i<t) and U_L monitored."""
    p, q = pq
    N = p * q
    sel = ops.selfridge(N)
    if "P" not in sel:
        return None
    D = sel["D"]
    L = N + 1
    t, d = ops.split2(L)
    hp, hq = p - int(sympy.legendre_symbol(D % p, p)), q - int(sympy.legendre_symbol(D % q, q))

    def f(h, j):
        gd = gcd(d, h)
        if j == 0:
            return Fr(gd - 1, h - 1)
        if j <= F.v2(h):
            return Fr(gd * 2 ** (j - 1), h - 1)
        return Fr(0)
    Sp = sum(f(hp, j) for j in range(t + 1))
    Sq = sum(f(hq, j) for j in range(t + 1))
    pred = Sp + Sq - Sp * Sq - sum(f(hp, j) * f(hq, j) for j in range(t + 1))
    inv4 = pow(4, -1, N)
    cnt = val = 0
    for P in range(N):
        Q = (P * P - D) * inv4 % N
        if gcd(Q, N) != 1:
            continue
        val += 1
        r = ops.strong_lucas(N, P, Q)
        hit = any(1 < gcd(r["U"][i], N) < N or 1 < gcd(r["V"][i], N) < N for i in range(t)) or \
            1 < gcd(r["U"][t], N) < N
        cnt += hit
    return {"p": p, "q": q, "D": D, "ok": pred * val == cnt, "count": cnt, "pred": str(pred * val)}


def check_torus_codex(args):
    """(7) with Q^K normalisation for rho = root of x^2 - Px + Q (norm Q), and z^m identity."""
    (p, q), Dl, trials, seed = args
    N = p * q
    rng = random.Random(seed)
    cp, cq = int(sympy.legendre_symbol(Dl % p, p)), int(sympy.legendre_symbol(Dl % q, q))
    J = cp * cq
    s = p + q
    table = {(1, 1): (s, s - 2, 0), (-1, -1): (-s, -s - 2, s), (1, -1): (q - p + 2, q - p + 2, p - 1),
             (-1, 1): (p - q + 2, p - q + 2, q - 1)}
    eN1, eNJ, K = table[(cp, cq)]
    bad = []
    for _ in range(trials):
        while True:
            a0, a1 = rng.randrange(N), rng.randrange(N)
            nu = (a0 * a0 - Dl * a1 * a1) % N
            if gcd(nu, N) == 1 and gcd(a1, N) == 1:
                break
        A = (a0, a1)
        z = qmul(A, A, N, Dl)
        inv = pow(nu, -1, N)
        z = (z[0] * inv % N, z[1] * inv % N)

        def zp(e):
            return qpow(z, e, N, Dl) if e >= 0 else qpow(conj(z, N), -e, N, Dl)
        if qpow(z, N + 1, N, Dl) != zp(eN1):
            bad.append("z^(N+1)")
        if qpow(z, N - J, N, Dl) != zp(eNJ):
            bad.append("z^(N-J)")
        # rho = A has norm nu (plays Q); claim rho^(N+1) = Q^K rho^(eN1) (negative power via conj/nu)
        lhs = qpow(A, N + 1, N, Dl)
        if eN1 >= 0:
            r = qpow(A, eN1, N, Dl)
        else:
            r = qpow(conj(A, N), -eN1, N, Dl)
            sc = pow(inv, -eN1, N)
            r = (r[0] * sc % N, r[1] * sc % N)
        QK = pow(nu, K, N)
        if lhs != (r[0] * QK % N, r[1] * QK % N):
            bad.append("Q^K")
    return {"p": p, "q": q, "D": Dl, "pattern": [cp, cq], "bad": sorted(set(bad))}


def check_z_from_UV(seed=0, n=200):
    rng = random.Random(seed)
    bad = 0
    for _ in range(n):
        p, q = sorted(rng.sample(list(sympy.primerange(5, 300)), 2))
        N = p * q
        sel = ops.selfridge(N)
        if "P" not in sel:
            continue
        D = sel["D"]
        P = rng.randrange(N)
        Q = (P * P - D) * pow(4, -1, N) % N
        if gcd(Q, N) != 1:
            continue
        mexp = rng.randrange(1, 10 ** 6)
        U, V = ops.lucas_UV(N, P, Q, mexp)
        # in Z/N[w], w^2 = D: alpha^m = (V + w U)/2, beta^m = (V - w U)/2, z^m = alpha^m / beta^m
        num, den = (V % N, U % N), (V % N, (-U) % N)
        nd = (V * V - D * U * U) % N
        if gcd(nd, N) != 1:
            continue
        zm = qmul(num, conj(den, N), N, D)
        zm = (zm[0] * pow(nd, -1, N) % N, zm[1] * pow(nd, -1, N) % N)
        alpha = (P * pow(2, -1, N) % N, pow(2, -1, N))
        za = qmul(alpha, alpha, N, D)
        za = (za[0] * pow(Q, -1, N) % N, za[1] * pow(Q, -1, N) % N)
        bad += qpow(za, mexp, N, D) != zm
    return bad


def check_weighted(bound=400, seed=0):
    """a^(Aq+Bp) = a^(AN+B) mod p always; mod q iff ord_q(a) | (A-B)(p-1)."""
    rng = random.Random(seed)
    bad = 0
    onesided = total = 0
    for p, q in semiprimes(bound):
        N = p * q
        for _ in range(5):
            a = rng.randrange(2, N - 1)
            if gcd(a, N) != 1:
                continue
            A, B = rng.randrange(1, 20), rng.randrange(1, 20)
            x, y = pow(a, A * q + B * p, N), pow(a, A * N + B, N)
            bad += (x - y) % p != 0
            total += 1
            onesided += (x - y) % q != 0
            pred_q = ((A - B) * (p - 1)) % sympy.n_order(a, q) == 0
            bad += pred_q != ((x - y) % q == 0)
    return {"bad": bad, "total": total, "collision_mod_p_only": onesided}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bound", type=int, default=3000)
    ap.add_argument("--lucas-bound", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--out", default="experiments/raw/primality_verify_codex.json")
    a = ap.parse_args()
    t0 = time.time()
    res = {"config": vars(a)}
    from src.primality.verify import torus_cases
    with Pool(a.workers) as pool:
        rs = pool.map(check_pair, semiprimes(a.bound), chunksize=4)
        res["pairs"] = {"count": len(rs), "fail": [r for r in rs if r["bad"]]}
        print("pairs", len(rs), "fail", len(res["pairs"]["fail"]), f"{time.time()-t0:.0f}s", flush=True)
        rs = [r for r in pool.map(check_lucas_codex, semiprimes(a.lucas_bound)) if r]
        res["lucas5"] = {"count": len(rs), "fail": [r for r in rs if not r["ok"]]}
        print("lucas(5)", len(rs), "fail", len(res["lucas5"]["fail"]), flush=True)
        rs = pool.map(check_torus_codex, torus_cases(800, per_pattern=2, seed=7))
        pats = {}
        for r in rs:
            pats[str(r["pattern"])] = pats.get(str(r["pattern"]), 0) + 1
        res["torus7"] = {"count": len(rs), "by_pattern": pats, "fail": [r for r in rs if r["bad"]]}
        print("torus(7)", len(rs), pats, "fail", len(res["torus7"]["fail"]), flush=True)
    res["z_from_UV_bad"] = check_z_from_UV()
    res["weighted"] = check_weighted()
    # Carmichael bound via the exact engine (engine itself verified exhaustively in verify.py)
    carm = []
    for c in CARMICHAEL_SMALL:
        ps = [int(x) for x in sympy.primefactors(c)]
        ex = F.mr_exact(ps)
        nofactor = 1 - ex["prefix"]
        carm.append({"N": c, "r": len(ps), "P_nofactor": str(nofactor), "bound": str(Fr(2) ** (1 - len(ps))),
                     "ok": nofactor <= Fr(2) ** (1 - len(ps)), "liar_eq_nofactor": ex["liar"] == nofactor})
    for kb in (8, 12, 16, 20):
        for sd in range(3):
            inst = carmichael_chernick(kb, sd)
            ps = inst.labels["factors"]
            ex = F.mr_exact(ps)
            nofactor = 1 - ex["prefix"]
            carm.append({"N": inst.N, "r": 3, "P_nofactor": str(nofactor), "bound": "1/4",
                         "ok": nofactor <= Fr(1, 4), "liar_eq_nofactor": ex["liar"] == nofactor})
    res["carmichael_bound"] = {"count": len(carm), "fail": [c for c in carm if not (c["ok"] and c["liar_eq_nofactor"])],
                               "rows": carm}
    res["wall_s"] = time.time() - t0
    print("z_from_UV bad", res["z_from_UV_bad"], "weighted", res["weighted"], "carm fail",
          len(res["carmichael_bound"]["fail"]), flush=True)
    with open(a.out, "w") as f:
        json.dump(res, f, indent=1, default=str)


if __name__ == "__main__":
    main()
