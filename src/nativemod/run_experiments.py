"""nativemod experiments (HARNESS), as pre-registered in notes/nativemod_prereg.md.

    .venv/bin/python -m src.nativemod.run_experiments e1 | e2 | e3 [--quick]
Raw output: experiments/raw/nativemod_<exp>.jsonl
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from fractions import Fraction
from math import gcd

import gmpy2

from src.common.instances import balanced_semiprime
from src.nativemod import ecalg
from src.nativemod import harness as H
from src.nativemod import reductions as R

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, f"nativemod_{name}.jsonl"), "w")


def _odd_semiprimes(limit):
    ps = [p for p in range(3, limit) if gmpy2.is_prime(p)]
    out = []
    for i, p in enumerate(ps):
        for q in ps[i + 1:]:
            if p * q > limit:
                break
            out.append((p * q, p, q))
    return sorted(out)


# ---------------------------------------------------------------------------------------------------------------
def e1(quick=False):
    f = _out("e1")
    summary = {}
    # E1a
    lim = 600 if quick else 4000
    mism = {k: 0 for k in (2, 4, 6, 8, 10, 14, 12, 16, 18, 20, 22, 26)}
    cnt = {k: 0 for k in mism}
    t0 = time.time()
    for N, p, q in _odd_semiprimes(lim):
        for k in (2, 4, 6, 8, 10, 14):
            C = H.hurwitz_moment(N, k)
            exp = 2 * (N + q) if k == 2 else -2 * (1 + p ** (k - 1))
            cnt[k] += 1
            if C != exp:
                mism[k] += 1
                f.write(json.dumps({"e": "E1a", "N": N, "k": k, "C": str(C), "exp": exp}) + "\n")
        if N <= (600 if quick else 1500):
            for k in (12, 16, 18, 20, 22, 26):
                C = H.hurwitz_moment(N, k)
                exp = -2 * H.trace_TN_level1(N, k) - 2 * (1 + p ** (k - 1))
                cnt[k] += 1
                if C != exp:
                    mism[k] += 1
                    f.write(json.dumps({"e": "E1a", "N": N, "k": k, "C": str(C), "exp": exp}) + "\n")
    summary["E1a"] = {"checked": cnt, "mismatches": mism, "limit": lim, "time": round(time.time() - t0, 1)}
    # E1b
    rng = random.Random("nativemod-E1b")
    # D1 (disclosed deviation): registered N <= 10^4; PARI's mfheckemat(N^2) overflows a 512 MB stack for
    # N > ~8000 (it evaluates q-expansions up to ~N^2), so the 25 semiprimes per level are drawn from N <= 5000.
    sp = [x for x in _odd_semiprimes(5000)]
    res = {}
    for M in (11, 37, 43, 53):
        cand = [x for x in sp if x[0] % M]
        pick = rng.sample(cand, 25)
        bad = [0, 0, 0]
        for N, p, q in pick:
            TN, Tp, Tq = H.hecke_full(M, N), H.hecke_full(M, p), H.hecke_full(M, q)
            if TN != Tp * Tq:
                bad[0] += 1
            cp = H.pari.charpoly(TN)
            if int(H.pari.subst(cp, H.pari("x"), H.sigma(N))) != 0:
                bad[1] += 1
            I = H.pari.matid(int(H.pari.matsize(TN)[0]))
            lhs = TN * TN - H.hecke_full(M, N * N) - N * I
            rhs = p * H.hecke_full(M, q * q) + q * H.hecke_full(M, p * p)
            if lhs != rhs:
                bad[2] += 1
        res[M] = {"n": 25, "TN=TpTq_fail": bad[0], "sigma_eig_fail": bad[1], "BC_identity_fail": bad[2],
                  "dim": int(H.pari.matsize(H.hecke_full(M, 2))[0])}
    summary["E1b"] = res
    # E1c
    B = 600 if quick else 3000
    t1, t2 = H.theta_counts(H.L11_1, B), H.theta_counts(H.L11_2, B)
    fco = H.coeffs_11a1(B)
    a1, a2 = H.aut_order(H.L11_1), H.aut_order(H.L11_2)
    i_fail = sum(1 for n in range(1, B + 1) if t1[n] - t2[n] != 4 * fco[n])
    ii_fail = 0
    for n in range(1, B + 1):
        if n % 11 == 0:
            continue
        lhs = Fraction(t1[n], a1) + Fraction(t2[n], a2)
        if lhs / (Fraction(1, a1) + Fraction(1, a2)) != Fraction(12, 5) * H.sigma(n):
            ii_fail += 1
    z4 = H.theta_counts([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], B, even=False)
    d4 = H.theta_counts(H.D4_HURWITZ, B)
    a22 = H.theta_counts(H.A2A2, B)
    iii = {"Z4": sum(1 for n in range(1, B + 1, 2) if z4[n] != 8 * H.sigma(n)),
           "D4": sum(1 for n in range(1, B + 1, 2) if d4[n] != 24 * H.sigma(n)),
           "A2A2": sum(1 for n in range(1, B + 1) if n % 3 and a22[n] != 12 * H.sigma(n))}
    summary["E1c"] = {"B": B, "aut": [a1, a2], "i_fail": i_fail, "ii_fail_registered": ii_fail, "iii_fail": iii}
    # E1d (descriptive)
    rows = []
    for n in (list(range(2, 601, 1)) if not quick else list(range(2, 201))):
        m = H.merel_size(n)
        rows.append((n, m, m / ((12 * math.log(2) / math.pi ** 2) * H.sigma(n) * math.log(n))))
    for n, m, r in rows:
        f.write(json.dumps({"e": "E1d", "n": n, "size": m, "ratio": round(r, 5)}) + "\n")
    tail = [r for n, m, r in rows if n >= 300] or [r for n, m, r in rows]
    summary["E1d"] = {"n_max": rows[-1][0], "ratio_mean_tail": round(sum(tail) / len(tail), 4),
                      "ratio_min_tail": round(min(tail), 4), "ratio_max_tail": round(max(tail), 4)}
    f.write(json.dumps({"summary": summary}, default=str) + "\n")
    f.close()
    print(json.dumps(summary, indent=1, default=str))


# ---------------------------------------------------------------------------------------------------------------
def _chi_m4(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)


def _one(n):
    return 1


def eis_value(p, q, k, chi, chip):
    N = p * q
    return sum(chi(N // d) * chip(d) * d ** (k - 1) for d in (1, p, q, N))


def e2(quick=False):
    f = _out("e2")
    fams = {"(2,1,1)": (2, _one, _one), "(4,1,1)": (4, _one, _one), "(12,1,1)": (12, _one, _one),
            "(3,1,chi-4)": (3, _one, _chi_m4), "(3,chi-4,1)": (3, _chi_m4, _one),
            "(2,chi-4,chi-4)": (2, _chi_m4, _chi_m4)}
    sizes = (64, 128) if quick else (64, 128, 256)
    nper = 10 if quick else 40
    summ = {"E2a": {}, "E2b": {}, "E2b_null": {}}
    for name, (k, chi, chip) in fams.items():
        ok = 0
        tot = 0
        for n in sizes:
            for seed in range(nper):
                inst = balanced_semiprime(n, seed)
                p, q = inst.factors
                a = eis_value(p, q, k, chi, chip)
                d = R.eis_factor(inst.N, a, k, chi, chip)
                good = d in (p, q)
                ok += good
                tot += 1
                f.write(json.dumps({"e": "E2a", "fam": name, "n": n, "seed": seed, "ok": good}) + "\n")
        summ["E2a"][name] = f"{ok}/{tot}"
    for curve in ("37a1", "11a1", "32a2"):
        ok = tot = 0
        wins = []
        for n in ((64,) if quick else (64, 128)):
            for seed in range(nper):
                inst = balanced_semiprime(n, seed)
                p, q = inst.factors
                if any(inst.N % r == 0 for r in (2, 3, 11, 37)):
                    continue
                ap, aq = H.ap(curve, p), H.ap(curve, q)
                t = (ap * q * pow(q, -1, p) + aq * p * pow(p, -1, q)) % inst.N
                st = {}
                d = R.crt_trace_factor(inst.N, t, stats=st)
                good = d in (p, q)
                ok += good
                tot += 1
                wins.append(st.get("windows"))
                f.write(json.dumps({"e": "E2b", "curve": curve, "n": n, "seed": seed, "ok": good,
                                    "stats": st}) + "\n")
        summ["E2b"][curve] = {"ok": f"{ok}/{tot}", "max_windows": max(w for w in wins if isinstance(w, int))
                              if any(isinstance(w, int) for w in wins) else wins[:3]}
    # null
    rng = random.Random("nativemod-E2b-null")
    for n in ((64,) if quick else (64, 128)):
        hits = 0
        for seed in range(nper):
            inst = balanced_semiprime(n, 1000 + seed)
            t = rng.randrange(inst.N)
            d = R.crt_trace_factor(inst.N, t)
            hits += d is not None
        summ["E2b_null"][n] = f"{hits}/{nper}"
    f.write(json.dumps({"summary": summ}) + "\n")
    f.close()
    print(json.dumps(summ, indent=1))


# ---------------------------------------------------------------------------------------------------------------
def cmh_z(S, T, strata):
    """Cochran-Mantel-Haenszel z for binary S (in {+1,-1}) vs binary T (0/1) within strata."""
    from collections import defaultdict
    tab = defaultdict(lambda: [0, 0, 0, 0])  # n11 (S=+1,T=1), n1. (S=+1), n.1 (T=1), n
    for s, t, h in zip(S, T, strata):
        c = tab[h]
        c[0] += (s == 1 and t == 1)
        c[1] += (s == 1)
        c[2] += (t == 1)
        c[3] += 1
    num = 0.0
    var = 0.0
    for a, r1, c1, n in tab.values():
        if n < 2:
            continue
        num += a - r1 * c1 / n
        var += r1 * (n - r1) * c1 * (n - c1) / (n * n * (n - 1))
    return num / math.sqrt(var) if var > 0 else float("nan")


def cluster_z(S, T, strata, P, Q):
    """Prime-cluster-robust version (post hoc): same numerator, variance from residual products summed over
    clusters of instances sharing a prime (dyadic-robust: pairs sharing p or q counted)."""
    from collections import defaultdict
    groups = defaultdict(list)
    for i, h in enumerate(strata):
        groups[h].append(i)
    resid_s = [0.0] * len(S)
    resid_t = [0.0] * len(S)
    for h, idx in groups.items():
        ms = sum(1 for i in idx if S[i] == 1) / len(idx)
        mt = sum(T[i] for i in idx) / len(idx)
        for i in idx:
            resid_s[i] = (S[i] == 1) - ms
            resid_t[i] = T[i] - mt
    u = [resid_s[i] * resid_t[i] for i in range(len(S))]
    num = sum(u)
    byprime = defaultdict(float)
    for i in range(len(S)):
        byprime[P[i]] += u[i]
        byprime[Q[i]] += u[i]
    var_cl = sum(v * v for v in byprime.values()) - sum(x * x for x in u)  # each i counted in 2 clusters
    return num / math.sqrt(var_cl) if var_cl > 0 else float("nan")


def e3(quick=False, check_all=False):
    f = _out("e3c" if check_all else "e3")
    curve = "37a1"
    summ = {}
    nper = 2000 if quick else 20000
    for l in (3, 5):
        psi = H.divpoly(curve, l)
        disc = int(H.pari.poldisc(H.pari.Polrev(psi)))
        for n in (32, 48):
            rows = []
            excluded = 0
            splits = []
            null_exp = 0.0
            for seed in range(nper):
                inst = balanced_semiprime(n, seed)
                N = inst.N
                p, q = inst.factors
                if gcd(N, 6 * 37 * l * disc) != 1:
                    excluded += 1
                    continue
                chk = ecalg.Checker(N) if check_all else None
                try:
                    g = ecalg.monic_division_poly(psi, l, N)
                    st = ecalg.statistics(N, g, check_all=check_all, checker=chk)
                except ecalg.Split as e:
                    # classification: a coefficient vanishing at exactly one prime with no stated mechanism
                    splits.append({"seed": seed, "d": e.d, "where": e.where, "class": "ACCIDENTAL"})
                    if chk is not None:
                        null_exp += chk.count * (1 / p + 1 / q)
                    continue
                if chk is not None:
                    null_exp += chk.count * (1 / p + 1 / q)
                ap, aq = H.ap(curve, p), H.ap(curve, q)
                rows.append({"seed": seed, "N": N, "p": p, "q": q, "S1": st["S1"], "S2": st["S2"],
                             "S3": st["S3"], "S4": st["S4"], "aN": (ap * aq) % l,
                             "eps": H.kron(-1, p), "pl": H.kron(p, l), "Nl": H.kron(N, l)})
            # tests
            strata = [r["N"] % (4 * l) for r in rows]
            P = [r["p"] for r in rows]
            Q = [r["q"] for r in rows]
            res = {"n_used": len(rows), "excluded": excluded, "splits": splits, "n_splits": len(splits),
                   "null_expected_splits": round(null_exp, 1),
                   "S4_rate": sum(r["S4"] for r in rows) / max(1, len(rows))}
            targets = {}
            for c in range(l):
                targets[f"T1_aN={c}"] = (lambda r, c=c: int(r["aN"] == c), lambda r: True)
            targets["T2_aN=0"] = (lambda r: int(r["aN"] == 0), lambda r: True)
            targets["T3_eps"] = (lambda r: int(r["eps"] == 1), lambda r: r["N"] % 4 == 1)
            targets["T4_pl"] = (lambda r: int(r["pl"] == 1), lambda r: r["Nl"] == 1)
            zs = {}
            for sname in ("S1", "S2", "S3"):
                for tname, (tf, filt) in targets.items():
                    sub = [r for r in rows if filt(r) and r[sname] != 0]
                    S = [r[sname] for r in sub]
                    T = [tf(r) for r in sub]
                    h = [r["N"] % (4 * l) for r in sub]
                    z = cmh_z(S, T, h)
                    zc = cluster_z(S, T, h, [r["p"] for r in sub], [r["q"] for r in sub])
                    zs[f"{sname}|{tname}"] = (round(z, 3), round(zc, 3), len(sub))
            # positive control: CRT-built c with (c|p) = +1 iff target T2, (c|q) random fixed
            if n == 32:
                rngc = random.Random(f"nativemod-pc-{l}")
                S, T, hh = [], [], []
                for r in rows:
                    tval = int(r["aN"] == 0)
                    p, q = r["p"], r["q"]
                    # residue / non-residue mod p by target; mod q a fixed residue 1
                    while True:
                        cp = rngc.randrange(1, p)
                        if H.kron(cp, p) == (1 if tval else -1):
                            break
                    c = (cp * q * pow(q, -1, p) + 1 * p * pow(p, -1, q)) % (p * q)
                    S.append(int(gmpy2.jacobi(c, p * q)))
                    T.append(tval)
                    hh.append(r["N"] % (4 * l))
                res["positive_control_z"] = round(cmh_z(S, T, hh), 2)
            res["z"] = zs
            res["max_abs_z_registered"] = max(abs(v[0]) for v in zs.values())
            res["max_abs_z_cluster"] = max(abs(v[1]) for v in zs.values() if v[1] == v[1])
            summ[f"l={l},n={n}"] = res
            f.write(json.dumps({"cell": f"l={l},n={n}", **{k: v for k, v in res.items()}}) + "\n")
            print(f"l={l} n={n}: used {len(rows)} excl {excluded} splits {len(splits)} (null {null_exp:.1f}) S4 {res['S4_rate']:.4f} "
                  f"max|z| {res['max_abs_z_registered']:.2f} (cluster {res['max_abs_z_cluster']:.2f}) "
                  f"pc {res.get('positive_control_z')}", flush=True)
    f.write(json.dumps({"summary": {k: {kk: vv for kk, vv in v.items() if kk != 'z'} for k, v in summ.items()}},
                       default=str) + "\n")
    f.close()


# ---------------------------------------------------------------------------------------------------------------
def posthoc(quick=False):
    """Post hoc (not pre-registered; labelled as such in the notes)."""
    f = _out("posthoc")
    B = 3000
    t1, t2, t3 = (H.theta_counts(G, B) for G in (H.L11_1, H.L11_2, H.L11_3))
    fco = H.coeffs_11a1(B)
    a = [H.aut_order(G) for G in (H.L11_1, H.L11_2, H.L11_3)]
    w = [Fraction(1, x) for x in a]
    siegel_fail = sum(1 for n in range(1, B + 1) if n % 11 and
                      (w[0] * t1[n] + w[1] * t2[n] + w[2] * t3[n]) / sum(w) != Fraction(12, 5) * H.sigma(n))
    brandt_fail = sum(1 for n in range(1, B + 1) if n % 11 and 3 * t1[n] + 2 * t2[n] != 12 * H.sigma(n))
    c_fail = {name: sum(1 for n in range(1, B + 1) if n % 11 and
                        Fraction(t[n]) != Fraction(12, 5) * H.sigma(n) + c * fco[n])
              for name, t, c in (("L1", t1, Fraction(8, 5)), ("L2", t2, Fraction(-12, 5)), ("L3", t3, Fraction(18, 5)))}
    iso = [[int(H.pari.qfisom(H.pari.matrix(4, 4, [x for r in G for x in r]),
                                H.pari.matrix(4, 4, [x for r in K for x in r])) != 0)
            for K in (H.L11_1, H.L11_2, H.L11_3)] for G in (H.L11_1, H.L11_2, H.L11_3)]
    mer = []
    for n in (101, 211, 503, 1009, 2003, 4001, 8009):
        m = H.merel_size(n)
        mer.append({"n": n, "size": m, "size/(n log^2 n)": round(m / (n * math.log(n) ** 2), 4),
                    "size/((12log2/pi^2) n log n)": round(m / ((12 * math.log(2) / math.pi ** 2) * n * math.log(n)), 4)})
    summ = {"aut": a, "isometry_matrix": iso, "siegel_3lattice_fail": siegel_fail,
            "brandt_3to2_fail": brandt_fail, "theta_eq_E_plus_c_f_fail": c_fail, "merel_primes": mer, "B": B}
    f.write(json.dumps(summ) + "\n")
    f.close()
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    which = sys.argv[1]
    quick = "--quick" in sys.argv
    if which == "e3c":
        e3(quick, check_all=True)
    else:
        {"e1": e1, "e2": e2, "e3": e3, "posthoc": posthoc}[which](quick)
