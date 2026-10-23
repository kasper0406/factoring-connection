"""HARNESS: pre-registered experiments E1-E5 of notes/tracedata_prereg.md.

Usage: .venv/bin/python -m src.tracedata.run_experiments e1|e2|e3a|e3b|e3c|e3d|e4|e5
"""
from __future__ import annotations

import json
import os
import random
import resource
import sys
import time
from math import gcd, sqrt

import gmpy2
import numpy as np

from src.common.instances import balanced_semiprime
from src.tracedata import decode as D
from src.tracedata import galois as G
from src.tracedata import harness as H
from src.tracedata import reductions as R

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, f"tracedata_{name}.jsonl"), "w")


def _rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20  # macOS: bytes


def _emit(f, rec):
    f.write(json.dumps(rec) + "\n")
    f.flush()
    print(json.dumps(rec), flush=True)


# ------------------------------------------------------------------------------------------------ E1
def e1(sizes=(32, 64, 128), seeds=range(30)):
    f = _out("e1")
    t0 = time.time()
    for name in ("C1_11a1", "C2_37a1", "C3_2_3", "CM1_32a2", "T7_26b1"):
        for n in sizes:
            cnt = dict(inst=0, skipped=0, bc2_ok=0, bc2_pred=0, bc3_ok=0, bc3_pred=0, mix_ok=0,
                       tor_ok=0, tor_n=0, cmqr_ok=0, cmqr_n=0, both_zero=0, aN_zero=0)
            for seed in seeds:
                inst = balanced_semiprime(n, seed)
                N, (p, q) = inst.N, inst.labels["factors"]
                if not H.good(N, name):
                    cnt["skipped"] += 1
                    continue
                cnt["inst"] += 1
                a, b = H.ap(name, p), H.ap(name, q)
                aN = a * b
                aN2 = H.a_pk(a, p, 2) * H.a_pk(b, q, 2)
                aN3 = H.a_pk(a, p, 3) * H.a_pk(b, q, 3)
                M = (p + 1 - a) * (q + 1 - b)
                both0 = (a == 0 and b == 0)
                cnt["both_zero"] += both0
                cnt["aN_zero"] += (aN == 0)
                g = R.bc2_factor(N, aN, aN2)
                ok = g is not None and N % g == 0
                cnt["bc2_ok"] += ok
                cnt["bc2_pred"] += (ok == (not both0))
                g = R.bc3_factor(N, aN, aN3)
                ok = g is not None and N % g == 0
                cnt["bc3_ok"] += ok
                cnt["bc3_pred"] += (ok == (aN != 0))
                g = R.mix_factor(N, sorted((a, b)), M)
                cnt["mix_ok"] += (g is not None and N % g == 0)
                if name in H.TORSION_M:
                    m = H.TORSION_M[name]
                    cnt["tor_n"] += 1
                    cnt["tor_ok"] += (R.s_mod_m_from_torsion(N, aN, m) == (p + q) % m)
                if name == "CM1_32a2" and N % 4 == 1:
                    cnt["cmqr_n"] += 1
                    cnt["cmqr_ok"] += (R.cm_split_bit(aN) == (p % 4 == 1))
            _emit(f, dict(part="abcef", curve=name, n=n, **cnt))
    # (d) CM2 on split-split N
    for n in sizes:
        ok = tried = 0
        seed = 0
        while tried < 30:
            inst = balanced_semiprime(n, 1000 + seed)
            seed += 1
            N, (p, q) = inst.N, inst.labels["factors"]
            if p % 4 != 1 or q % 4 != 1:
                continue
            tried += 1
            a, b = H.ap("CM1_32a2", p), H.ap("CM1_32a2", q)
            g3 = H.cm_weight3_ap_fast(p) * H.cm_weight3_ap_fast(q)
            fac = R.linked_pair_factor(N, a * b, g3, c=1, e=2)
            ok += (fac is not None and N % fac == 0)
        _emit(f, dict(part="d_cm2", n=n, tried=tried, ok=ok, seeds_scanned=seed))
    _emit(f, dict(time_s=time.time() - t0, rss_mb=_rss_mb()))
    f.close()


# ------------------------------------------------------------------------------------------------ E2
def e2():
    f = _out("e2")
    t0 = time.time()
    tot = 0.0
    for l in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
        Io = G.info_ordered(l)
        Ip = G.info_product(l)
        tot += Io
        _emit(f, dict(l=l, I_ord=Io, l2_I_ord=l * l * Io, I_prod=Ip, l3_I_prod=l ** 3 * Ip, cum_I_ord=tot))
    pt = {}
    for l in (3, 5, 7):
        pt[f"M_mod_{l}"] = G.pt_floor(G.dist_M_mod_l(l))
        pt[f"aN_mod_{l}"] = G.pt_floor(G.dist_aN_mod_l(l))
    pt["r_psi3"] = G.pt_floor(G.dist_r_psi3())
    pt["r_S3_cubic"] = G.pt_floor(G.dist_r_S3())
    _emit(f, dict(pt_floor=pt, pt_split_lower_bound={k: v / 2 for k, v in pt.items()},
                  adv_psi3_beta=G.adv_psi3_beta(), adv_aN_beta={l: G.adv_aN_beta(l) for l in (3, 5, 7)},
                  inv_ln2=1 / np.log(2), time_s=time.time() - t0))
    f.close()


# ------------------------------------------------------------------------------------------------ E3
ELLS_A = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
KS_A = [250, 1000, 4000, 10000]


def _sim_ordered(l, K, trials, rng):
    n = G.class_counts(l).astype(float)
    hit = 0
    for _ in range(trials):
        d = int(rng.integers(1, l))
        c = rng.multinomial(K, n[:, d] / n[:, d].sum())
        ll = c @ D._logn(l)
        ll[0] = -np.inf
        hit += int(np.argmax(ll)) == d
    return hit / trials


def e3a(nprimes=200, trials=2000):
    f = _out("e3a")
    t0 = time.time()
    rng = random.Random("tracedata-e3a-1")
    P = H.prime_list(rng, 30, nprimes)
    Kmax = max(KS_A)
    T = {}
    for i, p in enumerate(P):
        T[p] = H.fam_traces(p, Kmax)
    _emit(f, dict(stage="traces", primes=len(P), K=Kmax, time_s=time.time() - t0, rss_mb=_rss_mb()))
    srng = np.random.default_rng(12345)
    for K in KS_A:
        for l in ELLS_A:
            real = sum(D.rank_p(l, [t % l for t in T[p][:K]])[0] == p % l for p in P) / len(P)
            null = sum(D.rank_p(l, [t % l for t in T[P[(i + 1) % len(P)]][:K]])[0] == p % l
                       for i, p in enumerate(P)) / len(P)
            sim = _sim_ordered(l, K, trials, srng)
            chance = 1 / (l - 1)
            se_null = sqrt(chance * (1 - chance) / len(P))
            _emit(f, dict(K=K, l=l, real=real, sim=sim, diff=real - sim, P3a=abs(real - sim) <= 0.10,
                          null=null, chance=chance, null_z=(null - chance) / se_null,
                          P3b=abs(null - chance) <= 3 * se_null))
    _emit(f, dict(time_s=time.time() - t0, rss_mb=_rss_mb()))
    f.close()


ELLS_B = [3, 5, 7, 11, 13]
KS_B = [1000, 10000, 40000]


def _sim_product(l, K, trials, rng, null=False):
    LL = D._logL(l)
    Lp = np.exp(LL)
    hit = 0
    for _ in range(trials):
        Nm = int(rng.integers(1, l))
        d = int(rng.integers(1, l))
        s = (d + Nm * pow(d, -1, l)) % l
        if null:
            Nm2, d2 = int(rng.integers(1, l)), int(rng.integers(1, l))
            pr = Lp[d2, Nm2 * pow(d2, -1, l) % l]
        else:
            pr = Lp[d, Nm * pow(d, -1, l) % l]
        c = rng.multinomial(K, pr / pr.sum())
        res = np.repeat(np.arange(l), c)
        hit += D.rank_s(l, Nm, res)[0] == s
    return hit / trials


def e3b(nN=100, trials=2000):
    f = _out("e3b")
    t0 = time.time()
    rng = random.Random("tracedata-e3b-2")
    pairs = []
    while len(pairs) < nN:
        p, q = H.prime_list(rng, 30, 2)
        pairs.append((p, q))
    Kmax = max(KS_B)
    A = {}
    for (p, q) in pairs:
        tp, tq = H.fam_traces(p, Kmax), H.fam_traces(q, Kmax)
        A[(p, q)] = [x * y for x, y in zip(tp, tq)]
    _emit(f, dict(stage="traces", N=len(pairs), K=Kmax, time_s=time.time() - t0, rss_mb=_rss_mb()))
    srng = np.random.default_rng(54321)
    for K in KS_B:
        for l in ELLS_B:
            real = null = 0
            for i, (p, q) in enumerate(pairs):
                N = p * q
                s = (p + q) % l
                real += D.rank_s(l, N % l, [a % l for a in A[(p, q)][:K]])[0] == s
                other = pairs[(i + 1) % len(pairs)]
                null += D.rank_s(l, N % l, [a % l for a in A[other][:K]])[0] == s
            real /= len(pairs)
            null /= len(pairs)
            sim = _sim_product(l, K, trials, srng)
            sim_null = _sim_product(l, K, trials, srng, null=True)
            se = sqrt(max(sim_null * (1 - sim_null), 1e-9) / len(pairs))
            _emit(f, dict(K=K, l=l, real=real, sim=sim, diff=real - sim, P3c=abs(real - sim) <= 0.12,
                          null=null, sim_null=sim_null, null_z=(null - sim_null) / se,
                          P3d=abs(null - sim_null) <= 3 * se))
    _emit(f, dict(time_s=time.time() - t0, rss_mb=_rss_mb()))
    f.close()


ELLS_C = [3, 5, 7, 11, 13, 17, 19, 23, 29]


def e3c(nN=30, K=40000):
    f = _out("e3c")
    t0 = time.time()
    rng = random.Random("tracedata-e3c-3")
    ok = 0
    for i in range(nN):
        p, q = H.prime_list(rng, 28, 2)
        N = p * q
        tr = H.fam_traces(p, K)
        per_l = {l: D.rank_p(l, [t % l for t in tr]).index(p % l) for l in ELLS_C}
        fac = D.factor_from_ordered(N, tr, ELLS_C, top=2)
        good = fac is not None and N % fac == 0
        ok += good
        _emit(f, dict(i=i, N=N, factor_found=good, rank_of_truth=per_l))
    _emit(f, dict(success=ok, of=nN, P3e=ok >= 24, time_s=time.time() - t0, rss_mb=_rss_mb()))
    f.close()


def e3d(nN=30, K=40000, nsim=200):
    f = _out("e3d")
    t0 = time.time()
    rng = random.Random("tracedata-e3d-4")
    ok = 0
    for i in range(nN):
        p, q = H.prime_list(rng, 28, 2)
        N = p * q
        tp, tq = H.fam_traces(p, K), H.fam_traces(q, K)
        aN = [x * y for x, y in zip(tp, tq)]
        per_l = {l: D.rank_s(l, N % l, [a % l for a in aN]).index((p + q) % l) for l in ELLS_C}
        fac = D.factor_from_products(N, aN, ELLS_C, top=3)
        good = fac is not None and N % fac == 0
        ok += good
        _emit(f, dict(i=i, N=N, factor_found=good, rank_of_truth=per_l))
    # Haar-simulated end-to-end rate: every l has the true s in its top 3
    srng = np.random.default_rng(777)
    sim_ok = 0
    for _ in range(nsim):
        all_in = True
        for l in ELLS_C:
            Lp = np.exp(D._logL(l))
            Nm, d = int(srng.integers(1, l)), int(srng.integers(1, l))
            s = (d + Nm * pow(d, -1, l)) % l
            pr = Lp[d, Nm * pow(d, -1, l) % l]
            c = srng.multinomial(K, pr / pr.sum())
            if s not in D.rank_s(l, Nm, np.repeat(np.arange(l), c))[:3]:
                all_in = False
                break
        sim_ok += all_in
    rate, srate = ok / nN, sim_ok / nsim
    _emit(f, dict(success=ok, of=nN, rate=rate, sim_rate=srate, P3f=abs(rate - srate) <= 0.25,
                  time_s=time.time() - t0, rss_mb=_rss_mb()))
    f.close()


# ------------------------------------------------------------------------------------------------ E4
def _bayes(counts):
    tot = sum(a + b for a, b in counts.values())
    return sum(max(a, b) for a, b in counts.values()) / tot - 0.5, tot


def e4(bb=16):
    f = _out("e4")
    t0 = time.time()
    model_classes = {(det, r) for (det, r) in G.psi3_table()}
    pred_aN = G.adv_aN_beta(3)
    pred_psi = G.adv_psi3_beta()
    for name in ("C2_37a1", "C5_92"):
        P = [p for p in H.census_primes(bb) if H.disc(name) % p]
        data = {p: (p % 3, H.ap(name, p) % 3, H.roots_psi3(name, p)) for p in P}
        # surjectivity proxy: all (det, #lines, trace) combinations of GL2(F_3) occur
        model_full = set()
        for g in G.gl2_elements(3):
            a, b, c, d = g
            model_full.add(((a * d - b * c) % 3, G.eigenlines_F3(g), (a + d) % 3))
        seen = {(v[0], v[2], v[1]) for v in data.values()}
        surj = bool(model_full <= seen)
        exc_b = sum(1 for (d, t, r) in data.values() if (r == 2 and d != 2))
        T = {}
        for p, v in data.items():
            T[v] = T.get(v, 0) + 1
        cA, cB = {}, {}
        for v1, c1 in T.items():
            for v2, c2 in T.items():
                if v1[0] != v2[0]:
                    continue
                w = c1 * c2 - (c1 if v1 == v2 else 0)
                beta = 0 if v1[0] == 1 else 1
                cA.setdefault(v1[1] * v2[1] % 3, [0, 0])[beta] += w
                cB.setdefault(v1[2] * v2[2], [0, 0])[beta] += w
        advA, tot = _bayes(cA)
        advB, _ = _bayes(cB)
        # P4b on pairs: r_N in {1,16} => p = 1 mod 3
        exc_pairs = sum(w for key, (w1, w2) in cB.items() if key in (1, 16) for w in [w2])
        _emit(f, dict(curve=name, primes=len(P), pairs=tot, surjective_mod3_proxy=surj,
                      adv_aN=advA, pred_aN=pred_aN, P4a_aN=bool(abs(advA - pred_aN) <= 0.03),
                      adv_psi3=advB, pred_psi3=pred_psi, P4a_psi3=bool(abs(advB - pred_psi) <= 0.03),
                      P4b_prime_exceptions=exc_b, P4b_pair_exceptions=exc_pairs,
                      rN_table={str(k): v for k, v in sorted(cB.items())}))
    _emit(f, dict(time_s=time.time() - t0, rss_mb=_rss_mb()))
    f.close()


# ------------------------------------------------------------------------------------------------ E5
def e5():
    from src.tracedata.xe_descent import xe_factor
    f = _out("e5")
    t0 = time.time()
    cases = [("11a1_census", -13392, -1080432, 54443, 62819), ("cx1_contrived", 27331, 2, 151, 181)]
    for (lab, a, b, p, q) in cases:
        N = p * q
        Ep = H.pari.ellinit([0, 0, 0, a, b])
        Ap, Aq = int(H.pari.ellcard(Ep, p)), int(H.pari.ellcard(Ep, q))
        M = Ap * Aq
        grp = {r: [int(x) for x in H.pari.ellgroup(Ep, r)] for r in (p, q)}
        res = xe_factor(N, a, b, M, random.Random(f"tracedata-e5-{lab}"), B=1000, max_trials=20)
        _emit(f, dict(case=lab, N=N, M=M, M_odd=M % 2 == 1, groups=grp, factor=res["factor"],
                      trials=res["trials"], outcomes=res["outcomes"], P5=res["factor"] is not None))
    _emit(f, dict(time_s=time.time() - t0, rss_mb=_rss_mb()))
    f.close()


# ------------------------------------------------------------------------------------------------ POST-HOC (labelled)
def _counts_mod(m):
    from itertools import product
    from math import gcd as _g
    n = np.zeros((m, m))
    for a, b, c, d in product(range(m), repeat=4):
        det = (a * d - b * c) % m
        if _g(det, m) == 1:
            n[(a + d) % m, det] += 1
    return n


def _info_mod(m):
    from math import gcd as _g
    n = _counts_mod(m)
    units = [d for d in range(m) if _g(d, m) == 1]
    I = W = 0.0
    for c in units:
        J = np.array([np.outer(n[:, d], n[:, c * pow(d, -1, m) % m]).ravel() for d in units])
        w = J.sum()
        I += w * G._mi(J)
        W += w
    return I / W


def posthoc():
    """PH1 identifiability classes; PH2 per-curve information by observation kind; PH3 prime-power moduli;
    PH4 E3a real - sim summary."""
    f = _out("posthoc")
    for l in (5, 7, 11, 13):
        for kind in ("single", "ordered", "unordered", "M", "product"):
            cls = {c: G.identifiable_classes(l, kind, c) for c in range(1, l)}
            _emit(f, dict(ph="PH1", l=l, kind=kind, n_classes={c: len(v) for c, v in cls.items()},
                          classes_c1=cls[1], classes_c2=cls[2]))
    for l in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        r = {k: G.info_kind(l, k) for k in ("single", "ordered", "unordered", "M", "product")}
        _emit(f, dict(ph="PH2", l=l, info_bits=r, l2_scaled={k: v * l * l for k, v in r.items()},
                      l3_product=r["product"] * l ** 3, l4_M=r["M"] * l ** 4))
    for m in (3, 9, 5, 25, 4, 8):
        _emit(f, dict(ph="PH3", m=m, I_ord=_info_mod(m)))
    rows = [json.loads(x) for x in open(os.path.join(RAW, "tracedata_e3a.jsonl")) if '"diff"' in x]
    d = [r["diff"] for r in rows if 0.05 < r["sim"] < 0.995]
    _emit(f, dict(ph="PH4", cells=len(d), mean_diff=float(np.mean(d)), min=min(d), max=max(d),
                  negative=sum(x < 0 for x in d)))
    f.close()


if __name__ == "__main__":
    {"e1": e1, "e2": e2, "e3a": e3a, "e3b": e3b, "e3c": e3c, "e3d": e3d, "e4": e4, "e5": e5, "posthoc": posthoc}[sys.argv[1]]()
