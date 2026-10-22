"""Experiments E1-E4 of notes/ecount_prereg.md.  HARNESS driver.

Usage: .venv/bin/python -m src.ecount.run_experiments {e1|e2|e3|e4} [args]
Raw output: experiments/raw/ecount_e{1,2,3,4}.jsonl
"""
from __future__ import annotations

import json
import os
import random
import resource
import sys
import time
from math import gcd, isqrt, sqrt

import gmpy2
import numpy as np

from src.ecount import harness as H
from src.ecount.zn_curve import (Counters, du_affine_twist, du_malleability, du_projective_twist_repaired,
                                 sr2_supersingular, tw_trial)

RAW = "experiments/raw"


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, f"ecount_{name}.jsonl"), "w")


def _rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20  # macOS: bytes


def run_tw_fixed(N, a, b, M, rng, trials):
    outs = []
    ctr = Counters()
    for _ in range(trials):
        fac, out = tw_trial(N, a, b, M, rng, ctr)
        outs.append(out if fac is None else (out + "+"))
    return outs, ctr


# ------------------------------------------------------------------------------------------------ E1
def e1(sizes=(32, 64, 96, 128), seeds=range(40), trials=64):
    f = _out("e1")
    t0 = time.time()
    summary = {}
    for curve, (a, b) in H.CURVES.items():
        for n in sizes:
            agg = dict(inst=0, solved=0, exc=0, trials=0, struct=0, bonus_hits=0, rho_sum=0.0, bonus_sum=0.0,
                       null_struct=0, null_bonus=0, first=[])
            for seed in seeds:
                N, p, q = H.instance(curve, n, seed)
                M = H.count(a, b, p, q)
                dp, dq = H.local_data(a, b, p), H.local_data(a, b, q)
                rho, bonus, exc = H.rho_exact(dp, dq)
                rng = random.Random(f"ecount-e1|{curve}|{n}|{seed}")
                outs, ctr = run_tw_fixed(N, a, b, M, rng, trials)
                rng0 = random.Random(f"ecount-e1-null|{curve}|{n}|{seed}")
                nouts, _ = run_tw_fixed(N, a, b, M + 2, rng0, trials)
                struct = sum(o == "split+" for o in outs)
                bon = sum(o == "delta_gcd+" for o in outs)
                first = next((i + 1 for i, o in enumerate(outs) if o.endswith("+")), None)
                rec = dict(curve=curve, n=n, seed=seed, N=N, M=M, rho=rho, bonus=bonus, exc=exc,
                           ap=dp["ap"], aq=dq["ap"], twist_group_p=dp["twist_group"], twist_group_q=dq["twist_group"],
                           struct_splits=struct, bonus_splits=bon, trials=trials, first_split=first,
                           killed=sum(o == "killed" for o in outs), jac=sum(o == "jacobi_ne_-1" for o in outs),
                           affine=sum(o == "affine" for o in outs),
                           null_struct=sum(o == "split+" for o in nouts), null_bonus=sum(o == "delta_gcd+" for o in nouts),
                           ctr=ctr.as_dict())
                f.write(json.dumps(rec) + "\n")
                agg["inst"] += 1
                agg["solved"] += first is not None
                agg["exc"] += exc
                agg["trials"] += trials
                agg["struct"] += struct
                agg["bonus_hits"] += bon
                agg["rho_sum"] += rho * trials
                agg["bonus_sum"] += bonus * trials
                agg["null_struct"] += rec["null_struct"]
                agg["null_bonus"] += rec["null_bonus"]
                agg["first"].append(first)
            T = agg["trials"]
            pred = (agg["rho_sum"] + agg["bonus_sum"]) / T
            emp = (agg["struct"] + agg["bonus_hits"]) / T
            se = sqrt(pred * (1 - pred) / T)
            firsts = [x for x in agg["first"] if x is not None]
            s = dict(curve=curve, n=n, inst=agg["inst"], solved=agg["solved"], exc=agg["exc"],
                     emp_rate=emp, rho_rate=agg["rho_sum"] / T, pred_rate=pred, z=(emp - pred) / se if se else 0,
                     bonus_emp=agg["bonus_hits"] / T, bonus_pred=agg["bonus_sum"] / T,
                     null_struct=agg["null_struct"], null_bonus=agg["null_bonus"],
                     mean_first=sum(firsts) / len(firsts) if firsts else None, max_first=max(firsts) if firsts else None)
            summary[f"{curve}|{n}"] = s
            print(json.dumps(s), flush=True)
    # designed negative: CM1 with p = q = 3 mod 4
    a, b = H.CURVES["CM1_x3mx"]
    for n in sizes:
        tw_s = sr2_s = 0
        for seed in range(20):
            N, p, q = H.instance_cm1_33(n, seed)
            M = H.count(a, b, p, q)
            assert M == (p + 1) * (q + 1)
            outs, _ = run_tw_fixed(N, a, b, M, random.Random(f"ecount-e1-cm33|{n}|{seed}"), trials)
            st = sum(o == "split+" for o in outs)
            tw_s += st > 0
            g = sr2_supersingular(N, M)
            sr2_s += g in (p, q)
            f.write(json.dumps(dict(curve="CM1_33", n=n, seed=seed, N=N, M=M, struct_splits=st,
                                    bonus_splits=sum(o == "delta_gcd+" for o in outs), sr2=g in (p, q))) + "\n")
        s = dict(curve="CM1_33", n=n, inst=20, tw_instances_with_struct_split=tw_s, sr2_solved=sr2_s)
        summary[f"CM1_33|{n}"] = s
        print(json.dumps(s), flush=True)
    f.close()
    print(json.dumps(dict(time_s=time.time() - t0, rss_mb=_rss_mb())))
    json.dump(summary, open(os.path.join(RAW, "ecount_e1_summary.json"), "w"), indent=1)


# ------------------------------------------------------------------------------------------------ E2
def e2(bbits=(12, 14, 16, 18)):
    f = _out("e2")
    t0 = time.time()
    for curve, (a, b) in H.CURVES.items():
        for bb in bbits:
            P = [p for p in H.primes_window(bb) if H.disc_part(a, b) % p]
            D = [H.local_data(a, b, p) for p in P]
            k = len(P)
            pa = np.array(P, dtype=np.int64)
            A = np.array([d["A"] for d in D], dtype=np.int64)
            At = np.array([d["At"] for d in D], dtype=np.int64)
            r = np.array([d["r"] for d in D], dtype=np.int64)
            g1 = np.array([d["twist_group"][0] for d in D], dtype=np.int64)
            g2 = np.array([d["twist_group"][1] if len(d["twist_group"]) > 1 else 1 for d in D], dtype=np.int64)
            ss = np.array([d["ap"] == 0 for d in D])
            pairs = exc = exc_ss2 = exc_ss1 = exc_ord = low = 0
            rhos = []
            low_nonexc = 0
            rho_min_nonexc = 1.0
            exc_examples = []
            for i in range(k - 1):
                j = slice(i + 1, k)
                M = A[i] * A[j]
                e_i = (M % g1[i]) == 0
                e_j = (M % g1[j]) == 0
                ex = e_i & e_j
                Hi = np.gcd(g1[i], M) * np.gcd(g2[i], M)
                Hj = np.gcd(g1[j], M) * np.gcd(g2[j], M)
                num = (A[i] - 1 - r[i]) * (At[j] - Hj) + (At[i] - Hi) * (A[j] - 1 - r[j])
                rho = num / (4.0 * pa[i] * pa[j])
                thr = 0.125 - 4.0 / np.sqrt(pa[i])
                lo = rho < thr
                pairs += k - 1 - i
                exc += int(ex.sum())
                both = ss[i] & ss[j]
                one = ss[i] ^ ss[j]
                exc_ss2 += int((ex & both).sum())
                exc_ss1 += int((ex & one).sum())
                exc_ord += int((ex & ~(ss[i] | ss[j])).sum())
                low += int(lo.sum())
                low_nonexc += int((lo & ~ex).sum())
                if (~ex).any():
                    rho_min_nonexc = min(rho_min_nonexc, float(rho[~ex].min()))
                rhos.append(rho)
                if ex.any() and len(exc_examples) < 5:
                    for jj in np.nonzero(ex & ~(ss[i] | ss[j]))[0][:5 - len(exc_examples)]:
                        q = int(pa[i + 1 + jj])
                        dq = D[i + 1 + jj]
                        exc_examples.append(dict(p=int(pa[i]), q=q, Ap=int(A[i]), Atp=int(At[i]),
                                                 Gp=D[i]["twist_group"], Aq=dq["A"], Atq=dq["At"], Gq=dq["twist_group"]))
            allr = np.concatenate(rhos)
            rec = dict(curve=curve, b=bb, primes=k, pairs=pairs, exc=exc, exc_frac=exc / pairs,
                       exc_both_supersingular=exc_ss2, exc_one_supersingular=exc_ss1, exc_ordinary=exc_ord,
                       supersingular_primes=int(ss.sum()), low_rho=low, low_rho_nonexc=low_nonexc,
                       low_rho_nonexc_frac=low_nonexc / pairs,
                       rho_median=float(np.median(allr)), rho_mean=float(allr.mean()), rho_min_nonexc=rho_min_nonexc,
                       exc_ordinary_examples=exc_examples)
            f.write(json.dumps(rec) + "\n")
            print(json.dumps({k2: v for k2, v in rec.items() if k2 != "exc_ordinary_examples"}), flush=True)
    f.close()
    print(json.dumps(dict(time_s=time.time() - t0, rss_mb=_rss_mb())))


# ------------------------------------------------------------------------------------------------ E3
def e3():
    f = _out("e3")
    t0 = time.time()
    # (a) affine twist reduction
    for curve in H.NONCM:
        a, b = H.CURVES[curve]
        for n in (32, 64, 128):
            ok2 = ok1 = 0
            for seed in range(20):
                N, p, q = H.instance(curve, n, seed)
                Ea = H.affine_count(a, b, p, q)
                d = 2
                while gmpy2.jacobi(d, N) != -1:
                    d += 1
                g2 = du_affine_twist(N, Ea, H.twist_count(a, b, d, p, q, True), d)
                d1 = 2
                while not (gmpy2.legendre(d1, p) == -1 and gmpy2.legendre(d1, q) == -1):
                    d1 += 1
                g1 = du_affine_twist(N, Ea, H.twist_count(a, b, d1, p, q, True), d1)
                ok2 += g2 in (p, q)
                ok1 += g1 in (p, q)
            rec = dict(part="a", curve=curve, n=n, inst=20, case2_ok=ok2, case1_ok=ok1)
            f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)
    # (b) doubly supersingular counterexample (affine)
    a, b = H.CURVES["CM1_x3mx"]
    for n in (32, 64, 128):
        eqN = 0
        for seed in range(20):
            N, p, q = H.instance_cm1_33(n, seed)
            vals = [H.affine_count(a, b, p, q)] + [H.twist_count(a, b, d, p, q, True) for d in (2, 3, 5, 6, 7)]
            eqN += all(v == N for v in vals)
        rec = dict(part="b", curve="CM1_33", n=n, inst=20, all_affine_counts_equal_N=eqN)
        f.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    # (d) projective twist, repaired
    for curve in H.NONCM:
        a, b = H.CURVES[curve]
        for n in (32, 64, 96):
            ok = 0
            ctr = Counters()
            t1 = time.time()
            for seed in range(20):
                N, p, q = H.instance(curve, n, seed)
                d = 2
                while gmpy2.jacobi(d, N) != -1:
                    d += 1
                g = du_projective_twist_repaired(N, H.count(a, b, p, q), H.twist_count(a, b, d, p, q, False), d, ctr)
                ok += g in (p, q)
            rec = dict(part="d", curve=curve, n=n, inst=20, ok=ok, ctr=ctr.as_dict(), time_s=time.time() - t1)
            f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)
    # (c) malleability
    for curve in H.NONCM:
        a, b = H.CURVES[curve]
        for n in (48, 64, 96):
            ok = 0
            dMs, tried, inwin = [], [], []
            ctr = Counters()
            t1 = time.time()
            for seed in range(20):
                N, p, q = H.instance(curve, n, seed)
                M = H.count(a, b, p, q)
                res = du_malleability(N, H.factorint(M), ctr)
                ok += res["factor"] in (p, q)
                dMs.append(res["divisors_total"])
                tried.append(res["tried"])
                inwin.append(res["divisors_in_window"])
            rec = dict(part="c", curve=curve, n=n, inst=20, ok=ok, dM_max=max(dMs), dM_median=sorted(dMs)[10],
                       tried_max=max(tried), tried_median=sorted(tried)[10], in_window_max=max(inwin),
                       ctr=ctr.as_dict(), time_s=time.time() - t1)
            f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)
    f.close()
    print(json.dumps(dict(time_s=time.time() - t0, rss_mb=_rss_mb())))


# ------------------------------------------------------------------------------------------------ E4
def _bayes(counts):
    """counts: dict key -> [n_plus, n_minus]."""
    tot = sum(a + b for a, b in counts.values())
    return sum(max(a, b) for a, b in counts.values()) / tot - 0.5, tot


def e4(bb=16):
    f = _out("e4")
    t0 = time.time()
    pred_l = {l: H.chebotarev_adv_l(l) for l in (3, 5, 7)}
    pred_2 = H.chebotarev_adv_2()
    for curve in H.NONCM:
        a, b = H.CURVES[curve]
        P = [p for p in H.primes_window(bb) if H.disc_part(a, b) % p]
        A = {p: H.card(a, b, p) for p in P}
        R = {p: H.roots_f(a, b, p) for p in P}
        for l in (3, 5, 7):
            # per-prime class (p mod l, A_p mod l); ordered pairs p != q
            T = {}
            for p in P:
                key = (p % l, A[p] % l)
                T[key] = T.get(key, 0) + 1
            counts = {}
            for (d1, m1), c1 in T.items():
                for (d2, m2), c2 in T.items():
                    s1, s2 = gmpy2.legendre(d1, l), gmpy2.legendre(d2, l)
                    if s1 != s2:
                        continue
                    w = c1 * c2 - (c1 if (d1, m1) == (d2, m2) else 0)
                    key = ((d1 * d2) % l, (m1 * m2) % l)
                    counts.setdefault(key, [0, 0])[0 if s1 == 1 else 1] += w
            adv, tot = _bayes(counts)
            Mzero = sum(v for (d, m), v in T.items() if m == 0) / len(P)
            rec = dict(curve=curve, l=l, primes=len(P), pairs=tot, adv=adv, pred_GL2=pred_l[l],
                       diff=adv - pred_l[l], frac_Ap_0_mod_l=Mzero)
            f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)
        # l = 2: beta = (Delta_f / p) on (Delta_f / N) = 1, predictor M mod 2 (= [r_p = r_q = 0])
        Df = -(4 * a ** 3 + 27 * b ** 2)
        T = {}
        for p in P:
            key = (int(gmpy2.legendre(Df, p)), int(A[p] % 2))
            T[key] = T.get(key, 0) + 1
        counts = {}
        for (s1, m1), c1 in T.items():
            for (s2, m2), c2 in T.items():
                if s1 != s2:
                    continue
                w = c1 * c2 - (c1 if (s1, m1) == (s2, m2) else 0)
                counts.setdefault((m1 * m2) % 2, [0, 0])[0 if s1 == 1 else 1] += w
        adv, tot = _bayes(counts)
        podd = counts.get(1, [0, 0])
        rec = dict(curve=curve, l=2, primes=len(P), pairs=tot, adv=adv, pred_S3=pred_2, diff=adv - pred_2,
                   frac_M_odd=(podd[0] + podd[1]) / tot, M_odd_beta_minus=podd[1],
                   roots_check=all((R[p] == 0) == (A[p] % 2 == 1) for p in P))
        f.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    f.close()
    print(json.dumps(dict(time_s=time.time() - t0, rss_mb=_rss_mb())))


if __name__ == "__main__":
    which = sys.argv[1]
    {"e1": e1, "e2": e2, "e3": e3, "e4": e4}.get(which, lambda: None)()


# ------------------------------------------------------------------------------------------------ POST-HOC (labelled)
def _adv_from_table(T, l, use_m=True):
    counts = {}
    for (d1, m1), c1 in T.items():
        for (d2, m2), c2 in T.items():
            s1, s2 = gmpy2.legendre(d1, l), gmpy2.legendre(d2, l)
            if s1 != s2:
                continue
            w = c1 * c2 - (c1 if (d1, m1) == (d2, m2) else 0)
            key = ((d1 * d2) % l, (m1 * m2) % l if use_m else 0)
            counts.setdefault(key, [0, 0])[0 if s1 == 1 else 1] += w
    return _bayes(counts)[0]


def posthoc(bb=16, reps=200):
    """PH1: E1 null with M' uniform in the Hasse window of M (replaces the mis-designed M + 2 control).
    PH2: E4 baselines (advantage from N mod l alone) and a synthetic-population null: populations of the same size
    with Frobenius classes iid Haar on GL2(F_l); distribution of the Bayes advantage."""
    f = _out("posthoc")
    # PH1
    for curve in H.NONCM:
        a, b = H.CURVES[curve]
        for n in (64, 128):
            hits = trials = 0
            for seed in range(40):
                N, p, q = H.instance(curve, n, seed)
                rng = random.Random(f"ecount-ph1|{curve}|{n}|{seed}")
                w = 8 * isqrt(isqrt(N)) * isqrt(N)
                Mw = N + rng.randrange(-w, w)
                outs, _ = run_tw_fixed(N, a, b, Mw, rng, 64)
                hits += sum(o == "split+" for o in outs)
                trials += 64
            rec = dict(ph="PH1", curve=curve, n=n, trials=trials, struct_splits=hits)
            f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)
    # PH2
    rng = random.Random("ecount-ph2")
    for l in (3, 5, 7):
        G = [(a_, b_, c_, d_) for a_ in range(l) for b_ in range(l) for c_ in range(l) for d_ in range(l)
             if (a_ * d_ - b_ * c_) % l]
        k = len([p for p in H.primes_window(bb) if H.disc_part(*H.CURVES["C2_37a1"]) % p])
        sims = []
        for _ in range(reps):
            T = {}
            for _i in range(k):
                a_, b_, c_, d_ = rng.choice(G)
                key = ((a_ * d_ - b_ * c_) % l, ((1 - a_) * (1 - d_) - b_ * c_) % l)
                T[key] = T.get(key, 0) + 1
            sims.append((_adv_from_table(T, l), _adv_from_table(T, l, use_m=False)))
        sims.sort()
        adv_s = sorted(s[0] for s in sims)
        base_s = sorted(s[1] for s in sims)
        rec = dict(ph="PH2-null", l=l, primes=k, reps=reps, adv_q05=adv_s[reps // 20], adv_median=adv_s[reps // 2],
                   adv_q95=adv_s[-reps // 20], base_median=base_s[reps // 2], base_q95=base_s[-reps // 20])
        f.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
        for curve in H.NONCM:
            a, b = H.CURVES[curve]
            P = [p for p in H.primes_window(bb) if H.disc_part(a, b) % p]
            T = {}
            for p in P:
                key = (p % l, H.card(a, b, p) % l)
                T[key] = T.get(key, 0) + 1
            adv = _adv_from_table(T, l)
            base = _adv_from_table(T, l, use_m=False)
            rank = sum(s < adv for s in adv_s) / reps
            rec = dict(ph="PH2", curve=curve, l=l, adv=adv, base_N_mod_l_only=base, increment=adv - base,
                       null_rank=rank)
            f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)
    f.close()


if __name__ == "__main__" and sys.argv[1] == "posthoc":
    posthoc()
