"""Run the pre-registered experiments E1-E3 (notes/geometry_prereg.md). Writes experiments/raw/geometry_e*.jsonl."""
import json
import os
import sys
import time

from src.common.instances import balanced_semiprime
from src.geometry import modcurve as mc
from src.geometry import curves as cv
from src.geometry import harness as H
from src.geometry import decoders as D

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")


def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def e1(out):
    res = {"a_N1": 0, "a_curves": 0, "a_old": 0, "manin": 0, "orbifold": 0}
    bad = []
    # (a)(i) N = 1
    for l in primes_upto(200):
        v = mc.trace_formula({}, l)
        res["a_N1"] += 1
        if v != 0:
            bad.append(("N1", l, v))
    # (a)(ii) genus-one levels
    for N, lab in ((11, "11a1"), (14, "14a1"), (15, "15a1"), (17, "17a1"), (19, "19a1"), (21, "21a1")):
        f = mc.factor_small(N)
        for l in primes_upto(200):
            if N % l == 0:
                continue
            a = mc.trace_formula(f, l)
            b = cv.ap(cv.CREMONA[lab], l)
            res["a_curves"] += 1
            if a != b:
                bad.append((lab, l, a, b))
    # (a)(iii) levels with old forms / two newforms
    combos = {22: [("11a1", 2)], 26: [("26a1", 1), ("26b1", 1)], 33: [("33a1", 1), ("11a1", 2)]}
    for N, parts in combos.items():
        f = mc.factor_small(N)
        for l in primes_upto(200):
            if N % l == 0:
                continue
            a = mc.trace_formula(f, l)
            b = sum(m * cv.ap(cv.CREMONA[lab], l) for lab, m in parts)
            res["a_old"] += 1
            if a != b:
                bad.append((N, l, a, b))
    # (b) Manin symbols, squarefree N with psi(N) <= 400, primes l <= 23
    t0 = time.time()
    for N in range(2, 400):
        f = mc.factor_small(N)
        if not mc.squarefree(f) or mc.psi(f) > 400:
            continue
        g = mc.genus(f)
        w = len(f)
        for l in primes_upto(23):
            if N % l == 0:
                continue
            ms = mc.manin_trace(N, l)
            pred = 2 * mc.trace_formula(f, l) + (l + 1) * (2 ** w - 1)
            res["manin"] += 1
            if ms["dim"] != 2 * g + 2 ** w - 1 or ms["trace"] != pred:
                bad.append(("manin", N, l, ms, pred))
    res["manin_seconds"] = round(time.time() - t0, 1)
    # (c) orbifold signature from explicit P^1(Z/N)
    for N in range(1, 2001):
        f = mc.factor_small(N)
        if N > 600 and not mc.squarefree(f):
            continue
        if N == 1:
            continue
        ob = mc.orbifold_from_p1(N)
        cf = {"mu": mc.psi(f), "nu2": mc.nu2(f), "nu3": mc.nu3(f), "nu_inf": mc.nu_inf(f), "genus": mc.genus(f)}
        res["orbifold"] += 1
        if ob != cf:
            bad.append(("orb", N, ob, cf))
    res["mismatches"] = len(bad)
    res["bad"] = bad[:20]
    out.write(json.dumps({"exp": "E1", **res}) + "\n")
    print("E1", res)


def e2(out):
    for n in (32, 64):
        stats = {k: 0 for k in ("G", "A", "P", "L", "K", "C")}
        qK = qC = 0
        for seed in range(200):
            inst = balanced_semiprime(n, seed)
            N = inst.N
            p, q = inst.labels["factors"]
            ok = lambda r: r is not None and r[0] * r[1] == N and 1 < r[0] < N
            stats["G"] += ok(D.from_genus(N, H.genus_oracle(p, q)))
            stats["A"] += ok(D.from_index(N, H.index_oracle(p, q)))
            stats["P"] += ok(D.from_pillai(N, H.pillai_oracle(p, q)))
            stats["L"] += ok(D.from_sigma(N, H.sigma_oracle(p, q)))
            r, qk = D.from_necklace(N, lambda k, rr: H.necklace_oracle(p, q, k, rr))
            stats["K"] += ok(r)
            qK = max(qK, qk)
            r, qc = D.from_catmap(N, lambda t, rr: H.catmap_oracle(p, q, t, rr))
            stats["C"] += ok(r)
            qC = max(qC, qc)
        rec = {"exp": "E2", "n": n, "instances": 200, "splits": stats, "max_queries_K": qK, "max_queries_C": qC}
        out.write(json.dumps(rec) + "\n")
        print(rec)


def e3(out):
    L = [2, 3, 5, 7, 11, 13]
    for n in (32, 64):
        tot = det = correct = consistent_true = 0
        nunk = []
        seed = 0
        while tot < 2000:
            inst = balanced_semiprime(n, seed)
            seed += 1
            N = inst.N
            if N % 4 != 1:
                continue
            p, q = inst.labels["factors"]
            counts = {l: H.x0_count_oracle(p, q, l) for l in L}
            eps, ncons, unk = D.decode_eps_from_x0_counts(N, counts)
            # check the true pattern is consistent (recompute)
            true = {D0: D.kron(D0, p) for D0 in unk}
            ok_true = True
            for l in L:
                terms = D.x0_count_terms(l)
                s = sum(w * 2 * (1 + true[D0]) for D0, w in terms if D.kron(D0, N) == 1)
                a_l = l + 1 - counts[l]
                if s != 2 * (mc.sigma1(l) - 4 - a_l):
                    ok_true = False
            consistent_true += ok_true
            tot += 1
            nunk.append(len(unk))
            if eps is not None:
                det += 1
                correct += eps == D.kron(-4, p)
        rec = {"exp": "E3", "n": n, "instances": tot, "eps_determined": det, "eps_correct": correct,
               "true_pattern_consistent": consistent_true, "mean_unknowns": sum(nunk) / len(nunk),
               "seeds_scanned": seed}
        out.write(json.dumps(rec) + "\n")
        print(rec)


if __name__ == "__main__":
    which = sys.argv[1:] or ["e1", "e2", "e3"]
    os.makedirs(RAW, exist_ok=True)
    for w in which:
        with open(os.path.join(RAW, f"geometry_{w}.jsonl"), "w") as out:
            {"e1": e1, "e2": e2, "e3": e3}[w](out)
