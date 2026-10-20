"""POST-HOC diagnostics (disclosed; not pre-registered).

E1_distinct: replays run_E1's exact RNG stream and records the distinct full multisets A = prefix + suffix behind
the hits, to test whether the two P1 failures come from repeated detection of the same few members (the prefix
space is small: draws with replacement).
E1_exhaustive: for cells whose full M-multiset space is small, counts ALL members
#{A : |P_A(mu) mod+- N| <= V} exactly and compares with C(2B+M, M)(2V+1)/N.

Usage: PYTHONPATH=. .venv/bin/python -m src.nonexact.posthoc
"""
from __future__ import annotations

import json
import math
import os
import random
from collections import Counter
from itertools import combinations_with_replacement
from math import comb, gcd

import numpy as np

from src.common.instances import balanced_semiprime
from src.nonexact import core as C
from src.nonexact.run_experiment import RAW, load_cfg


def e1_distinct(cfg):
    c1 = cfg["E1"]
    rng = random.Random(c1["rng_seed"])
    out = []
    for cell in c1["cells"]:
        n, V = cell["n"], cell["V"]
        N = balanced_semiprime(n, c1["instance_seed"]).N
        mu = rng.randrange(2, N - 1)
        for B in cell["B"]:
            for s in cell["s"]:
                model = C.lin_cost_model(N, B, s, V)
                T = math.ceil(c1["target_hits"] / model["h"])
                Mtot = math.ceil(n * math.log(2) / math.log(2 * B)) + 2
                if T > c1["T_max"]:
                    continue
                ms = C.s_multisets_e(B, s)
                E = np.array([e for _, e in ms], dtype=np.int64)
                hitsets = Counter()
                hits = 0
                for _ in range(T):
                    pref = [rng.randint(-B, B) for _ in range(Mtot - s)]
                    c = 1
                    for a in pref:
                        c = c * (mu - a) % N
                    if gcd(c, N) != 1:
                        continue
                    w = C.lin_weights(c, mu, N, s)
                    if gcd(w[s], N) != 1:
                        continue
                    acc = np.full(E.shape[0], w[0] % N, dtype=np.int64)
                    for i in range(s):
                        acc = np.mod(acc + C._mulmod_vec(np.mod(E[:, i], N), w[i + 1], N), N)
                    cen = np.where(acc > N // 2, acc - N, acc)
                    idx = np.nonzero(np.abs(cen) <= V)[0]
                    for t in idx:
                        hits += 1
                        A = tuple(sorted(pref + list(ms[t][0])))
                        hitsets[A] += 1
                n_space = comb(2 * B + Mtot, Mtot)
                out.append({"n": n, "N": N, "mu": mu, "V": V, "B": B, "s": s, "M": Mtot, "T": T, "hits": hits, "distinct_members_hit": len(hitsets),
                            "max_multiplicity": max(hitsets.values()) if hitsets else 0,
                            "multiset_space": n_space,
                            "expected_members_in_space": n_space * (2 * V + 1) / N,
                            "top": [[list(a), m] for a, m in hitsets.most_common(3)]})
                print(out[-1]["n"], B, s, "hits", hits, "distinct", len(hitsets), "maxmult", out[-1]["max_multiplicity"],
                      "space", n_space, "exp members in space", round(out[-1]["expected_members_in_space"], 2), flush=True)
    return out


def _all_products_mod(mu: int, N: int, B: int, M: int) -> np.ndarray:
    """P_A(mu) mod N for all M-multisets A of [-B, B], by dynamic programming over the symbols (exact, int64,
    N < 2^40 via 20-bit limb multiplication)."""
    xs = [(mu - a) % N for a in range(-B, B + 1)]
    # layers[j] = array of products of j-multisets using symbols processed so far
    layers = [np.array([1], dtype=np.int64)] + [np.array([], dtype=np.int64) for _ in range(M)]
    for x in xs:
        new = [l.copy() for l in layers]
        # multisets that use symbol x at least once: build by repeated multiplication (unbounded multiplicity)
        for j in range(1, M + 1):
            if new[j - 1].size:
                add = C._mulmod_vec(new[j - 1], x, N)
                new[j] = np.concatenate([new[j], add])
        layers = new
    return layers[M]


def e1_exhaustive(distinct_rows, max_space=16_000_000):
    """Exact member count over the whole M-multiset space for each (n, B) with a small space."""
    done = set()
    out = []
    for r in distinct_rows:
        key = (r["n"], r["B"])
        if key in done or r["multiset_space"] > max_space:
            continue
        done.add(key)
        N, mu, V, B, M = r["N"], r["mu"], r["V"], r["B"], r["M"]
        vals = _all_products_mod(mu, N, B, M)
        assert vals.size == r["multiset_space"]
        cen = np.where(vals > N // 2, vals - N, vals)
        members = int(np.count_nonzero(np.abs(cen) <= V))
        out.append({"n": r["n"], "B": B, "M": M, "V": V, "space": int(vals.size), "members": members,
                    "expected": r["multiset_space"] * (2 * V + 1) / N})
        print("exhaustive", out[-1], flush=True)
    return out


def main():
    cfg = load_cfg()
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "exhaustive":
        rows = [json.loads(l) for l in open(os.path.join(RAW, "nonexact_E1_distinct_posthoc.jsonl"))]
        ex = e1_exhaustive(rows)
        with open(os.path.join(RAW, "nonexact_E1_exhaustive_posthoc.jsonl"), "w") as f:
            for r in ex:
                f.write(json.dumps(r) + "\n")
        return
    rows = e1_distinct(cfg)
    with open(os.path.join(RAW, "nonexact_E1_distinct_posthoc.jsonl"), "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    main()


def e3_reasons(cfg):
    """POST-HOC: replay run_E3 exactly and report, per violating pair, which P9 clause failed."""
    import itertools
    from sympy import primefactors
    from flint import fmpz_poly, nmod_poly
    from src.nonexact.run_experiment import _fibre_pairs, _monic_factor_of_degree
    c3 = cfg["E3"]
    rng = random.Random(c3["rng_seed"])
    H = c3["H"]
    out = []
    for (k, M, R) in c3["cells"]:
        pairs, total_pairs = _fibre_pairs(k, M, R, c3["pairs_cap"], rng)
        d = max(k, 2)
        reasons = Counter()
        n_pairs = 0
        for (A, A2) in pairs:
            P, P2 = C.from_roots(A), C.from_roots(A2)
            D = P - P2
            cD = C.coeffs(D, k)
            cand = []
            for delta in itertools.product(range(-2 * H, 2 * H + 1), repeat=k):
                diff = [cD[i] - delta[i] for i in range(k)]
                if not any(diff):
                    continue
                g = 0
                for x in diff:
                    g = gcd(g, x)
                for N in primefactors(abs(g)):
                    if N > max(4 * H + 1, c3["N_min"]):
                        cand.append((N, list(delta)))
            if not cand:
                continue
            N, delta = cand[rng.randrange(len(cand))]
            built = None
            for _ in range(c3["r_tries"]):
                r = [rng.randint(max(-H, x - H), min(H, x + H)) for x in delta]
                f = nmod_poly([int(x) % N for x in C.coeffs(P - fmpz_poly(r), M + 1)], N)
                g = _monic_factor_of_degree(f, d)
                if g is not None:
                    built = (r, fmpz_poly([int(x) for x in g.coeffs()]))
                    break
            if built is None:
                continue
            r, F = built
            r2 = [r[i] - delta[i] for i in range(k)]
            n_pairs += 1
            rels, E = C.ih_relations(A, A2, r, r2, N, k)
            v_bound = 2 * H * sum(R ** l for l in range(k))
            why = []
            if not all((a - b) % N == 0 for a, b in zip(C.coeffs(C.red_mod(P, F), d), r + [0] * (d - k))):
                why.append("member_P")
            if not all((a - b) % N == 0 for a, b in zip(C.coeffs(C.red_mod(P2, F), d), r2 + [0] * (d - k))):
                why.append("member_P2")
            for (b, u, v, nt) in rels:
                if (u - v) % N:
                    why.append("a_congruence")
                if nt != (int(E(b)) != 0) or nt != (u != v):
                    why.append("b_nontrivial")
                if abs(v) > v_bound:
                    why.append("c_vbound")
                if any(abs(b - a) > 2 * R for a in (A if b in A2 else A2)):
                    why.append("d_piece")
            distinct = {(u, v) if u > 0 else (-u, -v) for (_, u, v, nt) in rels if nt}
            n_supp = len(set(A) ^ set(A2))
            lb = C.ih_distinct_lower_bound(M, k, D.degree(), n_supp)
            good = sum(1 for (_, u, v, nt) in rels if nt)
            absvals = {abs(u) for (_, u, v, nt) in rels if nt}
            lb2 = -(-(n_supp - (k - 1)) // (2 * max(1, D.degree()))) if D.degree() > 0 else 1
            if len(distinct) < lb:
                why.append("e_distinct_lb")
            if len(distinct) < lb2:
                why.append("e_distinct_lb_corrected")
            for w in set(why):
                reasons[w] += 1
        out.append({"k": k, "M": M, "R": R, "pairs": n_pairs, "reasons": dict(reasons)})
        print(out[-1], flush=True)
    return out


def e2_baseline2(cfg):
    """POST-HOC: replay run_E2's RNG stream exactly; besides the registered per-coordinate baseline, draw a
    uniform-height baseline (every coefficient of Q uniform in [-||Q||_inf, ||Q||_inf]) from a separate RNG."""
    from flint import fmpz_poly
    from src.nonexact.run_experiment import balanced_semiprime as bs
    c2 = cfg["E2"]
    rng = random.Random(c2["rng_seed"])
    rng2 = random.Random(777)
    B, M = c2["B"], c2["M"]
    out = []
    for n_bits in c2["n"]:
        N = bs(n_bits, c2["instance_seed"]).N
        for d in c2["d"]:
            mu, mu2 = rng.randrange(2, N - 1), rng.randrange(2, N - 1)
            F = fmpz_poly([(-mu) % N, 1]) if d == 1 else C.quadratic_presentation(mu, mu2, N)
            for (j, k) in (c2["jk_d1"] if d == 1 else c2["jk_d2"]):
                qroots, b1n, b1r, b2n, b2r, nvec = 0, 0, 0, 0, 0, 0
                for _ in range(c2["T"]):
                    A = [rng.randint(-B, B) for _ in range(M)]
                    Lb = C.multiplier_lattice(C.from_roots(A), F, N, j, k)
                    Lb.sort(key=lambda v: sum(x * x for x in v))
                    for v in Lb:
                        cl = C.classify_vector(v, j, k, B)
                        if cl["Q_deg"] >= 1:
                            nvec += 1
                            qroots += int(cl["Q_int_roots"] > 0)
                            q = v[: j + 1]
                            for _b in range(c2["baseline_draws"]):
                                rq = [rng.randint(-abs(x), abs(x)) for x in q]
                                if any(rq[1:]):
                                    b1n += 1
                                    b1r += int(C.int_roots_count(rq) > 0)
                            h = C.height(q)
                            for _b in range(c2["baseline_draws"]):
                                rq = [rng2.randint(-h, h) for _ in q]
                                if any(rq[1:]):
                                    b2n += 1
                                    b2r += int(C.int_roots_count(rq) > 0)
                out.append({"n": n_bits, "d": d, "j": j, "k": k, "Q_with_int_root": qroots,
                            "baseline1_rate": b1r / b1n if b1n else None, "baseline2_rate": b2r / b2n if b2n else None,
                            "vectors_degQ": nvec, "observed_rate": qroots / nvec if nvec else None})
                print(out[-1], flush=True)
    return out


def rescore(cfg):
    """POST-HOC rescoring after review: P3 with the FROZEN statistic T(1 + p_bar)/hits (prereg text),
    P8 root-rate factors against the registered baseline, and the E1 best-s summary per (n, B)."""
    out = {"P3": [], "P8": [], "E1_best": []}
    rows = [json.loads(l) for l in open(os.path.join(RAW, "nonexact_E1.jsonl"))]
    best = {}
    for r in rows:
        if not r.get("run"):
            continue
        box = r["box_subset"]
        pbar = (sum(box) / len(box)) if (box and r["pred_p"] >= 1) else r["pred_p"]
        frozen = r["T"] * (1 + pbar) / r["hits"]
        model = r["pred_cost_per_hit"]
        out["P3"].append({"n": r["n"], "B": r["B"], "s": r["s"], "frozen_meas": frozen,
                          "frozen_over_model": frozen / model, "pass": 0.6 <= frozen / model <= 1.5})
        key = (r["n"], r["B"])
        if key not in best or frozen < best[key]["meas"]:
            best[key] = {"n": r["n"], "B": r["B"], "s": r["s"], "meas": frozen, "model": model,
                         "brute": r["brute_cost_per_hit"], "per_lattice_yield": r["hits"] / r["T"]}
    out["E1_best"] = list(best.values())
    for l in open(os.path.join(RAW, "nonexact_E2.jsonl")):
        r = json.loads(l)
        be = r["baseline_expected_count"]
        if be is None or be < 10:
            continue
        obs = r["Q_with_int_root"]
        f = be / obs if obs else float("inf")
        out["P8"].append({"n": r["n"], "d": r["d"], "j": r["j"], "obs": obs, "baseline_expected": be,
                          "factor": f, "pass": (1 / 3) <= (obs / be) <= 3})
    return out


def e3_theorem_check(cfg):
    """POST-HOC: replay run_E3 and evaluate the corrected Prop IH (hypothesis ||delta|| sum R^l < N;
    good points >= M - 2k + 3; distinct >= ceil((M-2k+3)/(2 deg D)))."""
    import itertools
    from sympy import primefactors
    from flint import fmpz_poly, nmod_poly
    from src.nonexact.run_experiment import _fibre_pairs, _monic_factor_of_degree
    c3 = cfg["E3"]
    rng = random.Random(c3["rng_seed"])
    H = c3["H"]
    out = []
    for (k, M, R) in c3["cells"]:
        pairs, _ = _fibre_pairs(k, M, R, c3["pairs_cap"], rng)
        d = max(k, 2)
        n_t, n_hyp, v_good, v_dist = 0, 0, 0, 0
        for (A, A2) in pairs:
            P, P2 = C.from_roots(A), C.from_roots(A2)
            D = P - P2
            cD = C.coeffs(D, k)
            cand = []
            for delta in itertools.product(range(-2 * H, 2 * H + 1), repeat=k):
                diff = [cD[i] - delta[i] for i in range(k)]
                if not any(diff):
                    continue
                g = 0
                for x in diff:
                    g = gcd(g, x)
                for N in primefactors(abs(g)):
                    if N > max(4 * H + 1, c3["N_min"]):
                        cand.append((N, list(delta)))
            if not cand:
                continue
            N, delta = cand[rng.randrange(len(cand))]
            built = None
            for _ in range(c3["r_tries"]):
                r = [rng.randint(max(-H, x - H), min(H, x + H)) for x in delta]
                f = nmod_poly([int(x) % N for x in C.coeffs(P - fmpz_poly(r), M + 1)], N)
                g = _monic_factor_of_degree(f, d)
                if g is not None:
                    built = r
                    break
            if built is None:
                continue
            r = built
            r2 = [r[i] - delta[i] for i in range(k)]
            n_t += 1
            if not C.ih_hypothesis(delta, R, N):
                continue
            n_hyp += 1
            rels, E = C.ih_relations(A, A2, r, r2, N, k)
            good = sum(1 for (_, u, v, nt) in rels if nt)
            distinct = {(u, v) if u > 0 else (-u, -v) for (_, u, v, nt) in rels if nt}
            v_good += int(good < C.ih_good_points_bound(M, k))
            v_dist += int(len(distinct) < C.ih_distinct_bound(M, k, D.degree()))
        out.append({"k": k, "M": M, "R": R, "tested": n_t, "hypothesis_holds": n_hyp,
                    "good_bound_violations": v_good, "distinct_bound_violations": v_dist})
        print(out[-1], flush=True)
    return out
