"""Pre-registered experiments E1, E1b, E2, E3 (notes/nonexact_prereg.md).

Usage: PYTHONPATH=. .venv/bin/python -m src.nonexact.run_experiment {E1,E2,E3,all}
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from itertools import combinations_with_replacement
from math import comb, gcd

import numpy as np
from flint import fmpz_poly, nmod_poly

from src.common.instances import balanced_semiprime
from src.nonexact import core as C

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")


def load_cfg():
    with open(os.path.join(ROOT, "experiments", "configs", "nonexact.json")) as f:
        return json.load(f)


def _write(name, rows):
    raw = os.environ.get("NONEXACT_RAW", RAW)
    os.makedirs(raw, exist_ok=True)
    with open(os.path.join(raw, f"nonexact_{name}.jsonl"), "w") as f:
        for r in rows:
            f.write(json.dumps(r, default=str) + "\n")


def _frac(x):
    return float(x)  # reporting only


# ------------------------------------------------------------------------------------------------ E1
def run_E1(cfg):
    c1 = cfg["E1"]
    rng = random.Random(c1["rng_seed"])
    rows, enum_rows = [], []
    for cell in c1["cells"]:
        n, V = cell["n"], cell["V"]
        N = balanced_semiprime(n, c1["instance_seed"]).N
        mu = rng.randrange(2, N - 1)
        for B in cell["B"]:
            for s in cell["s"]:
                model = C.lin_cost_model(N, B, s, V)
                h, p = model["h"], model["p"]
                T = math.ceil(c1["target_hits"] / h)
                Mtot = math.ceil(n * math.log(2) / math.log(2 * B)) + 2
                row = {"n": n, "N": N, "B": B, "s": s, "V": V, "M": Mtot, "T_needed": T,
                       "pred_p": _frac(p), "pred_h": _frac(h), "pred_cost_per_hit": _frac(model["cost_per_hit"]),
                       "brute_cost_per_hit": _frac(model["brute_cost_per_hit"]), "n_suffix": comb(2 * B + s, s)}
                if T > c1["T_max"]:
                    row["run"] = False
                    rows.append(row)
                    print("skip", row, flush=True)
                    continue
                t0 = time.time()
                ms = C.s_multisets_e(B, s)
                E = np.array([e for _, e in ms], dtype=np.int64)
                hits = 0
                nonunit = 0
                cs = []
                for _ in range(T):
                    pref = [rng.randint(-B, B) for _ in range(Mtot - s)]
                    c = 1
                    for a in pref:
                        c = c * (mu - a) % N
                    if gcd(c, N) != 1:
                        nonunit += 1  # would reveal a factor (ACCIDENTAL, prob ~ M B / p)
                        continue
                    w = C.lin_weights(c, mu, N, s)
                    if gcd(w[s], N) != 1:
                        nonunit += 1
                        continue
                    hits += C.lin_split_hits(w, N, V, E)
                    if len(cs) < max(c1["box_subset"], c1["enum_subset"]):
                        cs.append(c)
                # exact box points on a subset
                box = []
                for c in cs[: c1["box_subset"]]:
                    w = C.lin_weights(c, mu, N, s)
                    box.append(C.lin_box_points(w, N, B, s, V))
                pbar = sum(box) / len(box) if box else None
                p_used = max(1.0, pbar if (pbar is not None and _frac(p) >= 1) else _frac(p))
                meas_cost = (T * p_used) / hits if hits else None
                row.update({"run": True, "T": T, "hits": hits, "nonunit": nonunit,
                            "obs_over_pred_hits": hits / (T * _frac(h)),
                            "box_subset": box, "box_sum": sum(box), "box_pred_sum": len(box) * _frac(p),
                            "meas_cost_per_hit": meas_cost,
                            "meas_over_model_cost": (meas_cost / _frac(model["cost_per_hit"])) if meas_cost else None,
                            "brute_over_meas": (_frac(model["brute_cost_per_hit"]) / meas_cost) if meas_cost else None,
                            "secs": round(time.time() - t0, 2)})
                # E1b: LLL + exact Fincke-Pohst enumerator on a subset
                if _frac(p) <= c1["enum_p_max"]:
                    mism, nodes = 0, []
                    for c in cs[: c1["enum_subset"]]:
                        w = C.lin_weights(c, mu, N, s)
                        cnt = {}
                        pts = C.lin_lattice_enumerate(w, N, B, s, V, cnt)
                        exact = C.lin_box_points(w, N, B, s, V)
                        ok = len(pts) == exact and len({tuple(x) for x in pts}) == len(pts)
                        for x in pts:
                            val = (w[0] + sum(x[i] * w[i + 1] for i in range(s)) - x[s]) % N
                            ok = ok and val == 0
                        mism += 0 if ok else 1
                        nodes.append((cnt.get("nodes", 0), exact))
                    enum_rows.append({"n": n, "B": B, "s": s, "V": V, "lattices": len(nodes), "mismatches": mism,
                                      "nodes_points": nodes,
                                      "median_nodes_over_1p": float(np.median([a / (1 + b) for a, b in nodes]))})
                rows.append(row)
                print({k: row[k] for k in ("n", "B", "s", "T", "hits", "obs_over_pred_hits", "meas_over_model_cost",
                                           "secs")}, flush=True)
    _write("E1", rows)
    _write("E1b", enum_rows)


# ------------------------------------------------------------------------------------------------ E2
def _gh(n, d, N):
    # (N^d Gamma(n/2+1))^{1/n} / sqrt(pi)  -- reporting only
    return math.exp((d * math.log(N) + math.lgamma(n / 2 + 1)) / n) / math.sqrt(math.pi)


def run_E2(cfg):
    c2 = cfg["E2"]
    rng = random.Random(c2["rng_seed"])
    B, M = c2["B"], c2["M"]
    rows = []
    for n_bits in c2["n"]:
        N = balanced_semiprime(n_bits, c2["instance_seed"]).N
        for d in c2["d"]:
            mu, mu2 = rng.randrange(2, N - 1), rng.randrange(2, N - 1)
            F = fmpz_poly([(-mu) % N, 1]) if d == 1 else C.quadratic_presentation(mu, mu2, N)
            for (j, k) in (c2["jk_d1"] if d == 1 else c2["jk_d2"]):
                dim = j + 1 + k
                ratios, resid, useful, useful_degQ, qroots, base_roots, base_n, dets_ok = [], [], 0, 0, 0, 0, 0, 0
                for _ in range(c2["T"]):
                    A = [rng.randint(-B, B) for _ in range(M)]
                    P = C.from_roots(A)
                    Lb = C.multiplier_lattice(P, F, N, j, k)
                    dets_ok += int(C.det_abs(Lb) == N ** d)
                    Lb.sort(key=lambda v: sum(x * x for x in v))
                    b1 = Lb[0]
                    ratios.append(math.sqrt(sum(x * x for x in b1)) / _gh(dim, d, N))
                    rs = [v[j + 1:] for v in Lb if any(v[j + 1:])]
                    r1 = b1[j + 1:] if any(b1[j + 1:]) else min(rs, key=C.height)
                    resid.append(math.log2(max(1, C.height(r1))) / math.log2(N))
                    for v in Lb:
                        cl = C.classify_vector(v, j, k, B)
                        useful += int(cl["useful"])
                        useful_degQ += int(cl["useful"] and cl["Q_deg"] >= 1)
                        if cl["Q_deg"] >= 1:
                            qroots += int(cl["Q_int_roots"] > 0)
                            # baseline: random Q with the same coordinate heights
                            q = v[: j + 1]
                            for _b in range(c2["baseline_draws"]):
                                rq = [rng.randint(-abs(x), abs(x)) for x in q]
                                if any(rq[1:]):
                                    base_n += 1
                                    base_roots += int(C.int_roots_count(rq) > 0)
                nvec_degQ = None
                rows.append({"n": n_bits, "N": N, "d": d, "j": j, "k": k, "dim": dim, "T": c2["T"],
                             "det_ok": dets_ok, "median_b1_over_gh": float(np.median(ratios)),
                             "median_resid_exponent": float(np.median(resid)), "pred_resid_exponent": d / dim,
                             "useful": useful, "useful_degQ_ge1": useful_degQ, "Q_with_int_root": qroots,
                             "baseline_int_root_rate": base_roots / base_n if base_n else None,
                             "baseline_expected_count": (base_roots / base_n) * (c2["T"] * dim) if base_n else None})
                print(rows[-1], flush=True)
    _write("E2", rows)


# ------------------------------------------------------------------------------------------------ E3
def _fibre_pairs(k, M, R, cap, rng):
    """Pairs of M-multisets in [-R, R] with equal e_1..e_{M-k} and deg(P - P') = k - 1 exactly."""
    groups = {}
    for A in combinations_with_replacement(range(-R, R + 1), M):
        key = tuple(sum(a ** t for a in A) for t in range(1, M - k + 1))
        groups.setdefault(key, []).append(A)
    pairs = []
    for g in groups.values():
        if len(g) < 2:
            continue
        for i in range(len(g)):
            for jj in range(i + 1, len(g)):
                D = C.from_roots(g[i]) - C.from_roots(g[jj])
                if D.degree() == k - 1:
                    pairs.append((g[i], g[jj]))
    rng.shuffle(pairs)
    return pairs[:cap], len(pairs)


def _monic_factor_of_degree(f: nmod_poly, d: int):
    """A monic factor of f (over F_N, N prime) of degree exactly d, as a product of irreducible factors."""
    _, facs = f.factor()
    items = []
    for g, e in facs:
        items += [g] * e
    # subset-sum on degrees (small)
    best = {0: nmod_poly([1], f.modulus())}
    for g in items:
        new = dict(best)
        for dd, h in best.items():
            t = dd + g.degree()
            if t <= d and t not in new:
                new[t] = h * g
        best = new
    return best.get(d)


def run_E3(cfg):
    from sympy import isprime, primefactors
    c3 = cfg["E3"]
    rng = random.Random(c3["rng_seed"])
    H = c3["H"]
    rows = []
    for (k, M, R) in c3["cells"]:
        pairs, total_pairs = _fibre_pairs(k, M, R, c3["pairs_cap"], rng)
        d = max(k, 2)
        tested, viol, rel_counts, lb_counts, no_N, no_F = 0, 0, [], [], 0, 0
        examples = []
        for (A, A2) in pairs:
            P, P2 = C.from_roots(A), C.from_roots(A2)
            D = P - P2
            cD = C.coeffs(D, k)
            # candidate (N, delta): N prime > max(4H+1, N_min) dividing cD - delta coefficientwise, delta != cD
            cand = []
            import itertools
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
                no_N += 1
                continue
            N, delta = cand[rng.randrange(len(cand))]
            # r with ||r||, ||r - delta|| <= H, and a monic F of degree d dividing P - r mod N
            built = None
            for _ in range(c3["r_tries"]):
                r = [rng.randint(max(-H, x - H), min(H, x + H)) for x in delta]
                f = nmod_poly([int(x) % N for x in C.coeffs(P - fmpz_poly(r), M + 1)], N)
                g = _monic_factor_of_degree(f, d)
                if g is not None:
                    built = (r, fmpz_poly([int(x) for x in g.coeffs()]))
                    break
            if built is None:
                no_F += 1
                continue
            r, F = built
            r2 = [r[i] - delta[i] for i in range(k)]
            ok_mem = all((a - b) % N == 0 for a, b in zip(C.coeffs(C.red_mod(P, F), d), r + [0] * (d - k)))
            ok_mem2 = all((a - b) % N == 0 for a, b in zip(C.coeffs(C.red_mod(P2, F), d), r2 + [0] * (d - k)))
            rels, E = C.ih_relations(A, A2, r, r2, N, k)
            v_bound = 2 * H * sum(R ** l for l in range(k))
            bad = (not ok_mem) or (not ok_mem2)
            for (b, u, v, nt) in rels:
                bad = bad or ((u - v) % N != 0) or (nt != (int(E(b)) != 0)) or (nt != (u != v)) or abs(v) > v_bound
                bad = bad or any(abs(b - a) > 2 * R for a in (A if b in A2 else A2))
            # sign-normalise: (u, v) and (-u, -v) are the same relation
            distinct = {(u, v) if u > 0 else (-u, -v) for (_, u, v, nt) in rels if nt}
            n_supp = len(set(A) ^ set(A2))
            lb = C.ih_distinct_lower_bound(M, k, D.degree(), n_supp)
            bad = bad or len(distinct) < lb
            if k == 1:
                bad = bad or len(distinct) != 1
            viol += int(bad)
            tested += 1
            rel_counts.append(len(distinct))
            lb_counts.append(lb)
            if len(examples) < 3:
                examples.append({"A": A, "A2": A2, "N": N, "F": C.coeffs(F), "r": r, "r2": r2,
                                 "relations": sorted(distinct)[:6]})
        rows.append({"k": k, "M": M, "R": R, "pairs_total": total_pairs, "pairs_used": len(pairs),
                     "tested": tested, "violations": viol, "no_N": no_N, "no_F": no_F,
                     "mean_distinct_rel": (sum(rel_counts) / len(rel_counts)) if rel_counts else None,
                     "min_distinct_rel": min(rel_counts) if rel_counts else None,
                     "mean_lower_bound": (sum(lb_counts) / len(lb_counts)) if lb_counts else None,
                     "examples": examples})
        print({x: rows[-1][x] for x in ("k", "M", "R", "pairs_total", "tested", "violations", "mean_distinct_rel",
                                         "mean_lower_bound")}, flush=True)
    _write("E3", rows)


def main(argv):
    cfg = load_cfg()
    which = argv[1] if len(argv) > 1 else "all"
    if which in ("E1", "all"):
        run_E1(cfg)
    if which in ("E2", "all"):
        run_E2(cfg)
    if which in ("E3", "all"):
        run_E3(cfg)


if __name__ == "__main__":
    main(sys.argv)
