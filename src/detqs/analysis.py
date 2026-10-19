"""HARNESS analysis for E1-E3 (frozen with the pre-registration).  Reads
experiments/raw/detqs_e*.  Writes experiments/analysis/detqs_summary.json."""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from typing import Dict, List

import numpy as np

from .harness import ANA, RAW, E1_M, E2_M, T_of, y_of

# model local pair ratios rho_l^model = avg_{l|k} R / avg_{l nmid k} R (frozen, from
# localmodel.pair_ratio_model at u(T/4), see prereg table)
MODEL_R = {}  # filled from prereg values at import of summarize


def _binom_two_sided(k: int, M: int, p: float) -> float:
    """Exact two-sided binomial p-value (sum of probs <= P(k)); log-space pmf
    (post-prereg fix: the float form overflowed at M = 4600; same test)."""
    def lp(j: int) -> float:
        return (math.lgamma(M + 1) - math.lgamma(j + 1) - math.lgamma(M - j + 1)
                + j * math.log(p) + (M - j) * math.log1p(-p))
    lk = lp(k)
    tot = 0.0
    for j in range(M + 1):
        lj = lp(j)
        if lj <= lk + 1e-9:
            tot += math.exp(lj)
    return min(1.0, tot)


def _var(xs: List[float]) -> float:
    m = sum(xs) / len(xs)
    return sum((x - m) ** 2 for x in xs) / len(xs)


def e1_stats(n: int) -> Dict:
    recs = [json.loads(l) for l in open(os.path.join(RAW, f"detqs_e1_n{n}.jsonl"))]
    out: Dict = {"n": n, "y": y_of(n), "T": T_of(n)}
    I = [r for r in recs if r["kind"] == "int" and not r.get("square")]
    S = [r for r in recs if r["kind"] == "sp" and not r.get("square")]
    Zi = [r["Z"] for r in I]
    Li = [r["L"] for r in I]
    Zbar = sum(Zi) / len(Zi)
    Lbar = sum(Li) / len(Li)
    V = _var([z / Zbar for z in Zi])
    VL = _var([l / Lbar for l in Li])
    corr = float(np.corrcoef(np.array(Zi, float), np.array(Li, float))[0, 1])
    c = sum(Zi) / sum(Li)
    ri = [z / (c * l) for z, l in zip(Zi, Li)]
    Vres = _var(ri)
    strata = defaultdict(list)
    for r in I:
        N = int(r["N"])
        strata[(N % 8, N % 3, N % 5)].append(r["Z"] / Zbar)
    num, den = 0.0, 0
    for v in strata.values():
        if len(v) >= 2:
            m = sum(v) / len(v)
            num += sum((x - m) ** 2 for x in v)
            den += len(v)
    Vstrat = num / den
    rs = [r["Z"] / (c * r["L"]) for r in S]
    out.update(dict(M_int=len(I), M_sp=len(S), Zbar=Zbar, V=V, poisson=1 / Zbar, VL=VL, corr=corr,
                    Vres=Vres, Vstrat=Vstrat, min_r_int=min(ri), min_r_sp=min(rs),
                    mean_r_sp=sum(rs) / len(rs), min_Z_over_fb_sp=min(r["Z"] / (r["fb"] + 2) for r in S),
                    h1_fail_sp=sum(1 for r in S if r["Z"] < r["fb"] + 2),
                    frac_Z_below_half_mean_int=sum(1 for z in Zi if z < Zbar / 2) / len(Zi)))
    out["P1a"] = V >= 10 / Zbar
    out["P1b"] = 0.4 * VL <= V <= 2.5 * VL
    out["P1c"] = corr >= 0.8
    out["P1d"] = Vres <= 0.25 * V
    out["P1e"] = Vstrat <= 0.6 * V
    out["P1f"] = min(ri) >= 0.5 and min(rs) >= 0.5
    out["P1g"] = out["h1_fail_sp"] == 0
    out["P1h"] = abs(out["mean_r_sp"] - 1) <= 0.15
    return out


def e3_stats(n: int, model: Dict[int, float]) -> Dict:
    d = np.load(os.path.join(RAW, f"detqs_e3_n{n}.npz"))
    A, G, cnt, T = d["A"], d["G"], int(d["cnt"]), int(d["T"])
    ks = np.arange(T // 8, T // 2)
    C = (A[ks] / cnt) / (G[ks] / cnt / cnt)
    out = {"n": n}
    for l in (2, 3, 5):
        a = C[ks % l == 0].mean()
        b = C[ks % l != 0].mean()
        rat = float(a / b)
        out[f"ratio_{l}"] = rat
        out[f"model_{l}"] = model[l]
        out[f"P3a_{l}"] = rat > 1
        out[f"P3b_{l}"] = 0.4 <= (rat - 1) / (model[l] - 1) <= 2.0
    return out


def e2_stats(n: int) -> Dict:
    recs = [json.loads(l) for l in open(os.path.join(RAW, f"detqs_e2_n{n}.jsonl"))]
    M = len(recs)
    dstar = [r["first_fruitful"] if r["first_fruitful"] is not None else 10 ** 6 for r in recs]
    out: Dict = {"n": n, "M": M}
    pv = {}
    for j in range(1, 6):
        k = sum(1 for d in dstar if d > j)
        pv[j] = (k, _binom_two_sided(k, M, 2.0 ** (-j)))
    out["tail_counts"] = {j: pv[j][0] for j in pv}
    out["tail_pvals"] = {j: pv[j][1] for j in pv}
    out["n_gt15"] = sum(1 for d in dstar if d > 15)
    jm = [r["jacobi_match"] for r in recs if r["jacobi_match"] is not None]
    out["jacobi_runs"] = len(jm)
    out["P2a"] = all(jm)
    out["P2b_count"] = sum(1 for j in pv if pv[j][1] >= 0.001)
    out["h1_fail"] = sum(1 for r in recs if r["Z"] < r["fb"] + 2)
    out["accidental"] = sum(1 for r in recs if r["class"] == "ACCIDENTAL")
    out["success"] = sum(1 for r in recs if r["class"] == "STRUCTURAL")
    small = [(s, o) for r in recs for s, o in zip(r["dep_sizes"], r["outcomes"]) if s <= 2]
    out["small_deps"] = len(small)
    out["P2f"] = all(o == 1 for s, o in small)
    out["mean_deps"] = sum(r["deps_tested"] for r in recs) / M
    return out


def summarize() -> None:
    from .prereg_values import MODEL_PAIR_RATIO
    res = {"E1": {}, "E2": {}, "E3": {}}
    for n in sorted(E1_M):
        if os.path.exists(os.path.join(RAW, f"detqs_e1_n{n}.jsonl")):
            res["E1"][n] = e1_stats(n)
            res["E3"][n] = e3_stats(n, MODEL_PAIR_RATIO[n])
    for n in sorted(E2_M):
        if os.path.exists(os.path.join(RAW, f"detqs_e2_n{n}.jsonl")):
            res["E2"][n] = e2_stats(n)
    # pooled E2
    allrec = []
    for n in res["E2"]:
        allrec += [json.loads(l) for l in open(os.path.join(RAW, f"detqs_e2_n{n}.jsonl"))]
    if allrec:
        M = len(allrec)
        dstar = [r["first_fruitful"] if r["first_fruitful"] is not None else 10 ** 6 for r in allrec]
        res["E2_pooled"] = {"M": M, "tail": {j: (sum(1 for d in dstar if d > j),
                                                _binom_two_sided(sum(1 for d in dstar if d > j), M, 2.0 ** (-j)))
                                             for j in range(1, 6)},
                            "n_gt15": sum(1 for d in dstar if d > 15)}
        res["E2_pooled"]["P2b_pooled"] = all(v[1] >= 0.001 for v in res["E2_pooled"]["tail"].values())
        res["E2_pooled"]["P2c"] = res["E2_pooled"]["n_gt15"] <= 2
    json.dump(res, open(os.path.join(ANA, "detqs_summary.json"), "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))
