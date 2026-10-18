"""Offline analysis for pm1bsgs (HARNESS side).

Reads experiments/raw/pm1bsgs_{phase1,smooth,rho,validate}.jsonl*, writes
experiments/analysis/pm1bsgs_*.json and prints the tables used in
notes/pm1bsgs_analysis.md.

  .venv/bin/python -m src.pm1bsgs.analysis all
"""
from __future__ import annotations

import gzip
import json
import math
import os
import sys
from multiprocessing import Pool
from typing import Dict, List, Tuple

import numpy as np

from src.common.instances import is_prime

from . import costmodel as cm
from . import harness as H
from . import operational as op
from .experiments import BDAY_B1, BDAY_KMAX, B1_GRID, NBASES, NS, RAW, ROOT

AN = os.path.join(ROOT, "experiments", "analysis")
C_LO, C_HI = 6, 36
RESTART_C0 = 8


def _r(path):
    op_ = gzip.open if path.endswith(".gz") else open
    with op_(path, "rt") as f:
        return [json.loads(l) for l in f]


def _dump(name, obj):
    os.makedirs(AN, exist_ok=True)
    with open(os.path.join(AN, f"pm1bsgs_{name}.json"), "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True)


# ------------------------------------------------------------------ parameter grids and costs
def powN_nominal(n: int) -> int:
    return (n - 1) + (n // 2 - 1)


def grid(method: str, n: int) -> List[Tuple[Dict, int]]:
    out = []
    if method in ("a_N", "a_order", "d_plain", "a_shift"):
        mode = "order" if method == "a_order" else "N"
        b1s = [1] if method == "d_plain" else B1_GRID
        for B1 in b1s:
            s1 = 1 + cm.stage1_cost(B1)
            for s in range(0, 22):
                for J in (1 << s, 1 << (s + 1)):
                    mb = 1 << s
                    c = s1 + cm.bsgs_stage2_cost(mb, J, "order", 0) + (powN_nominal(n) if mode == "N" else 0)
                    out.append(({"B1": B1, "m": mb, "J": J, "mode": mode}, c))
    elif method == "a_N_gp":
        # same algorithm as a_N, but charged under the geometric-progression model:
        # build F with GP roots in 2 M(m), evaluate at GP points in 3 M(m) per block
        # (Bostan-Schost / chirp-z style; constants ASSUMED, not implemented here)
        for B1 in B1_GRID:
            s1 = 1 + cm.stage1_cost(B1)
            for s in range(0, 22):
                for J in (1 << s, 1 << (s + 1)):
                    mb = 1 << s
                    nb = -(-J // mb)
                    c = (s1 + powN_nominal(n) + (mb - 1) + op.pow_cost(mb) + (J - 1)
                         + 2 * op.mcost(mb) + nb * (3 * op.mcost(mb) + 1) + J)
                    out.append(({"B1": B1, "m": mb, "J": J, "mode": "N"}, c))
    elif method == "b_prime":
        for B1 in B1_GRID[1:]:
            s1 = 1 + cm.stage1_cost(B1)
            out.append(({"B1": B1, "B2": B1}, s1))
            k1 = B1.bit_length() - 1
            for k2 in range(k1 + 1, 31):
                B2 = 1 << k2
                out.append(({"B1": B1, "B2": B2}, s1 + cm.prime_stage2_cost(B1, B2)))
    elif method == "c_bday":
        for B1 in BDAY_B1:
            s1 = 1 + cm.stage1_cost(B1)
            for s in range(1, 22):
                K = 1 << s
                out.append(({"B1": B1, "K": K}, s1 + cm.birthday_stage2_cost(K)))
    else:
        raise ValueError(method)
    return out


METHODS = ("a_N", "a_order", "a_shift", "a_N_gp", "b_prime", "c_bday", "d_plain")


# ------------------------------------------------------------------ outcomes
class Inst:
    def __init__(self, rec):
        self.n = rec["n"]
        self.seed = rec["seed"]
        self.N = int(rec["N"])
        self.p, self.q = int(rec["p"]), int(rec["q"])
        self.P = H.PrimeInfo(self.p, {l: v for l, v in rec["pfac"]})
        self.Q = H.PrimeInfo(self.q, {l: v for l, v in rec["qfac"]})
        self.bases = [(int(a), int(b)) for a, b in rec["bases"]]
        self.bday0 = rec["bday_first_pure_base0"]
        self._st = {}

    def stage1(self, bseed: int, B1: int):
        key = (bseed, B1)
        if key not in self._st:
            op_, oq = self.bases[bseed]
            bs = H.BaseState(0, op_, oq, H.factor_of_order(op_, self.P.fac) if op_ else {},
                             H.factor_of_order(oq, self.Q.fac) if oq else {})
            self._st[key] = H.stage1_outcome(bs, B1)
        return self._st[key]


def success(inst: Inst, bseed: int, method: str, params: Dict, bday_cache=None) -> bool:
    B1 = params["B1"]
    st = inst.stage1(bseed, B1)
    if st[0] != "ok":
        return st[0] == "stage1"
    dp, dq = st[1], st[2]
    if method in ("a_N", "a_order", "a_N_gp", "d_plain"):
        return H.predict_bsgs(dp, dq, params["m"], params["J"], params["mode"], inst.N)
    if method == "a_shift":
        T = op.shift_target(inst.N, bseed)
        lo, hi = H.window(params["m"], params["J"], "N", T)
        return H.interval_pure_hit(dp, dq, lo, hi)
    if method == "b_prime":
        if params["B2"] <= B1:
            return False
        return H.predict_prime(dp, dq, B1, params["B2"])
    if method == "c_bday":
        if bseed == 0:
            fp = inst.bday0[str(B1)]
        else:
            fp = bday_cache[(inst.seed, bseed, B1)]
        return isinstance(fp, int) and fp <= params["K"]
    raise ValueError(method)


GRIDS: Dict = {}


def _succ_matrix_worker(args):
    rec, method = args
    inst = Inst(rec)
    g = GRIDS[(method, inst.n)]
    return [success(inst, 0, method, p) for p, _ in g]


# ------------------------------------------------------------------ main analysis
def budget_and_restarts(recs_by_n: Dict[int, List[Dict]]):
    out = {}
    for n in NS:
        for method in METHODS:
            GRIDS[(method, n)] = grid(method, n)
    pool = Pool(32)
    for n in NS:
        recs = recs_by_n[n]
        half = len(recs) // 2
        train, hold = recs[:half], recs[half:]
        out[n] = {}
        for method in METHODS:
            g = GRIDS[(method, n)]
            costs = np.array([c for _, c in g], dtype=np.float64)
            S = np.array(pool.map(_succ_matrix_worker, [(r, method) for r in recs], chunksize=8), dtype=bool)
            Str, Sho = S[:half], S[half:]
            mtr = Str.mean(0)
            curve = []
            best_by_c = {}
            for c in range(C_LO, C_HI + 1):
                ok = costs <= 2.0 ** c
                if not ok.any():
                    curve.append({"c": c, "param": None, "train": 0.0, "holdout": 0.0})
                    continue
                idx = np.where(ok)[0]
                # tie-break on train success: prefer the LARGER cost (more coverage), so a
                # saturated schedule keeps growing its window on restarts
                j = idx[np.lexsort((-costs[idx], -mtr[idx]))[0]]
                best_by_c[c] = int(j)
                curve.append({"c": c, "param": g[j][0], "cost": float(costs[j]), "train": float(mtr[j]),
                              "holdout": float(Sho[:, j].mean())})
            # restarts on holdout: level c uses best_by_c[c] with base seed (c - C0 + 1)
            jobs = [(r, method, best_by_c, [g[best_by_c[c]][0] if c in best_by_c else None for c in range(C_LO, C_HI + 1)],
                     costs.tolist()) for r in hold]
            res = pool.map(_restart_worker, jobs, chunksize=4)
            tot = np.array([x[0] for x in res])
            cens = np.array([x[1] for x in res])
            out[n][method] = {"curve": curve, "restart_costs": tot.tolist(), "censored": cens.tolist(),
                              "mean_cost": float(tot.mean()), "median_cost": float(np.median(tot)),
                              "frac_censored": float(cens.mean()),
                              "levels_used": [x[2] for x in res]}
            print(n, method, f"mean={tot.mean():.3g} med={np.median(tot):.3g} cens={cens.mean():.3f}", flush=True)
    pool.close()
    return out


def _restart_worker(args):
    rec, method, best_by_c, params_by_c, costs = args
    inst = Inst(rec)
    total = 0.0
    level = 0
    for c in range(RESTART_C0, C_HI + 1):
        if c not in best_by_c:
            continue
        j = best_by_c[c]
        params = params_by_c[c - C_LO]
        bseed = 1 + (c - RESTART_C0)
        total += costs[j]
        level += 1
        cache = None
        if method == "c_bday":
            st = inst.stage1(bseed, params["B1"])
            cache = {}
            if st[0] == "ok":
                fp = H.birthday_first_pure(inst.N, bseed, st[1], st[2], params["K"])
                cache[(inst.seed, bseed, params["B1"])] = fp if fp is not None else "none"
            else:
                cache[(inst.seed, bseed, params["B1"])] = st[0]
        if success(inst, bseed, method, params, cache):
            return total, False, level
    return total, True, level


def fit_slope(ns, ys):
    x = np.array(ns, dtype=float)
    y = np.log2(np.array(ys, dtype=float))
    A = np.vstack([x, np.ones_like(x)]).T
    sl, ic = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(sl), float(ic)


def rho_stats():
    recs = _r(os.path.join(RAW, "pm1bsgs_rho.jsonl"))
    out = {}
    for n in NS:
        xs = [r["eq_ops"] for r in recs if r["n"] == n]
        out[n] = {"mean": float(np.mean(xs)), "median": float(np.median(xs)), "count": len(xs),
                  "found": sum(r["found"] for r in recs if r["n"] == n)}
    return out


# ------------------------------------------------------------------ smooth-part statistics
def dickman_rho(umax=40.0, h=1e-3):
    n = int(umax / h) + 1
    u = np.arange(n) * h
    rho = np.ones(n)
    k1 = int(round(1 / h))
    for i in range(k1 + 1, n):
        # rho'(u) = -rho(u-1)/u, trapezoid
        rho[i] = rho[i - 1] - h * 0.5 * (rho[i - 1 - k1] / u[i - 1] + rho[i - k1] / u[i])
    return u, rho


def smooth_stats():
    recs = _r(os.path.join(RAW, "pm1bsgs_smooth.jsonl.gz"))
    u, rho = dickman_rho()
    h = u[1] - u[0]
    cdf_model = np.exp(-0.5772156649015329) * np.cumsum(rho) * h  # P(u_S <= x)
    out = {"by_n": {}, "model_cdf_points": {}}
    for x in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        out["model_cdf_points"][x] = float(np.interp(x, u, cdf_model))
    for n in NS:
        rs = [r for r in recs if r["n"] == n]
        row = {"count": len(rs)}
        for kB in (8, 12, 16):
            B = 1 << kB
            us = []
            for r in rs:
                S = 1
                for l, v in r["fac"]:
                    S *= l ** min(v, H.e_l(l, B))
                us.append(math.log(S) / math.log(B))
            us = np.array(us)
            row[f"B=2^{kB}"] = {"mean_u": float(us.mean()), "median_u": float(np.median(us)),
                                "cdf": {x: float((us <= x).mean()) for x in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)},
                                "E_S^-1/2": float(np.mean([B ** (-x / 2) for x in us]))}
        # Sophie-Germain-type: (r-1)/2 prime ; largest prime factor statistics
        sg = np.mean([len(r["fac"]) == 2 and r["fac"][0] == [2, 1] for r in rs])
        lpf = [math.log(max(l for l, _ in r["fac"])) / (r["bits"] * math.log(2)) for r in rs]
        row["frac_(r-1)/2_prime"] = float(sg)
        row["HL_prediction_C2/ln r"] = float(0.6601618158 / ((n // 2) * math.log(2)))
        row["median_log_P1/log_r"] = float(np.median(lpf))
        row["frac_S_2^16_le_2^4"] = float(np.mean([_S(r, 1 << 16) <= 16 for r in rs]))
        out["by_n"][n] = row
    return out


def _S(r, B):
    S = 1
    for l, v in r["fac"]:
        S *= l ** min(v, H.e_l(l, B))
    return S


def validate_stats():
    recs = _r(os.path.join(RAW, "pm1bsgs_validate.jsonl"))
    out = {"runs": len(recs), "disagree": sum(not r["agree"] for r in recs)}
    for m in ("bsgs", "prime", "bday"):
        rs = [r for r in recs if r["method"] == m]
        out[m] = {"runs": len(rs), "found": sum(r["found"] for r in rs),
                  "structural": sum(r["cls"] == "STRUCTURAL" for r in rs),
                  "accidental": sum(r["cls"] == "ACCIDENTAL" for r in rs),
                  "pred_common2": sum(r["predicted"] == "common2" for r in rs),
                  "pred_stage1_N": sum(r["predicted"] == "stage1_N" for r in rs),
                  "backtracks": sum(1 for r in rs if r["counts"].get("backtrack_gcd", 0) > 0)}
    return out


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("smooth", "all"):
        s = smooth_stats()
        _dump("smooth", s)
        print(json.dumps(s, indent=1)[:4000])
    if what in ("asym", "all"):
        a = asymptotics()
        a["empirical_ES"] = empirical_ES()
        _dump("asymptotics", a)
        print(json.dumps(a, indent=1))
    if what in ("basedep", "all"):
        b = base_dependence()
        _dump("basedep", b)
        print(json.dumps(b, indent=1))
    if what in ("fail", "all"):
        _dump("failure", failure_stats())
    if what in ("validate", "all"):
        v = validate_stats()
        _dump("validate", v)
        print(v)
    if what in ("budget", "all"):
        recs = _r(os.path.join(RAW, "pm1bsgs_phase1.jsonl.gz"))
        by_n = {n: sorted([r for r in recs if r["n"] == n], key=lambda r: r["seed"]) for n in NS}
        res = budget_and_restarts(by_n)
        rho = rho_stats()
        fits = {}
        for m in METHODS:
            fits[m] = {"mean": fit_slope(NS, [res[n][m]["mean_cost"] for n in NS]),
                       "median": fit_slope(NS, [res[n][m]["median_cost"] for n in NS]),
                       "mean_n>=56": fit_slope(NS[3:], [res[n][m]["mean_cost"] for n in NS[3:]])}
        fits["rho"] = {"mean": fit_slope(NS, [rho[n]["mean"] for n in NS]),
                       "median": fit_slope(NS, [rho[n]["median"] for n in NS]),
                       "mean_n>=56": fit_slope(NS[3:], [rho[n]["mean"] for n in NS[3:]])}
        _dump("budget", {str(n): {m: {k: v for k, v in res[n][m].items() if k not in ("restart_costs", "censored", "levels_used")}
                                  for m in METHODS} for n in NS})
        _dump("restarts", {"per_n": {str(n): {m: {"mean": res[n][m]["mean_cost"], "median": res[n][m]["median_cost"],
                                                  "frac_censored": res[n][m]["frac_censored"],
                                                  "costs": res[n][m]["restart_costs"]} for m in METHODS} for n in NS},
                           "rho": {str(n): rho[n] for n in NS}, "fits_log2cost_vs_n": fits})
        print(json.dumps(fits, indent=1))



# ------------------------------------------------------------------ asymptotics (floating point: model evaluation only)
def buchstab_omega(umax=40.0, h=1e-3):
    n = int(umax / h) + 1
    u = np.arange(n) * h
    w = np.zeros(n)
    k1 = int(round(1 / h))
    k2 = 2 * k1
    w[k1:k2 + 1] = 1.0 / u[k1:k2 + 1]
    # (u w(u))' = w(u-1)
    uw = u * w
    for i in range(k2 + 1, n):
        uw[i] = uw[i - 1] + h * 0.5 * (w[i - 1 - k1] + w[i - k1])
        w[i] = uw[i] / u[i]
    return u, w


def asymptotics():
    u, rho = dickman_rho(60.0)
    _, om = buchstab_omega(60.0)
    h = u[1] - u[0]
    gam = 0.5772156649015329

    def rho_at(x):
        return float(np.interp(x, u, rho))

    def G(w, v):
        """P(log S_B / log B >= v) for an integer of size B^w (rho*omega model,
        atom rho(w) at full smoothness)."""
        i0, i1 = int(round(max(v, 0) / h)), int(round((w - 1) / h))
        mass = 0.0
        if i1 > i0:
            j = np.arange(i0, i1 + 1)
            vals = rho[j] * np.interp(w - u[j], u, om)
            mass = float(np.trapezoid(vals, dx=h))
        return mass + (rho_at(w) if v <= w else 0.0)

    out = {"mass_check": {w: G(w, 0.0) for w in (3.0, 5.0, 10.0, 20.0)}}
    # c(delta): cost exponent (in log p) at success probability delta, cost max(B, sqrt X)
    # optimum is B = sqrt X (so the window B^2 only admits R = 1 or a prime): G(1/c, 1/c - 2)
    def c_for(delta_p):
        lo, hi = 0.02, 0.5
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            w = 1.0 / mid
            if G(w, w - 2.0) >= delta_p:
                hi = mid
            else:
                lo = mid
        return hi
    cd = {}
    for d in (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99):
        cd[d] = {"single_prime_c": c_for(d), "semiprime_c": c_for(1 - math.sqrt(1 - d))}
    out["c_of_delta"] = cd

    # prime-only stage 2 at balanced cost: success ~ p-1 fully B-smooth -> rho(1/c) (log factors ignored)
    def c_prime(dp):
        lo, hi = 0.02, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if rho_at(1.0 / mid) >= dp:
                hi = mid
            else:
                lo = mid
        return hi
    out["c_of_delta_prime_only"] = {d: {"single_prime_c": c_prime(d), "semiprime_c": c_prime(1 - math.sqrt(1 - d))}
                                    for d in (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)}
    # general (beta, xi): show optimum has beta = xi/2 (scan)
    scan = []
    for beta in np.linspace(0.05, 0.5, 46):
        for xi in np.linspace(0.0, 1.0, 101):
            w = 1.0 / beta
            v = (1.0 - xi) / beta
            scan.append((float(beta), float(xi), max(float(beta), float(xi) / 2), G(w, v)))
    best = {}
    for d in (0.1, 0.5, 0.9):
        feas = [s for s in scan if s[3] >= d]
        b = min(feas, key=lambda s: (s[2], -s[3]))
        best[d] = {"beta": b[0], "xi": b[1], "cost_exp": b[2], "succ": b[3]}
    out["scan_best_single_prime"] = best
    # E[S_B^{-1/2}] for p-1, p prime (independent valuations, P(l^k | p-1) = 1/phi(l^k))
    from .operational import primes_upto
    es = {}
    for kB in (4, 8, 12, 16, 20, 24, 30, 40, 60):
        B = 2 ** kB
        logprod = 0.0
        lim = min(B, 1 << 22)
        for l in primes_upto(lim):
            e = H.e_l(l, B)
            tot = 0.0
            # P(min(v,e) = k)
            for k in range(0, e + 1):
                pk = (1.0 if k == 0 else 1.0 / ((l - 1) * l ** (k - 1)))
                pk1 = (1.0 / ((l - 1) * l ** k)) if k < e else 0.0
                tot += (pk - pk1) * l ** (-k / 2)
            logprod += math.log(tot)
        if B > lim:  # primes in (2^22, B]: e = 1, factor (1 - 1/(l-1)(1 - l^{-1/2})); Mertens tail
            logprod += -(math.log(math.log(B)) - math.log(math.log(lim)))  # sum 1/l over (lim,B]
        es[kB] = {"E_S^-1/2": math.exp(logprod), "times_lnB": math.exp(logprod) * kB * math.log(2)}
    out["E_S_minus_half_model"] = es
    return out


def empirical_ES():
    """E[S^{-1/2}] (single prime) and E[max(S_p,S_q)^{-1/2}] (semiprime) from phase-1 factorisations."""
    recs = _r(os.path.join(RAW, "pm1bsgs_phase1.jsonl.gz"))
    out = {}
    for n in NS:
        rs = [r for r in recs if r["n"] == n]
        row = {}
        for kB in (4, 8, 12, 16, 20):
            B = 1 << kB
            sp = np.array([_S({"fac": r["pfac"]}, B) for r in rs], dtype=float)
            sq = np.array([_S({"fac": r["qfac"]}, B) for r in rs], dtype=float)
            row[kB] = {"E_S^-1/2": float(np.mean(sp ** -0.5)),
                       "E_max^-1/2": float(np.mean(np.maximum(sp, sq) ** -0.5)),
                       "P_S=2": float(np.mean(sp == 2))}
        out[n] = row
    return out



def failure_stats(per_n: int = 200):
    """Exact gcd-failure probabilities (over a uniform base) on phase-1 instances:
    user's single test at x* = N mod d_p; P(d_p = d_q); and the full-pipeline
    outcome distribution for B1 = 2^10 with a window covering all of R."""
    recs = _r(os.path.join(RAW, "pm1bsgs_phase1.jsonl.gz"))
    out = {}
    for n in NS:
        rs = sorted([r for r in recs if r["n"] == n], key=lambda r: r["seed"])[:per_n]
        ux, eq, s1N, com, succ, nohit, used = [], [], [], [], [], [], 0
        for r in rs:
            inst = Inst(r)
            B1 = 1 << 10
            u_ = H.user_first_x_failure_prob(inst.P, inst.Q, inst.N, B1)
            e_ = H.prob_equal_orders(inst.P, inst.Q, B1, exclude_one=True)
            if u_ is not None:
                ux.append(float(u_))
                eq.append(float(e_))
            X = max(inst.p, inst.q)
            s = (X.bit_length() + 1) // 2
            params = {"B1": B1, "m": 1 << s, "J": 1 << s, "mode": "N"}
            ex = H.exact_outcome_probs(inst.P, inst.Q, inst.N, "bsgs", params)
            if ex is None:
                continue
            used += 1
            s1N.append(float(ex.get("stage1_N", 0)))
            com.append(float(ex.get("common2", 0)))
            succ.append(float(ex.get("stage1", 0) + ex.get("success2", 0)))
            nohit.append(float(ex.get("nohit", 0)))
        out[n] = {"instances": len(rs), "user_first_x_fail_mean": float(np.mean(ux)),
                  "user_first_x_fail_max": float(np.max(ux)), "user_first_x_n_used": len(ux),
                  "user_first_x_fail_median": float(np.median(ux)),
                  "P_equal_orders_mean": float(np.mean(eq)), "P_equal_orders_max": float(np.max(eq)),
                  "pipeline_used": used, "stage1_N_mean": float(np.mean(s1N)), "common2_mean": float(np.mean(com)),
                  "success_mean": float(np.mean(succ)), "nohit_mean": float(np.mean(nohit)),
                  "min_success": float(np.min(succ))}
        print(n, out[n], flush=True)
    return out



def base_dependence():
    """How much does a fresh base change the outcome at a fixed parameter?
    For a_N and b_prime at the train-best parameter for budgets 2^18 and 2^22:
    fraction of holdout instances whose success over 10 bases is all-or-nothing."""
    recs = _r(os.path.join(RAW, "pm1bsgs_phase1.jsonl.gz"))
    bud = json.load(open(os.path.join(AN, "pm1bsgs_budget.json")))
    out = {}
    for n in (48, 64, 80):
        rs = sorted([r for r in recs if r["n"] == n], key=lambda r: r["seed"])
        hold = rs[len(rs) // 2:]
        out[n] = {}
        for method in ("a_N", "b_prime"):
            for c in (18, 22):
                params = [x for x in bud[str(n)][method]["curve"] if x["c"] == c][0]["param"]
                allsame, rate = 0, []
                for r in hold:
                    inst = Inst(r)
                    v = [success(inst, b, method, params) for b in range(10)]
                    allsame += (all(v) or not any(v))
                    rate.append(np.mean(v))
                out[n][f"{method}@2^{c}"] = {"param": params, "frac_all_or_nothing": allsame / len(hold),
                                            "mean_success": float(np.mean(rate))}
    return out


if __name__ == "__main__":
    main()
