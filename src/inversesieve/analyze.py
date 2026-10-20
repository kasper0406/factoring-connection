"""Frozen analysis for notes/inversesieve_prereg.md (verdicts P1-P9).  Floats are used only for summary
statistics of exact counts.

Run:  .venv/bin/python -m src.inversesieve.analyze  -> experiments/analysis/inversesieve_summary.json
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
OUT = os.path.join(ROOT, "experiments", "analysis", "inversesieve_summary.json")


def _load(name):
    p = os.path.join(RAW, f"inversesieve_{name}.jsonl")
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


def ks_pvalue(D, n_eff):
    """Asymptotic Kolmogorov p-value Q(sqrt(n_eff) D) (two-sided)."""
    lam = (math.sqrt(n_eff) + 0.12 + 0.11 / math.sqrt(n_eff)) * D
    if lam < 1e-9:
        return 1.0
    s = sum(2 * (-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 200))
    return max(0.0, min(1.0, s))


def ks_2samp(a, b):
    a, b = np.sort(np.asarray(a, float)), np.sort(np.asarray(b, float))
    allv = np.concatenate([a, b])
    cdfa = np.searchsorted(a, allv, side="right") / len(a)
    cdfb = np.searchsorted(b, allv, side="right") / len(b)
    D = float(np.max(np.abs(cdfa - cdfb)))
    return D, ks_pvalue(D, len(a) * len(b) / (len(a) + len(b)))


def ks_uniform(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if n == 0:
        return None, None
    i = np.arange(1, n + 1)
    D = float(max(np.max(i / n - x), np.max(x - (i - 1) / n)))
    return D, ks_pvalue(D, n)


def rr_ci(a, n1, b, n2):
    """Risk ratio (a/n1)/(b/n2) with a 95% log-normal CI (0.5 continuity correction)."""
    a_, b_ = a + 0.5, b + 0.5
    rr = (a_ / (n1 + 0.5)) / (b_ / (n2 + 0.5))
    se = math.sqrt(1 / a_ - 1 / (n1 + 0.5) + 1 / b_ - 1 / (n2 + 0.5))
    return rr, rr * math.exp(-1.96 * se), rr * math.exp(1.96 * se)


def within(x, lo, hi):
    return x is not None and lo <= x <= hi


def analyze_e1(recs):
    by = defaultdict(list)
    for r in recs:
        by[(r["n"], r["shape"])].append(r)
    out = {}
    for key, rs in sorted(by.items()):
        c = {}
        c["N_count"] = len(rs)
        c["U_mean"] = float(np.mean([r["U_size"] for r in rs]))
        c["U_var_over_mean"] = float(np.var([r["U_size"] for r in rs]) / max(1e-9, c["U_mean"]))
        c["U_empty"] = sum(1 for r in rs if r["U_size"] == 0)
        c["P1_ratio"] = sum(r["U_size"] for r in rs) / sum(r["E_model"] for r in rs)
        c["P1"] = within(c["P1_ratio"], 0.85, 1.15)
        for tag in ("q", "1", "4"):
            obs = np.array([r["cover"][tag]["obs"] for r in rs], float)
            nm = np.array([np.mean(r["cover"][tag]["null"]) if r["cover"][tag]["null"] else 0 for r in rs])
            diff = obs - nm
            sd = float(np.std(diff, ddof=1)) if len(diff) > 1 else 0.0
            z = float(np.mean(diff) / (sd / math.sqrt(len(diff)))) if sd > 0 else 0.0
            c[f"cover_{tag}_ratio"] = float(obs.sum() / max(1e-9, nm.sum()))
            c[f"cover_{tag}_z"] = z
        c["P2"] = within(c["cover_1_ratio"], 0.97, 1.03) and abs(c["cover_1_z"]) <= 4
        ap_o = sum(r["ap"]["obs"] for r in rs if "ap" in r)
        ap_n = sum(np.mean(r["ap"]["null"]) for r in rs if "ap" in r)
        c["ap_obs"], c["ap_null"] = ap_o, float(ap_n)
        c["P3_ratio"] = ap_o / ap_n if ap_n > 0 else None
        c["P3"] = "uninformative" if ap_n < 20 else within(c["P3_ratio"], 0.75, 1.33)
        mo = sum(v["obs_missing"] for r in rs for v in r["occ"].values())
        mn = sum(np.mean(v["null_missing"]) for r in rs for v in r["occ"].values())
        co = sum(v["obs_chi2"] for r in rs for v in r["occ"].values())
        cn = sum(np.mean(v["null_chi2"]) for r in rs for v in r["occ"].values())
        c["occ_missing_obs"], c["occ_missing_null"] = mo, float(mn)
        c["occ_chi2_ratio"] = co / cn if cn > 0 else None
        c["P4"] = (mo <= mn + 3 * math.sqrt(mn) + 3) and (cn == 0 or within(c["occ_chi2_ratio"], 0.9, 1.1))
        po = sum(sum(r["pair_obs"]) for r in rs)
        pm = sum(sum(r["pair_model"]) for r in rs)
        pu = sum(sum(r["pair_uniform"]) for r in rs)
        c["pair_ratio_model"], c["pair_ratio_uniform"] = po / pm, po / pu
        g_obs, g_null = [], []
        for r in rs:
            k, I = r["U_size"], r["b0"] - r["a0"] + 1
            g_obs += [g * k / I for g in r["gaps"]]
            g_null += [g * k / I for g in r["null_gaps"]]
        D, pv = ks_2samp(g_obs, g_null) if g_obs and g_null else (None, None)
        c["gap_KS_D"], c["gap_KS_p"] = D, pv
        c["P5"] = within(c["pair_ratio_model"], 0.85, 1.15) and (pv is None or pv > 1e-4)
        eo = sum(r["energy"]["obs"] for r in rs if "energy" in r)
        en = sum(np.mean(r["energy"]["null"]) for r in rs if "energy" in r)
        c["energy_ratio"] = eo / en
        c["P6"] = within(c["energy_ratio"], 0.9, 1.1)
        hu = sum(r["halves"]["U"] for r in rs)
        hp = sum(r["halves"]["U1"] * r["halves"]["U2"] / r["halves"]["I"] for r in rs)
        c["halves_ratio"] = hu / hp
        c["P5b"] = within(c["halves_ratio"], 0.85, 1.15)
        col_o = sum(r["crt"]["1"]["obs"]["collisions"] for r in rs)
        col_n = sum(np.mean([z["collisions"] for z in r["crt"]["1"]["null"]]) for r in rs if r["crt"]["1"]["null"])
        c["crt1_collision_ratio"] = col_o / col_n if col_n > 0 else None
        c["P4b"] = "uninformative" if col_n < 20 else within(c["crt1_collision_ratio"], 0.8, 1.25)
        lo = sum(r["lehman"]["in"] for r in rs)
        le = sum(r["U_size"] * r["lehman"]["measure"] / (r["b0"] - r["a0"] + 1) for r in rs)
        c["lehman_U_ratio"] = lo / le if le > 0 else None
        c["P8c"] = within(c["lehman_U_ratio"], 0.8, 1.25)
        tr_o = np.array([[t[1] for t in r["traj"]] for r in rs], float).sum(0)
        tr_m = np.array([r["traj_model"] for r in rs], float).sum(0)
        c["traj_ratio"] = (tr_o / tr_m).tolist()
        out[f"{key[0]}{key[1]}"] = c
    return out


def analyze_e2(recs):
    er, ef, mf_eq, er_plain = [], [], True, []
    for r in recs:
        for x in r["primes"]:
            er.append(x["ER"] / np.mean(x["ER_symnull"]))
            ef.append(x["EF"] / np.mean(x["EF_symnull"]))
            er_plain.append(x["ER"] / np.mean(x["ER_null"]))
            mf_eq &= (x["MF"] == x["MR"])
    res = {"ER_over_symnull_median": float(np.median(er)), "EF_over_symnull_median": float(np.median(ef)),
           "ER_over_null_median": float(np.median(er_plain)),
           "ER_range": [float(min(er)), float(max(er))], "EF_range": [float(min(ef)), float(max(ef))],
           "MF_equals_MR_all": bool(mf_eq), "count": len(er)}
    res["P7"] = within(res["ER_over_symnull_median"], 0.85, 1.2) and res["EF_over_symnull_median"] >= 1.3
    by = defaultdict(list)
    for r in recs:
        for x in r["primes"]:
            by[f"{r['n']}{r['shape']}"].append(x["EF"] / np.mean(x["EF_symnull"]))
    res["EF_ratio_by_config"] = {k: float(np.median(v)) for k, v in by.items()}
    return res


def analyze_e3(recs):
    unc = [r for r in recs if r["kind"] == "uncovered"]
    ctl = [r for r in recs if r["kind"] == "control"]
    res = {"n_uncovered": len(unc), "n_control": len(ctl),
           "control_uncovered_frac": sum(r["uncovered"] for r in ctl) / max(1, len(ctl))}
    for feat in ("lehman_in", "Lstar_le_X"):
        a, b = sum(r[feat] for r in unc), sum(r[feat] for r in ctl)
        rr, lo, hi = rr_ci(a, len(unc), b, len(ctl))
        res[feat] = {"unc_rate": a / len(unc), "ctl_rate": b / len(ctl), "RR": rr, "CI": [lo, hi]}
    D, pv = ks_uniform([r["nn_pct"] for r in unc])
    res["nn_pct_KS"] = [D, pv]
    res["near_ratio"] = sum(r["near"] for r in unc) / sum(r["near_exp"] for r in unc)
    held = {}
    for r in recs:
        if "held" in r:
            held[(r["n"], r["shape"], r["idx"])] = r["held"]
    hq = sum(h["Q"] for h in held.values())
    hp = sum(h["Q1"] * h["pred_ratio"] for h in held.values())
    res["held_ratio"] = hq / hp
    D2, p2 = ks_2samp([r["dq"] for r in unc], [r["dq"] for r in ctl])
    res["dq_KS"] = [D2, p2]
    D3, p3 = ks_2samp([math.log(r["Lstar"]) for r in unc], [math.log(r["Lstar"]) for r in ctl])
    res["Lstar_KS"] = [D3, p3]
    jump = lambda r: r["cfB2"] / r["cfB"]
    D4, p4 = ks_2samp([math.log(jump(r)) for r in unc], [math.log(jump(r)) for r in ctl])
    res["cfjump_KS"] = [D4, p4]
    li = res["lehman_in"]
    res["P8"] = (within(li["RR"], 0.7, 1.4) and li["CI"][0] <= 1 <= li["CI"][1]
                 and (pv is None or pv > 0.01) and within(res["near_ratio"], 0.8, 1.25)
                 and within(res["held_ratio"], 0.85, 1.15))
    by = defaultdict(lambda: [0, 0, 0, 0])
    for r in unc:
        by[f"{r['n']}{r['shape']}"][0] += r["lehman_in"]; by[f"{r['n']}{r['shape']}"][1] += 1
    for r in ctl:
        by[f"{r['n']}{r['shape']}"][2] += r["lehman_in"]; by[f"{r['n']}{r['shape']}"][3] += 1
    res["lehman_by_config"] = {k: [v[0] / max(1, v[1]), v[2] / max(1, v[3]), v[1], v[3]] for k, v in sorted(by.items())}
    return res


def analyze_e4(recs):
    worst = [r for r in recs if "hstar" in r]
    ctl = defaultdict(list)
    for r in recs:
        if r.get("kind") == "e4control":
            ctl[r["n"]].append(r)
    res = {"records": len(worst), "all_uncovered": all(r["uncovered"] for r in worst),
           "all_inU": all(r["inU"] for r in worst)}
    big = [r for r in worst if r["J"] >= 16]
    small = [r for r in worst if r["J"] <= 2 and r["signed"]]
    pct = []
    for r in big:
        if not ctl.get(r["n"]):
            continue
        cs = np.array([c["Lstar"] for c in ctl[r["n"]]])
        pct.append(float((np.count_nonzero(cs < r["Lstar"]) + 0.5 * np.count_nonzero(cs == r["Lstar"])) / len(cs)))
    res["big_J_count"] = len(big)
    res["Lstar_pct_KS"] = list(ks_uniform(pct))
    res["Lstar_pct_mean"] = float(np.mean(pct)) if pct else None
    nn = [r["nn_pct"] for r in big if r["nn_pct"] is not None]
    res["nn_pct_KS"] = list(ks_uniform(nn))
    res["U_over_model_bigJ"] = float(sum(r["U_size"] for r in big) / sum(r["E_model"] for r in big))
    res["lehman_in_bigJ"] = sum(r["lehman_in"] for r in big) / max(1, len(big))
    res["lehman_in_ctl"] = {n: sum(c["lehman_in"] for c in v) / len(v) for n, v in ctl.items() if v}
    res["smallJ_signed_censored"] = [(r["n"], r["Y"], r["J"], r["p_eq_q_mod_l1"], r["hstar"]) for r in small]
    res["P9"] = (res["all_uncovered"] and res["all_inU"]
                 and (res["Lstar_pct_KS"][1] is None or res["Lstar_pct_KS"][1] > 0.01)
                 and (res["nn_pct_KS"][1] is None or res["nn_pct_KS"][1] > 0.01))
    # the J <= 2 signed (1,1)-excluded artefact: censored records should have p = q (mod l1)
    res["P9b"] = all(r["p_eq_q_mod_l1"] for r in small if r["hstar"] == max(x["hstar"] for x in small if (x["n"], x["Y"]) == (r["n"], r["Y"])))
    return res


def main():
    s = {}
    e1, e2, e3, e4 = _load("e1"), _load("e2"), _load("e3"), _load("e4")
    if e1:
        s["e1"] = analyze_e1(e1)
    if e2:
        s["e2"] = analyze_e2(e2)
    if e3:
        s["e3"] = analyze_e3(e3)
    if e4:
        s["e4"] = analyze_e4(e4)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(s, open(OUT, "w"), indent=1, default=float)
    print(json.dumps(s, indent=1, default=float))


if __name__ == "__main__":
    main()
