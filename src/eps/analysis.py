"""Analysis of the eps/QR_2 battery (reads harness output; floats are used only here, for statistics).

Primary test (frozen in notes/eps_prereg.md section 5): the stratified permutation (Mantel/CMH-type
correlation) statistic.  Strata h = N mod 2^K_FLOOR.  For statistic x and target y,
    U = sum_h sum_{i in h} (x_i - xbar_h)(y_i - ybar_h),
    V = sum_h  [ sum_{i in h} (x_i - xbar_h)^2 ] [ sum_{i in h} (y_i - ybar_h)^2 ] / (n_h - 1),
    z = U / sqrt(V).
V is the exact variance of U under within-stratum permutation of y.  Any function of N mod 2^K_FLOOR has
x_i = xbar_h for every i, so U = V = 0: such statistics are reported DEGENERATE (z := 0) and can never fire.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
from typing import Dict, List, Optional, Tuple

import random

import numpy as np

from src.eps.operational import MENU, MENU_NAMES, S_BONFERRONI

K_FLOOR = 8
Z_THRESH = 6.0
Z_POWER = 0.8416   # one-sided normal quantile for 80% power
Z_REPLICATE = 4.0


def load(path: str) -> Dict:
    with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline())
        assert header["menu"] == MENU_NAMES, "menu mismatch"
        recs = [json.loads(line) for line in fh]
    X = np.array([r["x"] for r in recs], dtype=np.float64)
    N = [int(r["N"]) for r in recs]
    d = {"header": header, "X": X,
         "N": N,
         "Nmod": np.array([n % (1 << 16) for n in N], dtype=np.int64),
         "eps": np.array([r["eps"] for r in recs], dtype=np.int64),
         "qr2": np.array([-1 if r["qr2"] is None else r["qr2"] for r in recs], dtype=np.int64),
         "p": [int(r["p"]) for r in recs], "q": [int(r["q"]) for r in recs],
         "seed": np.array([r["seed"] for r in recs], dtype=np.int64),
         "t": np.array([r["t"] for r in recs]), "tg": [r["tg"] for r in recs],
         "first": [r["first"] for r in recs], "nonunits": [r["nonunits"] for r in recs],
         "gcd_checks": np.array([r["gcd_checks"] for r in recs], dtype=np.int64),
         "cf_nonunit": np.array([r.get("cf_nonunit", 0) for r in recs], dtype=np.int64),
         "cf_fb_ok": np.array([r.get("cf_fb_ok", 1) for r in recs], dtype=np.int64)}
    return d


def strat_z(x: np.ndarray, y: np.ndarray, strata: np.ndarray) -> Dict:
    """Stratified permutation z of x against y (both 1-D float arrays), strata integer labels."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    U = 0.0
    V = 0.0
    sxx_tot = 0.0
    max_dev2 = 0.0
    neff = 0
    for h in np.unique(strata):
        m = strata == h
        nh = int(m.sum())
        if nh < 2:
            continue
        xd = x[m] - x[m].mean()
        yd = y[m] - y[m].mean()
        sxx = float(xd @ xd)
        syy = float(yd @ yd)
        U += float(xd @ yd)
        V += sxx * syy / (nh - 1)
        sxx_tot += sxx
        if sxx > 0:
            max_dev2 = max(max_dev2, float((xd ** 2).max()))
        neff += nh - 1
    if V <= 1e-12 * max(1.0, sxx_tot):
        return {"z": 0.0, "degenerate": True, "U": U, "V": V, "neff": neff, "lindeberg": None}
    return {"z": U / math.sqrt(V), "degenerate": False, "U": U, "V": V, "neff": neff,
            "lindeberg": max_dev2 / sxx_tot}


def raw_z(x: np.ndarray, y: np.ndarray) -> float:
    """Unstratified point-biserial z = r sqrt(n-1) (for comparison only)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.std() == 0 or y.std() == 0:
        return 0.0
    r = float(np.corrcoef(x, y)[0, 1])
    return r * math.sqrt(len(x) - 1)


def mde_rho(neff: int, z_thr: float = Z_THRESH) -> float:
    return (z_thr + Z_POWER) / math.sqrt(neff)


def floor_baseline(Nmod: np.ndarray, y: np.ndarray, ks=range(2, 13)) -> List[Dict]:
    """Best predictor of y from N mod 2^k: majority label per class, 2-fold cross-validated
    (train on even positions, test on odd and vice versa).  Ties and unseen classes -> overall majority."""
    out = []
    y = np.asarray(y)
    idx = np.arange(len(y))
    for k in ks:
        cls = Nmod % (1 << k)
        correct = 0
        for tr, te in ((idx % 2 == 0, idx % 2 == 1), (idx % 2 == 1, idx % 2 == 0)):
            glob = 1 if (y[tr] == 1).mean() >= 0.5 else 0
            pred = {}
            for c in np.unique(cls[tr]):
                yy = y[tr][cls[tr] == c]
                fr = (yy == 1).mean()
                pred[c] = 1 if fr > 0.5 else (0 if fr < 0.5 else glob)
            pt = np.array([pred.get(c, glob) for c in cls[te]])
            correct += int(((pt == 1) == (y[te] == 1)).sum())
        acc = correct / len(y)
        out.append({"k": k, "classes": int(len(np.unique(cls))), "cv_acc": acc,
                    "se": math.sqrt(0.25 / len(y))})
    return out


def chi2_floor(Nmod: np.ndarray, y: np.ndarray, k: int = K_FLOOR) -> Dict:
    cls = Nmod % (1 << k)
    tot1 = (y == 1).mean()
    chi2 = 0.0
    dof = -1
    for c in np.unique(cls):
        m = cls == c
        n1 = float((y[m] == 1).sum())
        e1 = m.sum() * tot1
        e0 = m.sum() * (1 - tot1)
        chi2 += (n1 - e1) ** 2 / e1 + ((m.sum() - n1) - e0) ** 2 / e0
        dof += 1
    # Wilson-Hilferty normal approximation of the chi-square tail
    zwh = ((chi2 / dof) ** (1 / 3) - (1 - 2 / (9 * dof))) / math.sqrt(2 / (9 * dof)) if dof > 0 else 0.0
    return {"k": k, "chi2": chi2, "dof": dof, "z_wh": zwh}


def pc_check(x: np.ndarray, y01: np.ndarray, strata: np.ndarray) -> Dict:
    """Positive control: indicator x in {0,1}, claimed one-sided x=1 => y=1."""
    n = len(x)
    fire = x == 1
    f = float(fire.mean())
    viol = int((fire & (y01 != 1)).sum())
    pi = float((y01 == 1).mean())
    zs = strat_z(x, y01, strata)
    rho_pred = math.sqrt(f * (1 - pi) / ((1 - f) * pi)) if 0 < f < 1 else 0.0
    z_pred = rho_pred * math.sqrt(n - 1)
    return {"fires": int(fire.sum()), "freq": f, "violations": viol, "pi": pi, "z_obs": zs["z"],
            "z_pred": z_pred, "detected": zs["z"] >= Z_THRESH,
            "consistent": abs(zs["z"] - z_pred) <= 3.0, "powered": z_pred >= 8.0}


def analyze(d: Dict, target: str, shuffle_seed: int = 0) -> Dict:
    """target 'eps' (on all of D_n) or 'qr2' (on N = 1 mod 8)."""
    if target == "eps":
        mask = np.ones(len(d["eps"]), dtype=bool)
        y01 = (d["eps"] == 1).astype(np.int64)
        pcs = [m for m, _, r in MENU if r == "pc_eps"]
    else:
        mask = d["qr2"] >= 0
        y01 = d["qr2"].copy()
        pcs = [m for m, _, r in MENU if r == "pc_qr2"]
    X = d["X"][mask]
    y = y01[mask]
    Nmod = d["Nmod"][mask]
    strata = Nmod % (1 << K_FLOOR)
    res = {"target": target, "n_inst": int(mask.sum()), "frac_pos": float(y.mean()),
           "strata": int(len(np.unique(strata)))}
    rows = []
    for j, (name, grp, role) in enumerate(MENU):
        x = X[:, j]
        zs = strat_z(x, y, strata)
        rows.append({"stat": name, "group": grp, "role": role, "z": zs["z"], "degenerate": zs["degenerate"],
                     "raw_z": raw_z(x, y), "lindeberg": zs["lindeberg"], "neff": zs["neff"],
                     "n_distinct": int(len(np.unique(x))), "mean": float(x.mean())})
    res["rows"] = rows
    neff = max(r["neff"] for r in rows)
    res["mde_rho"] = mde_rho(neff)
    res["pc"] = {m: pc_check(X[:, MENU_NAMES.index(m)], y, strata) for m in pcs}
    res["floor"] = floor_baseline(Nmod, y)
    res["floor_chi2"] = chi2_floor(Nmod, y)
    # planted pipeline control: x = y w.p. 0.525 (rho = 0.05), from random.Random(f"planted|{seed}")
    seeds = d["seed"][mask]
    flip = np.array([random.Random(f"planted|{int(s)}").random() >= 0.525 for s in seeds])
    x_pl = np.where(flip, 1 - y, y)
    zpl = strat_z(x_pl, y, strata)
    res["planted"] = {"z": zpl["z"], "z_exp": 0.05 * math.sqrt(zpl["neff"]),
                      "ok": abs(zpl["z"] - 0.05 * math.sqrt(zpl["neff"])) <= 3.0}
    # calibration: permute y within strata (fixed seed) and recompute every z
    rng = np.random.default_rng(shuffle_seed)
    ys = y.copy()
    for h in np.unique(strata):
        idx = np.where(strata == h)[0]
        ys[idx] = ys[rng.permutation(idx)]
    zsh = [strat_z(X[:, j], ys, strata) for j in range(len(MENU))]
    zz = np.array([s["z"] for s in zsh if not s["degenerate"]])
    res["shuffle"] = {"n_nondeg": int(len(zz)), "mean_z2": float((zz ** 2).mean()),
                      "max_abs_z": float(np.abs(zz).max())}
    # hits: menu or nc statistics (not positive controls) with |z| >= threshold
    res["hits"] = [r["stat"] for r in rows if r["role"] in ("menu", "nc") and abs(r["z"]) >= Z_THRESH]
    nondeg = np.array([r["z"] for r in rows if r["role"] in ("menu", "nc") and not r["degenerate"]])
    res["null_summary"] = {"n_nondeg": int(len(nondeg)), "mean_z2": float((nondeg ** 2).mean()),
                           "max_abs_z": float(np.abs(nondeg).max())}
    return res


def audit(d: Dict, stat: str, target: str) -> Dict:
    """Leakage-audit steps 2-4 of notes/eps_prereg.md section 6 for one statistic."""
    from src.eps.operational import compute_stats   # operational recomputation from N only
    j = MENU_NAMES.index(stat)
    mask = np.ones(len(d["eps"]), dtype=bool) if target == "eps" else d["qr2"] >= 0
    y = ((d["eps"] == 1).astype(np.int64) if target == "eps" else d["qr2"])[mask]
    x = d["X"][mask, j]
    Ns = [N for N, m in zip(d["N"], mask) if m]
    # step 2: recompute (first 2000 instances, to bound cost) from N alone
    sub = range(min(2000, len(Ns)))
    mism = sum(1 for i in sub if compute_stats(Ns[i])[0][stat] != x[i])
    Nm = np.array([N % (1 << 12) for N in Ns])
    N105 = np.array([N % 105 for N in Ns])
    ps = [p for p, m in zip(d["p"], mask) if m]
    qs = [q for q, m in zip(d["q"], mask) if m]
    ratio = np.array([16 * q // p for p, q in zip(ps, qs)], dtype=float)
    lead = d["X"][mask, MENU_NAMES.index("nc_lead16")]
    sp = (d["seed"][mask] % 2).astype(float)
    s8 = Nm % 256
    # CF-group diagnostic (prereg section 6.6): restrict to instances whose period does not close within 1024
    first_q1 = np.array([f.get("Q1", 0) for f, m in zip(d["first"], mask) if m])
    nc = first_q1 == 0
    z_noclose = strat_z(x[nc], y[nc], s8[nc])["z"] if nc.sum() > 100 else None
    return {"stat": stat, "target": target, "recompute_checked": len(sub), "recompute_mismatch": mism,
            "z_period_not_closed": z_noclose, "n_period_not_closed": int(nc.sum()),
            "z_strata8": strat_z(x, y, s8)["z"],
            "z_strata12": strat_z(x, y, Nm)["z"],
            "z_strata8x105": strat_z(x, y, s8 * 105 + N105)["z"],
            "z_vs_ratio": strat_z(x, ratio, s8)["z"], "z_vs_lead16": strat_z(x, lead, s8)["z"],
            "z_vs_seedparity": strat_z(x, sp, s8)["z"]}


def load_curve(path: str) -> Dict:
    with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline())
        recs = [json.loads(line) for line in fh]
    return {"header": header, "recs": recs}


def curve_summary(c: Dict, Ts=(64, 256, 1024, 4096, 16384)) -> Dict:
    recs = c["recs"]
    out = {"n": c["header"]["n"], "count": len(recs), "frac_N1mod8": float(np.mean([r["qr2"] is not None for r in recs])),
           "frac_eps_pos": float(np.mean([r["eps"] == 1 for r in recs]))}
    for k in ("pc1", "pc2", "pc3", "pcq2", "pcq3"):
        sub = recs if not k.startswith("pcq") else [r for r in recs if r["qr2"] is not None]
        firsts = np.array([r["first"][k] for r in sub])
        if k.startswith("pcq"):
            viol = sum(1 for r in sub if r["first"][k] > 0 and r["qr2"] != 1)
        else:
            viol = sum(1 for r in sub if r["first"][k] > 0 and r["eps"] != 1)
        out[k] = {"violations": viol, "n": len(sub),
                  "freq": {T: float(((firsts > 0) & (firsts <= T)).mean()) for T in Ts},
                  "fires": {T: int(((firsts > 0) & (firsts <= T)).sum()) for T in Ts}}
    out["cf_nonunit"] = int(sum(r["cf_nonunit"] for r in recs))
    out["cf_fb_not_ok"] = int(sum(1 for r in recs if not r["cf_fb_ok"]))
    return out


def nonunit_summary(d: Dict) -> Dict:
    ev = sum(len(x) for x in d["nonunits"])
    checks = int(d["gcd_checks"].sum())
    exp = float(sum(c * (1.0 / p + 1.0 / q) for c, p, q in zip(d["gcd_checks"], d["p"], d["q"])))
    proper = sum(1 for x, N in zip(d["nonunits"], d["N"]) for w, g in x if 1 < g < N)
    return {"events": ev, "proper": proper, "gcd_checks": checks, "null_expected": exp,
            "cf_nonunit": int(d["cf_nonunit"].sum()), "cf_fb_not_ok": int((d["cf_fb_ok"] == 0).sum())}


def cost_summary(d: Dict) -> Dict:
    t = d["t"]
    groups = {}
    for g in d["tg"][0]:
        groups[g] = float(np.mean([x[g] for x in d["tg"]]))
    return {"mean_s": float(t.mean()), "max_s": float(t.max()), "total_cpu_h": float(t.sum() / 3600),
            "group_mean_s": groups}


def pc_frequency(d: Dict, Ts=(16, 32, 64, 128, 256, 512, 1024)) -> Dict:
    out = {}
    for k in ("pc1", "pc2", "pc3", "pcq2", "pcq3", "Q1"):
        firsts = np.array([f[k] for f in d["first"]])
        sub = firsts
        if k.startswith("pcq"):
            sub = firsts[d["qr2"] >= 0]
        out[k] = {T: float(((sub > 0) & (sub <= T)).mean()) for T in Ts}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--out", default="experiments/analysis/eps_battery_summary.json")
    ap.add_argument("--curve", action="store_true", help="paths are supplementary PC-curve files")
    ap.add_argument("--audit", nargs="*", default=[], help="stat:target pairs to audit (section 6)")
    a = ap.parse_args(argv)
    if a.curve:
        cs = [curve_summary(load_curve(p)) for p in a.paths]
        for c in cs:
            print(c["n"], {k: (c[k]["violations"], c[k]["freq"][1024], c[k]["freq"][16384])
                           for k in ("pc1", "pc2", "pc3", "pcq2", "pcq3")})
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        with open(a.out, "w") as fh:
            json.dump(cs, fh, indent=1)
        return
    summary = {}
    for path in a.paths:
        d = load(path)
        n = d["header"]["n"]
        key = os.path.basename(path).replace(".jsonl.gz", "")
        summary[key] = {"n": n, "header": {k: v for k, v in d["header"].items() if k != "menu"},
                        "eps": analyze(d, "eps"), "qr2": analyze(d, "qr2"),
                        "nonunits": nonunit_summary(d), "cost": cost_summary(d),
                        "pc_frequency": pc_frequency(d),
                        "audits": [audit(d, *sa.split(":")) for sa in a.audit]}
        e, q = summary[key]["eps"], summary[key]["qr2"]
        print(f"{key}: n={n} eps hits={e['hits']} qr2 hits={q['hits']} "
              f"null mean z^2 eps={e['null_summary']['mean_z2']:.2f} qr2={q['null_summary']['mean_z2']:.2f}")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(summary, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
