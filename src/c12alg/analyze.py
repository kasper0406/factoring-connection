"""c12alg analysis of E1a, E1b, E2 exactly as registered in notes/c12alg_prereg.md.  HARNESS.

Usage:
  .venv/bin/python -m src.c12alg.analyze E1a experiments/raw/c12alg_e1a.jsonl
  .venv/bin/python -m src.c12alg.analyze E1b experiments/raw/c12alg_e1b.jsonl
  .venv/bin/python -m src.c12alg.analyze E2  experiments/raw/c12alg_e2.jsonl
Each writes <input>_summary.json next to the input.
"""
from __future__ import annotations

import json
import random
import sys
from collections import defaultdict
from fractions import Fraction
from math import gcd, sqrt
from statistics import median


def e1a(path):
    R = [json.loads(l) for l in open(path)]
    S = {"n": len(R)}
    S["P1a_tv_le_bound"] = [sum(r["tv"] <= r["bound"] + 1e-15 for r in R), len(R)]
    S["P1a_exact"] = all(Fraction(r["tv_num"], r["tv_den"]) <= Fraction(r["bound"]).limit_denominator(10 ** 9) + Fraction(1, 10 ** 12)
                         for r in R)
    S["P1b_xmax_le_bound"] = [sum(r["xmax"] <= r["xmax_bound"] for r in R), len(R)]
    by = defaultdict(list)
    for r in R:
        by[(r["d"], r["V"], r["H"])].append(r)
    S["configs"] = {}
    ok_c = True
    for k, rs in sorted(by.items()):
        d, V, H = k
        nonvac = [r for r in rs if r["bound"] < 1]
        ratios = [r["tv"] * (2 * H + 1) / (d * V * V) for r in nonvac]
        med = median(ratios) if ratios else None
        S["configs"][str(k)] = {"n": len(rs), "nonvacuous": len(nonvac), "median_tv": median(r["tv"] for r in rs),
                                "max_tv": max(r["tv"] for r in rs), "bound": rs[0]["bound"],
                                "median_tv_over_dV2_per_H": med,
                                "max_tv_over_bound": max(r["tv"] / r["bound"] for r in rs)}
        if med is not None and not (0.05 <= med <= 4):
            ok_c = False
    S["P1c_band_all_configs"] = ok_c
    S["P1c_monotone"] = S["configs"]["(3, 2, 120)"]["median_tv"] < S["configs"]["(3, 2, 40)"]["median_tv"]
    return S


def _boot(fields, stat, B=2000, seed=1):
    rng = random.Random(seed)
    keys = sorted(fields)
    vals = []
    for _ in range(B):
        samp = [rng.choice(keys) for _ in keys]
        vals.append(stat([r for k in samp for r in fields[k]]))
    vals.sort()
    return [vals[int(0.025 * B)], vals[int(0.975 * B) - 1]]


def e1b(path):
    R = [json.loads(l) for l in open(path)]
    fields = defaultdict(list)
    for r in R:
        fields[r["field"]].append(r)
    S = {"n_pairs": len(R), "n_fields": len(fields)}

    def pooled(key):
        return lambda rs: sum(r["J"] for r in rs) / sum(r[key] for r in rs)

    def slope(rs):
        xs = [r["J_full"] / (r["M1"] * r["M2"]) - 1 for r in rs]
        ys = [r["J"] / (r["M1"] * r["M2"]) - 1 for r in rs]
        return sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)

    def resid(rs):
        return sum((r["J"] - r["J_full"]) / (r["M1"] * r["M2"]) for r in rs) / len(rs)

    def slope_gcd(rs):
        xs = [r["J_gcd"] / (r["M1"] * r["M2"]) - 1 for r in rs]
        ys = [r["J"] / (r["M1"] * r["M2"]) - 1 for r in rs]
        return sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)

    S["P1d_pooled_J_over_Jgcd"] = [pooled("J_gcd")(R), _boot(fields, pooled("J_gcd"))]
    S["P1e_pooled_J_over_Jfull"] = [pooled("J_full")(R), _boot(fields, pooled("J_full"))]
    S["P1f_slope_full"] = [slope(R), _boot(fields, slope)]
    S["P1g_mean_resid"] = [resid(R), _boot(fields, resid)]
    S["descr_slope_gcd_only"] = slope_gcd(R)
    S["descr_pooled_J_over_M1M2"] = sum(r["J"] for r in R) / sum(r["M1"] * r["M2"] for r in R)
    ex = [r["J"] / (r["M1"] * r["M2"]) - 1 for r in R]
    S["descr_excess_range"] = [min(ex), max(ex)]
    S["descr_arch_range"] = [min(r["arch"] for r in R), max(r["arch"] for r in R)]
    p = S
    S["verdicts"] = {
        "P1d": 0.95 <= p["P1d_pooled_J_over_Jgcd"][0] <= 1.08,
        "P1e": 0.97 <= p["P1e_pooled_J_over_Jfull"][0] <= 1.03 and 0.95 <= p["P1e_pooled_J_over_Jfull"][1][0]
        and p["P1e_pooled_J_over_Jfull"][1][1] <= 1.05,
        "P1f": 0.8 <= p["P1f_slope_full"][0] <= 1.2,
        "P1g": abs(p["P1g_mean_resid"][0]) <= 0.03,
    }
    return S


def _mobius(n):
    res, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    if m > 1:
        res = -res
    return res


def p_content(W: int, H: int, d: int):
    """Pr[gcd(w, c_0..c_{d-1}) > 1], w uniform in (W, 2W], c uniform in [-H, H]^d.  Exact: enumerate c, and for
    g = gcd(c) (g = 0 when c = 0, then gcd(w, 0) = w > 1) count w in (W, 2W] with gcd(w, g) > 1 by Moebius over the
    squarefree divisors of g.  (v1 truncated the Moebius sum at e <= 2H + 1; corrected after codex c12alg_2 #13.)"""
    import itertools
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def frac_coprime(g: int) -> Fraction:
        ps = [p for p in range(2, g + 1) if g % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
        tot = Fraction(0)
        for mask in range(1 << len(ps)):
            e = 1
            for i, p in enumerate(ps):
                if mask >> i & 1:
                    e *= p
            sign = -1 if bin(mask).count("1") % 2 else 1
            tot += sign * Fraction((2 * W) // e - W // e, W)
        return tot
    num = Fraction(0)
    n = 0
    for c in itertools.product(range(-H, H + 1), repeat=d):
        g = 0
        for ci in c:
            g = gcd(g, ci)
        n += 1
        num += 1 if g == 0 else (1 - frac_coprime(g))
    return float(num / n), float(Fraction(1, (2 * H + 1) ** d))


def e2(path):
    R = [json.loads(l) for l in open(path)]
    S = {"blocks": {}}
    viol = defaultdict(int)
    dep_sq = dep_nonsq = 0
    dd = defaultdict(lambda: [0, 0])
    cont_obs = cont_den = 0
    cont_pred_num = 0.0
    cont_var = 0.0
    for r in R:
        b = S["blocks"].setdefault(r["block"], {"n": 0, "ok": 0, "errors": []})
        b["n"] += 1
        b["ok"] += bool(r["ok"])
        if "error" in r:
            b["errors"].append([r["spec"], r["error"]])
        for L in r.get("log", []):
            if L["step"] != "split":
                continue
            c = L["ctr"]
            for k in ("ROW_IDENTITY_VIOLATION", "CONGRUENCE_VIOLATION", "CU_VIOLATION", "FORCED_RESIDUE_VIOLATION",
                      "IR_DW_CLAIM_VIOLATION"):
                viol[k] += c.get(k, 0)
            viol["cu_checked"] += c.get("cu_checked", 0)
            dep_sq += c.get("dep_square", 0)
            dep_nonsq += c.get("dep_nonsquare", 0)
            if r["block"] == "E2b":
                for dep in L["deps"]:
                    dd[L["omega"]][0] += 1
                    dd[L["omega"]][1] += bool(dep["fruitful"])
            # content: denominator = fields that passed gcd and degree checks (sampled - rej_gcdN - degree drops)
            den = c.get("fields_sampled", 0) - c.get("rej_gcdN", 0) - c.get("rej_degree_drop", 0)
            obs = c.get("rej_content", 0)
            W = c["params"]["W"]
            H = c["params"]["H"]
            d = c["params"]["d"]
            pc, _ = p_content(W, H, d)
            cont_obs += obs
            cont_den += den
            cont_pred_num += pc * den
            cont_var += pc * (1 - pc) * den
    S["P2a"] = {k: [v["ok"], v["n"]] for k, v in S["blocks"].items()}
    S["P2b_violations"] = dict(viol)
    S["P2c_squares"] = [dep_sq, dep_sq + dep_nonsq]
    z = (cont_obs - cont_pred_num) / sqrt(cont_var) if cont_var else None
    S["P2d_content"] = {"observed": cont_obs, "fields": cont_den, "predicted": cont_pred_num, "z": z}
    S["P2e_DD"] = {}
    for om, (n, f) in sorted(dd.items()):
        p0 = 1 - 2.0 ** (1 - om)
        zz = (f - n * p0) / sqrt(n * p0 * (1 - p0)) if n else None
        S["P2e_DD"][str(om)] = {"tested": n, "fruitful": f, "frac": f / n if n else None, "pred": p0, "z": zz}
    S["verdicts"] = {
        "P2a": all(v[0] == v[1] for v in S["P2a"].values()),
        "P2b": all(viol[k] == 0 for k in ("ROW_IDENTITY_VIOLATION", "CONGRUENCE_VIOLATION", "CU_VIOLATION",
                                           "FORCED_RESIDUE_VIOLATION", "IR_DW_CLAIM_VIOLATION")),
        "P2c": dep_nonsq == 0,
        "P2d": z is not None and abs(z) <= 3,
        "P2e": all(abs(v["z"]) <= 3 for v in S["P2e_DD"].values() if v["z"] is not None),
    }
    return S


if __name__ == "__main__":
    which, path = sys.argv[1], sys.argv[2]
    S = {"E1a": e1a, "E1b": e1b, "E2": e2}[which](path)
    out = path.replace(".jsonl", "_summary.json")
    json.dump(S, open(out, "w"), indent=1, default=str)
    print(json.dumps(S, indent=1, default=str))
