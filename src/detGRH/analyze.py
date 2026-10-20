"""detGRH E1 analysis: evaluates the pre-registered predictions (notes/detGRH_prereg.md).

  PYTHONPATH=. .venv/bin/python -m src.detGRH.analyze experiments/raw/detGRH_e1.jsonl \
      experiments/analysis/detGRH_e1_summary.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict

from src.common.metrics import read_jsonl


def main(raw: str, out: str) -> None:
    recs = list(read_jsonl(raw))
    s: dict = {}
    a = [r for r in recs if r["exp"] == "E1a"][0]
    s["E1a"] = {"MR1_tested": a["tested_MR1"], "MR1_violations": len(a["viol_MR1"]),
                "MR2_tested": a["tested_MR2"], "MR2_violations": len(a["viol_MR2"]),
                "P-MR1": len(a["viol_MR1"]) == 0, "P-MR2": len(a["viol_MR2"]) == 0}
    inst = [r for r in recs if r["exp"] == "E1bcd" and "factor_exit" not in r]
    s["instances"] = len(inst)
    s["factor_exits"] = sum(1 for r in recs if r["exp"] == "E1bcd" and "factor_exit" in r)
    gen_ok = [r for r in inst if r["gen_small_primes"]]
    s["gen_small_primes"] = f"{len(gen_ok)}/{len(inst)}"
    full = [r for r in gen_ok if r["full_rank"]]
    s["full_rank"] = f"{sum(r['full_rank'] for r in inst)}/{len(inst)}"
    z1 = [r for r in full if r["miller_factor_proper"]]
    s["P-Z1"] = {"factored": f"{len(z1)}/{len(full)}", "pass": len(z1) == len(full)}
    c_f2 = [r for r in inst if r["ctrl_f2_fruitful_basis"] == 0]
    c_z = [r for r in full if r["ctrl_miller_proper"]]
    s["P-Z2"] = {"ctrl_f2_fruitless": f"{len(c_f2)}/{len(inst)}", "ctrl_z_factored": f"{len(c_z)}/{len(full)}",
                 "pass": len(c_f2) == len(inst) and len(c_z) == len(full)}
    s["f2_route_fruitful"] = f"{sum(r['f2_fruitful'] for r in inst)}/{len(inst)}"
    bases = defaultdict(int)
    for r in full:
        bases[r["miller_base"]] += 1
    s["miller_base_hist"] = dict(sorted(bases.items()))
    # E1c
    g = [r for r in inst if r["GB_full"]]
    fl1 = [r for r in g if r["dimW"] != r["omega"] - 1]
    s["P-FL1"] = {"tested": len(g), "violations": len(fl1), "pass": not fl1}
    withc = [r for r in g if "c" in r]
    fl2_t = [r for r in withc if r["c"] < r["dimW"]]
    fl2_v = [r for r in fl2_t if not r["f2_fruitful"]]
    ge = [r for r in withc if r["c"] >= r["dimW"]]
    s["P-FL2"] = {"premise_true": len(fl2_t), "violations": len(fl2_v), "pass": not fl2_v,
                  "fruitful_rate_when_c_ge_dimW": f"{sum(r['f2_fruitful'] for r in ge)}/{len(ge)}"}
    chist = defaultdict(int)
    for r in withc:
        chist[r["c"]] += 1
    s["c_hist"] = dict(sorted(chist.items()))
    # E1d
    viol = sum(r["cy_violations"] for r in inst)
    s["P-CY1"] = {"violations": viol, "nontrivial_defects": sum(r["defect_nontriv"] for r in inst),
                  "gcd_proper": sum(r["defect_gcd_proper"] for r in inst), "pass": viol == 0}
    cy2, cy3 = {}, {}
    for grp, sel in (("semiprime", lambda r: r["omega"] == 2), ("tri", lambda r: r["omega"] == 3)):
        rr = [r for r in inst if sel(r)]
        sp = sum(r["smooth_pairs"] for r in rr)
        dn = sum(r["defect_nontriv"] for r in rr)
        cn = sum(r["chi_n"] for r in rr)
        cd = sum(r["chi_nondiag"] for r in rr)
        lo, hi = (0.30, 0.70) if grp == "semiprime" else (0.55, 0.95)
        cy2[grp] = {"pairs": sp, "nontrivial": dn, "rate": round(dn / sp, 4) if sp else None,
                    "window": [lo, hi], "pass": bool(sp) and lo <= dn / sp <= hi}
        cy3[grp] = {"seeds": cn, "nondiag": cd, "rate": round(cd / cn, 4) if cn else None,
                    "window": [lo, hi], "pass": bool(cn) and lo <= cd / cn <= hi}
    s["P-CY2"], s["P-CY3"] = cy2, cy3
    # per class/bits table
    tab = {}
    for r in inst:
        k = f"{r['cls']}/{r['bits']}"
        t = tab.setdefault(k, {"n": 0, "full_rank": 0, "miller": 0, "f2": 0, "pairs": 0, "defects": 0,
                               "chi_n": 0, "chi_nondiag": 0, "kernel_rel": 0})
        t["n"] += 1
        t["full_rank"] += r["full_rank"]
        t["miller"] += r.get("miller_factor_proper", False)
        t["f2"] += r["f2_fruitful"]
        t["pairs"] += r["smooth_pairs"]
        t["defects"] += r["defect_nontriv"]
        t["chi_n"] += r["chi_n"]
        t["chi_nondiag"] += r["chi_nondiag"]
        t["kernel_rel"] += r["n_kernel"]
    s["table"] = tab
    json.dump(s, open(out, "w"), indent=1)
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
