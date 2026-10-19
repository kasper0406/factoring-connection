"""Analysis of thmA E1/E1b (statistics exactly as registered in notes/thmA_prereg.md).  HARNESS.

Usage: .venv/bin/python -m src.thmA.analyze experiments/raw/thmA_e1.jsonl [summary.json]
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from math import sqrt


def main(path: str, out: str | None = None):
    recs = [json.loads(l) for l in open(path)]
    S: dict = {"blocks": {}}
    dd = defaultdict(lambda: [0, 0])            # omega -> [tested square deps, fruitful]
    squares = nonsq = 0
    row_viol = cong_viol = 0
    acc_obs, acc_null = 0, 0.0
    steps = defaultdict(lambda: defaultdict(float))
    split_kinds = defaultdict(int)
    for r in recs:
        b = S["blocks"].setdefault(r["block"], {"n": 0, "ok": 0, "errors": [], "by_family": {}})
        b["n"] += 1
        b["ok"] += bool(r["ok"])
        fam = r["spec"]["family"]
        bf = b["by_family"].setdefault(fam, [0, 0])
        bf[0] += 1
        bf[1] += bool(r["ok"])
        if "error" in r:
            b["errors"].append([r["spec"], r["error"]])
        for L in r.get("log", []):
            if L["step"] != "split":
                split_kinds[(r["block"], L["step"])] += 1
                continue
            c = L["ctr"]
            split_kinds[(r["block"], L["factor_kind"].split("_")[0])] += 1
            row_viol += c.get("ROW_IDENTITY_VIOLATION", 0)
            cong_viol += c.get("CONGRUENCE_VIOLATION", 0)
            squares += c.get("dep_square", 0)
            nonsq += c.get("dep_nonsquare", 0)
            for dep in L["deps"]:
                dd[L["omega"]][0] += 1
                dd[L["omega"]][1] += bool(dep["fruitful"])
            if r["block"] == "E1":
                lam = (3 * c.get("r_trials", 0) + 2 * c.get("sq_candidates", 0)
                       + 4 * c.get("fields_sampled", 0)) * L["sum_inv_p"]
                acc_null += lam
                acc_obs += L["factor_kind"].startswith("ACCIDENTAL")
            key = f'{r["block"]}|w{L["omega"]}|{c["n_bits"] // 8 * 8}b'
            st = steps[key]
            st["splits"] += 1
            for k in ("r_trials", "shape_ok", "descent_fail", "sq_searches", "sq_success", "rows",
                      "fields_sampled", "rows_last", "cols_last", "kernel_dim_last", "lv6_rounds", "seconds"):
                st[k] += c.get(k, 0)
            st["accidental_continued"] += len(c.get("accidental_continued", []))
    S["P1_E1_correct"] = [S["blocks"]["E1"]["ok"], S["blocks"]["E1"]["n"]]
    S["P2_row_identity_violations"] = row_viol
    S["P3_congruence_violations"] = cong_viol
    S["P4_square_fraction"] = [squares, squares + nonsq, squares / max(1, squares + nonsq)]
    p5 = {}
    for om, (n, k) in sorted(dd.items()):
        p = 1 - 2 ** (1 - om)
        z = (k - n * p) / sqrt(n * p * (1 - p)) if n else None
        p5[om] = {"tested": n, "fruitful": k, "rate": k / n if n else None, "pred": p, "z": z,
                  "testable": n >= 30, "holds": (n >= 30 and abs(z) <= 3)}
    S["P5_DD"] = p5
    S["P6_accidental"] = {"observed": acc_obs, "null_lambda": acc_null,
                          "threshold": acc_null + 3 * sqrt(acc_null) + 1,
                          "holds": acc_obs <= acc_null + 3 * sqrt(acc_null) + 1}
    pre = [r for r in recs if r["block"] == "E1" and r["spec"]["family"] in
           ("unbal", "small_x_semi", "ppow", "semi_sq", "prime")]
    S["P7_preprocessing"] = [sum(r["ok"] for r in pre), len(pre)]
    S["split_kinds"] = {f"{a}|{b}": v for (a, b), v in sorted(split_kinds.items())}
    desc = {}
    for key, st in sorted(steps.items()):
        n = st["splits"]
        desc[key] = {"splits": int(n),
                     "r_trials/split": round(st["r_trials"] / n),
                     "shape_ok_rate": round(st["shape_ok"] / max(1, st["r_trials"]), 4),
                     "descent_success_given_shape": round(st["rows"] / max(1, st["shape_ok"]), 4),
                     "sq_success_rate(gamma_hat)": round(st["sq_success"] / max(1, st["sq_searches"]), 4),
                     "fields/split": round(st["fields_sampled"] / n, 1),
                     "rows_last/split": round(st["rows_last"] / n), "cols/split": round(st["cols_last"] / n),
                     "kernel_dim/split": round(st["kernel_dim_last"] / n, 1),
                     "lv6_rounds/split": round(st["lv6_rounds"] / n, 2),
                     "accidental_continued/split": round(st["accidental_continued"] / n, 2),
                     "seconds/split": round(st["seconds"] / n, 2)}
    S["descriptive"] = desc
    txt = json.dumps(S, indent=1, default=str)
    print(txt)
    if out:
        open(out, "w").write(txt)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
