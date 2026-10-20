"""POST-HOC (not pre-registered; added after codex rank3_2 #10 noted that every pre-registered P2 check was vacuous,
error term >= box population). Re-checks the Lemma R2 / Theorem R1 count bound on the E1 (B) and (A) rank-3
instances with boxes large enough that error < box population, so a failing main term would be detectable.

Run: .venv/bin/python -m src.rank3.e1_posthoc   (writes experiments/raw/rank3_E1_posthoc.jsonl)
"""
from __future__ import annotations

import bisect
import json

from src.rank3.forms import box_abs_values, count_bound_terms, delta_F

OUT = "experiments/raw/rank3_E1_posthoc.jsonl"


def main() -> None:
    rows = [json.loads(l) for l in open("experiments/raw/rank3_E1.jsonl")]
    viol = checks = nonvac = 0
    with open(OUT, "w") as fh:
        for r_ in rows:
            if r_["r"] != 3 or (r_["d"], r_["d2"]) != (3, 3):
                continue
            f, g, P, Q = r_["f"], r_["g"], r_["P"], r_["Q"]
            info = delta_F(f, P, g, Q, 3)
            E = 40
            vals = box_abs_values(f, P, g, Q, 3, E)
            out = []
            for k in range(max(0, vals[0].bit_length() - 1), vals[-1].bit_length() + 1):
                X = 1 << k
                c = bisect.bisect_right(vals, X)
                main, err = count_bound_terms(3, 6, E, X, info)
                ok = bool(c <= main.upper() + err)
                informative = bool(c > err)  # count exceeds the error term: the main term is actually tested
                checks += 1
                viol += not ok
                nonvac += informative
                out.append({"log2X": k, "count": c, "main_mid": float(main.mid()), "err": err, "ok": ok,
                            "informative": informative, "ratio_count_minus_err_over_main":
                            (c - err) / float(main.mid()) if informative else None})
            rec = {"kind": r_["kind"], "seed": r_["seed"], "N": r_["N"], "E": E, "box": len(vals), "checks": out}
            fh.write(json.dumps(rec) + "\n")
            inf = [x["ratio_count_minus_err_over_main"] for x in out if x["informative"]]
            print(r_["kind"], r_["seed"], "box", len(vals), "err", out[0]["err"], "informative X", len(inf),
                  "max (count-err)/main", max(inf) if inf else None, flush=True)
    print(f"post-hoc P2: {checks} checks, {viol} violations, {nonvac} informative (count > error)")


if __name__ == "__main__":
    main()
