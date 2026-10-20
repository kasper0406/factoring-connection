"""Analysis of experiments/raw/dlog_e1.jsonl against notes/dlog_prereg.md P1-P7 (HARNESS)."""
from __future__ import annotations

import json
import sys
from collections import Counter



def _binom_p(k, n, q):
    """Exact two-sided binomial p-value (sum of outcome probabilities <= that of k), q = 1/l rational."""
    from fractions import Fraction as F
    from math import comb
    q = F(q).limit_denominator(100)
    pk = [comb(n, i) * q ** i * (1 - q) ** (n - i) for i in range(n + 1)]
    return float(sum(x for x in pk if x <= pk[k]))


def main(path: str, out: str | None = None):
    recs = [json.loads(l) for l in open(path)]
    res = {}
    for blk in ("E1", "E2"):
        R = [r for r in recs if r["block"] == blk]
        res[f"P1_{blk}"] = f"{sum(r['correct'] for r in R)}/{len(R)}"
        res[f"P1_{blk}_Xequal"] = f"{sum(r['X_equal'] for r in R)}/{len(R)}"
    viol = Counter()
    for r in recs:
        for k, v in r["ctr"].items():
            if k.isupper() or k.startswith(("ROW_IDENTITY", "TORSION", "KERNEL", "NONUNIT", "WRONG")):
                viol[k] += v if isinstance(v, int) else 1
    res["violations"] = dict(viol)
    res["torsion_relations_total"] = sum(1 for r in recs for d in r["deps"] if "W_ok" in d)
    res["torsion_relations_gamma_nontrivial"] = sum(1 for r in recs for d in r["deps"] if d.get("gamma_nontrivial"))
    res["unit_rank_le_bound"] = sum(1 for r in recs if r["ctr"].get("unit_rank_implied", 0) <= sum(r["ctr"].get("r1r2", [9])) - 1)
    # P5
    pooled = {}
    for r in recs:
        if r["block"] != "E1":
            continue
        for d in r["deps"]:
            for l, vals in d.get("law", {}).items():
                c = pooled.setdefault(l, [0, 0])
                c[0] += sum(1 for v in vals if v == 0)
                c[1] += len(vals)
    res["P5"] = {l: {"zero": k, "n": n, "rate": round(k / n, 4), "pred": round(1 / int(l), 4),
                     "p": _binom_p(k, n, 1 / int(l))} for l, (k, n) in pooled.items() if n}
    # P6
    E2 = [r for r in recs if r["block"] == "E2"]
    a = sum(r["ctr"].get("basis_gamma_nontrivial_auxpm1", [0, 0])[0] for r in E2)
    b = sum(r["ctr"].get("basis_gamma_nontrivial_auxpm1", [0, 0])[1] for r in E2)
    res["P6_nontorsion_auxpm1"] = f"{a}/{b}"
    res["P6_unit_rank_full"] = f"{sum(1 for r in E2 if r['ctr'].get('unit_rank_implied') == sum(r['ctr'].get('r1r2', [9])) - 1)}/{len(E2)}"
    res["P7_ker2_minus_kerZ"] = Counter(r["ctr"].get("ker_mod2_minus_kerZ") for r in recs)
    res["r1r2"] = Counter(str(r["ctr"].get("r1r2")) for r in recs)
    res["per_run"] = [{k: r["ctr"].get(k) for k in ("rows_last", "pair_rows", "cols", "kernel_dim", "basis_nontrivial_gamma",
                                                       "unit_relations", "unit_rank_implied", "r1r2", "ybits", "fields_sampled",
                                                       "T_trials", "rounds", "atom_bad_prime_discarded", "ker_mod2_minus_kerZ")}
                      | {"block": r["block"], "bits": r["inst"]["bits"], "seed": r["inst"]["seed"], "wall": r["wall"]} for r in recs]
    s = json.dumps(res, indent=1, default=str)
    print(s)
    if out:
        open(out, "w").write(s)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
