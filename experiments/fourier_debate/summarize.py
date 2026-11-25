"""Pooled first-hit costs and 95% paired whole-form bootstrap intervals."""
import json
from collections import Counter
import numpy as np

from common import HERE, REPS, SIZES

ARMS = ["pairs", "incidence", "blind", "matched"]
FIELDS = ["cpu_ms", "wall_ms", "branches", "tested", "failed"]
B = 10000


def interval(point, draws):
    assert np.isfinite(draws).all()
    low, high = np.quantile(draws, [0.025, 0.975])
    return [float(point), float(low), float(high)]


def main():
    data = [json.loads(s) for s in (HERE / "heldout.jsonl").read_text().splitlines()]
    assert len(data) == len(SIZES)*REPS
    rng = np.random.default_rng(252611)
    result = {}
    bands = {}
    for bits in SIZES:
        allrows = [r for r in data if r["bits"] == bits]
        rows = [r for r in allrows if "generation_error" not in r]
        n = len(rows)
        hit = np.array([bool(r["matched_points"]) for r in rows])
        W = rng.multinomial(n, np.full(n, 1/n), size=B)
        successes = hit.sum()
        denom = W @ hit
        assert denom.min() > 0
        full = np.array([r["unrestricted"] for r in rows])
        witnesses = np.array([len(r["matched_points"]) for r in rows])
        report = dict(forms=n, generation_errors=len(allrows)-n, B=rows[0]["B"],
                      successes=int(successes), empty=n-int(successes),
                      exactly_one=int((witnesses == 1).sum()), eligible_total=int(witnesses.sum()),
                      unrestricted_total=int(full.sum()), matched_mean=interval(witnesses.mean(), W@witnesses/n),
                      arms={}, ratios={})
        result[bits] = report
        points = {}
        draws = {}
        for arm in ARMS:
            entries = [r["arms"][arm] for r in rows]
            ops = [e["counts"] for e in entries]
            spatial = arm != "matched"
            node_key = "nodes" if spatial else "branches"
            fail_key = "dead_leaves" if spatial else "failed_slopes"
            vals = np.array([[e["cpu_ns"]/1e6, e["wall_ns"]/1e6, c[node_key],
                              c.get("norm_tests", 0), c.get(fail_key, 0)] for e, c in zip(entries, ops)])
            per = vals.sum(axis=0) / successes
            boot = (W @ vals) / denom[:, None]
            points[arm], draws[arm] = per, boot
            totals = sum((Counter(c) for c in ops), Counter())
            report["arms"][arm] = dict(
                per_hit={k: interval(v, bs) for k, v, bs in zip(FIELDS, per, boot.T)},
                total_cpu_s=float(vals[:,0].sum()/1000),
                total_wall_s=float(vals[:,1].sum()/1000), counters=totals)
        for name, numer, ref, field in (
            ("nodes_guided_over_blind", "incidence", "blind", "branches"),
            ("tests_guided_over_blind", "incidence", "blind", "tested"),
            ("cpu_incidence_over_matched", "incidence", "matched", "cpu_ms"),
            ("cpu_pairs_over_matched", "pairs", "matched", "cpu_ms"),
            ("cpu_incidence_over_blind", "incidence", "blind", "cpu_ms"),
            ("wall_incidence_over_matched", "incidence", "matched", "wall_ms"),
        ):
            i = FIELDS.index(field)
            p = points[numer][i] / points[ref][i]
            ds = draws[numer][:,i] / draws[ref][:,i]
            report["ratios"][name] = interval(p, ds)
            bands[bits, name] = ds
    result["largest_over_smallest"] = {}
    for name in result[SIZES[0]]["ratios"]:
        lo, hi = SIZES[0], SIZES[-1]
        a, b = (result[size]["ratios"][name][0] for size in (hi, lo))
        result["largest_over_smallest"][name] = interval(a/b, bands[hi,name]/bands[lo,name])
    (HERE / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    for b in SIZES:
        r = result[b]
        print(b, r["B"], r["successes"], r["exactly_one"], r["eligible_total"], r["unrestricted_total"])
        for arm, rec in r["arms"].items():
            print(arm, {k: [round(v, 3) for v in rec["per_hit"][k]] for k in ["branches", "cpu_ms","tested"]}, "total CPU s", round(rec["total_cpu_s"],3))
        print({k:[round(v,4) for v in vs] for k,vs in r["ratios"].items()})
    print("Trend", result["largest_over_smallest"])


if __name__ == "__main__":
    main()
