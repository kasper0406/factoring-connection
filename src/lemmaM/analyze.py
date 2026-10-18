"""Analysis of E1/E2 raw JSONL (frozen with the pre-registration)."""
from __future__ import annotations
import json
import math
import sys
from collections import defaultdict


def e1(path):
    recs = [json.loads(l) for l in open(path)]
    out = {}
    acc = sum(1 for r in recs if r.get("accidental_factor"))
    out["accidental_factor_instances"] = acc
    for mode in ["pure", "mixed"]:
        for grp in ["D", "Q"]:
            for cls in ["semi", "tri"]:
                rows = [r for r in recs if r.get("mode") == mode and "res" in r
                        and (r["cls"] == "tri") == (cls == "tri")]
                p0 = 0.75 if cls == "tri" else 0.5
                n = sum(r["res"][grp][0] for r in rows)
                k = sum(r["res"][grp][1] for r in rows)
                if n == 0:
                    continue
                # cluster-robust z: instance-level residuals
                resid = [(r["res"][grp][1] - p0 * r["res"][grp][0]) for r in rows]
                var = sum(x * x for x in resid)
                z = sum(resid) / math.sqrt(var) if var > 0 else 0.0
                zb = (k - p0 * n) / math.sqrt(n * p0 * (1 - p0))
                out[f"{mode}|{grp}|{cls}"] = dict(n=n, fruitful=k, rate=round(k / n, 4),
                                                   pred=p0, z_cluster=round(z, 2), z_binom=round(zb, 2))
    out["nonsquare"] = sum(r["res"]["nonsquare"] for r in recs if "res" in r)
    out["flip_ok"] = sum(r["res"]["flip_ok"] for r in recs if "res" in r)
    out["flip_fail"] = sum(r["res"]["flip_fail"] for r in recs if "res" in r)
    # P5: r-trials per accepted relation, plain vs pure, by bits
    by = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in recs:
        if r.get("mode") in ("plain", "pure") and "counters" in r:
            c = r["counters"]
            by[r["bits"]][r["mode"]][0] += c["r_trials"]
            by[r["bits"]][r["mode"]][1] += c["desc_accepted"]
    p5 = {}
    for b, d in sorted(by.items()):
        tp = d["plain"][0] / max(1, d["plain"][1])
        td = d["pure"][0] / max(1, d["pure"][1])
        p5[b] = dict(trials_per_rel_plain=round(tp, 2), trials_per_rel_descent=round(td, 2),
                     ratio=round(tp / td, 2))
    out["P5"] = p5
    sq = defaultdict(lambda: [0, 0, 0])
    for r in recs:
        if r.get("mode") == "pure" and "counters" in r:
            c = r["counters"]
            sq[r["bits"]][0] += c["sq_attempts"]
            sq[r["bits"]][1] += c["sq_success"]
            sq[r["bits"]][2] += c["desc_accepted"]
    out["special_q"] = {b: dict(attempts=v[0], success=v[1], per_descent=round(v[0] / max(1, v[2]), 2))
                        for b, v in sorted(sq.items())}
    return out


def e2(path):
    recs = [json.loads(l) for l in open(path)]
    out = dict(capped=sum(1 for r in recs if r.get("capped")),
               identity_fail=sum(1 for r in recs if "P" in r and not (r["ok_Q"] and r["ok_A"] and r["ok_unit"])))
    tab = defaultdict(lambda: defaultdict(int))
    for r in recs:
        if "P" in r:
            tab[r["cls"]][r["unit"]] += 1
    out["unit_classes"] = {k: dict(v) for k, v in tab.items()}
    xs, ys = [], []
    for r in recs:
        if "P" in r:
            xs.append(r["bits"])
            ys.append(math.log2(r["P"]))
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    out["log2P_slope"] = round(slope, 3)   # float used only for the fit (allowed)
    return out


if __name__ == "__main__":
    which, path = sys.argv[1], sys.argv[2]
    print(json.dumps(e1(path) if which == "e1" else e2(path), indent=1))
