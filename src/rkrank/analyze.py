"""Analysis of rkrank E1 (predictions P0-P5 of notes/rkrank_prereg.md).
Usage: .venv/bin/python -m src.rkrank.analyze experiments/raw/rkrank_e1.jsonl
"""
from __future__ import annotations

import json
import sys


def main():
    recs = [json.loads(l) for l in open(sys.argv[1])]
    print(f"instances: {len(recs)}")
    P0 = P1a = P1b_n = P1b_ok = P2a_n = P2a_ok = P2b_n = P2b_ok = 0
    P1a_n = 0
    c3max_bad = []
    v_in = v_n = vneg_in = vneg_n = 0
    P4_n = P4_ok = 0
    P5_bad = []
    allrecs = 0
    desc = []
    for r in recs:
        om = r["omega"]
        for L in r["levels"]:
            allrecs += 1
            P0 += L["gen_in_perp"]
            ex_peel = L["dim_perp"] - om - L["n_peel"]
            ex_occ = L["dim_perp"] - om - L["n_occ"]["0"]
            desc.append((r["cls"], r["bits"], r["seed"], L["level"], ex_peel, ex_occ, L["ck"]["1"],
                         L["ck"]["3"], L["n_peel"], L["sign_in_perp"], L["const_real"]))
            if ex_peel < 0:
                P1a += 1
            if L["level"] >= 1.5:
                P1a_n += 1
                if ex_peel > 4:
                    c3max_bad.append(("P1a", r["cls"], r["bits"], r["seed"], L["level"], ex_peel))
                P2a_n += 1
                P2a_ok += L["ck"]["3"] == 0
                if L["ck"]["3"] > 2:
                    c3max_bad.append(("P2a", r["cls"], r["bits"], r["seed"], L["level"], L["ck"]["3"]))
                for v in L.get("vallee", []):
                    if v["inP3"]:
                        v_n += 1
                        v_in += v["inW"]
                        if v["neg"] and not L["sign_in_perp"]:
                            vneg_n += 1
                            vneg_in += v["inW"]
            if L["level"] >= 2.0:
                P1b_n += 1
                P1b_ok += ex_peel == 0
                P2b_n += 1
                P2b_ok += L["ck"]["1"] == 0
            P4_n += 1
            P4_ok += L["sign_in_perp"] == (1 in L["const_real"])
            if "delta" in L:
                d = L["delta"]
                vals = [d["liouville_joint"], d["liouville_rat"], d["liouville_alg"]] + d["dense"]
                if not all(0.45 <= x <= 0.55 for x in vals):
                    P5_bad.append((r["cls"], r["bits"], r["seed"], [round(x, 3) for x in vals]))
    print(f"level-records: {allrecs}")
    print(f"P0  W ⊆ Λ/2: {P0}/{allrecs}")
    print(f"P1  excess_peel < 0 (sanity): {P1a}; levels>=1.5 records {P1a_n}; "
          f"excess_peel=0 at levels>=2.0: {P1b_ok}/{P1b_n}")
    print(f"P2  c3=0 at levels>=1.5: {P2a_ok}/{P2a_n};  c1=0 at levels>=2.0: {P2b_ok}/{P2b_n}")
    print("    violations of the 'all records' bounds:", c3max_bad)
    print(f"P3  small-square vectors supported on P3 in W: {v_in}/{v_n}; "
          f"negative z (no sign obstruction): {vneg_in}/{vneg_n}")
    print(f"P4  sign obstruction <=> constant-1 real column: {P4_ok}/{P4_n}")
    print(f"P5  delta out of [0.45,0.55]: {P5_bad}")
    print("descriptive (cls bits seed level ex_peel ex_occ c1 c3 n_peel sign_in_perp const_real):")
    for d in sorted(desc):
        print("   ", d)
    tops = [L for r in recs for L in r["levels"] if "delta" in L]
    sp = [x for L in tops for x in L["delta"]["sparse"]]
    if sp:
        print(f"sparse single-prime delta: min {min(sp):.4f} max {max(sp):.4f}")
    print("levels reached:", sorted({(r['bits'], max(L['level'] for L in r['levels'])) for r in recs}))
    print("secs:", [round(r["secs"]) for r in recs], "accidental gcd:", sum(r["accidental_gcd"] for r in recs))


if __name__ == "__main__":
    main()
