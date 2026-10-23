"""E2 summary (v2): confirmed hits per target with the corrected statistics (singleton fix, dyadic-robust z).
Primary z: n = 32 dyadic-robust (z_cr), n = 40 permutation (z); both reported."""
import json

out = {}
for n in (32, 40):
    m = json.load(open(f"experiments/analysis/leakmap_e2_main_n{n}.json"))
    f = json.load(open(f"experiments/analysis/leakmap_e2_fresh_n{n}.json"))
    key = "z_cr" if n == 32 else "z"
    for t in m["targets"]:
        hits = []
        for st, d in m["targets"][t]["stats"].items():
            if abs(d[key]) >= 5:
                fz = f["targets"][t]["stats"][st][key]
                hits.append(dict(stat=st, z=d[key], z_fresh=fz, confirmed=bool(fz * d[key] > 0 and abs(fz) >= 3.5)))
        mx = max(abs(d[key]) for d in m["targets"][t]["stats"].values())
        mxo = max(abs(d["z_cr" if key == "z" else "z"]) for d in m["targets"][t]["stats"].values())
        mxf = max(abs(d[key]) for d in f["targets"][t]["stats"].values())
        tg = m["targets"][t]
        out[f"{n}_{t}"] = dict(n=n, target=t, rows=tg["n_rows"], hits=len(hits),
                               confirmed=sum(h["confirmed"] for h in hits), max_primary=mx, max_other=mxo,
                               max_fresh=mxf, mde=tg["mde_rho"], mde_eff=tg.get("mde_rho_eff"),
                               deff=tg.get("max_design_effect"), detail=hits)
        print(n, t, "hits", len(hits), "conf", out[f"{n}_{t}"]["confirmed"], "max", round(mx, 2),
              "other", round(mxo, 2), "fresh", round(mxf, 2), "mde", round(tg["mde_rho"], 4),
              "mde_eff", round(tg.get("mde_rho_eff") or 0, 4), "deff", round(tg.get("max_design_effect") or 0, 2))
json.dump(out, open("experiments/analysis/leakmap_e2_summary.json", "w"), indent=1)
