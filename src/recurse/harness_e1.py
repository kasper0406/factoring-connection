"""HARNESS for experiment recurse E1 (notes/recurse_prereg.md).  Knows the factorisations (labels) and
uses sympy.factorint on inner split inputs ONLY to label omega (diagnostic).  The operational pipeline
(src/recurse/pipeline_rec.py) receives only N, the public config and a seed.

Usage: .venv/bin/python -m src.recurse.harness_e1 experiments/configs/recurse_e1.json experiments/raw/recurse_e1.jsonl
"""
from __future__ import annotations

import json
import sys

import sympy

from src.thmA.harness_e1 import make_instance
from . import pipeline_rec as pr


def main(cfg_path: str, out_path: str):
    cfg = json.load(open(cfg_path))
    with open(out_path, "w") as fh:
        for spec in cfg["instances"]:
            N, facs = make_instance(spec)
            rec = {"spec": spec, "N": N, "bits": N.bit_length(), "true_factors": facs, "runs": {}}
            trs = {}
            for o in cfg["oracles"]:
                try:
                    res = pr.factor(N, cfg["pipeline"], spec["seed"], o)
                    trs[o] = pr.outer_transcript(res)
                    splits = []
                    for s in res["stats"].get("splits", []):
                        om = len(sympy.factorint(s["n"]))           # HARNESS-side label
                        splits.append({"depth": s["depth"], "bits": s["n"].bit_length(), "omega": om,
                                       "kind": s["factor_kind"],
                                       "fruitful": [d["fruitful"] for d in s["deps"]],
                                       "row_viol": s["ctr"].get("ROW_IDENTITY_VIOLATION", 0),
                                       "shape_tests": s["ctr"].get("shape_tests", 0),
                                       "seconds": s["seconds"]})
                    rec["runs"][o] = {"ok": res["factors"] == facs, "seconds": res["seconds"],
                                      "max_depth": res["stats"].get("max_depth", 0), "splits": splits}
                except Exception as ex:
                    rec["runs"][o] = {"ok": False, "error": repr(ex)}
            rec["outer_identical"] = (len(trs) == 2 and trs["pb"] == trs["rec"])
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            print(spec, {o: (r.get("ok"), r.get("seconds"), r.get("max_depth")) for o, r in rec["runs"].items()},
                  rec["outer_identical"], flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
