"""HARNESS for c12alg E2/E2b (notes/c12alg_prereg.md).  Knows the factorisations (labels); the operational
pipeline (src/c12alg/dw_pipeline.factor_dw = thmA pipeline with DW fields) receives only N, config and seed.

Usage: .venv/bin/python -m src.c12alg.harness_e2 experiments/configs/c12alg_e2.json experiments/raw/c12alg_e2.jsonl
"""
from __future__ import annotations

import json
import sys
import time

from src.thmA.harness_e1 import make_instance
from . import dw_pipeline


def main(cfg_path: str, out_path: str):
    cfg = json.load(open(cfg_path))
    with open(out_path, "w") as fh:
        for block in cfg["blocks"]:
            pcfg = dict(cfg["pipeline"]); pcfg.update(block.get("pipeline_override", {}))
            for spec in block["instances"]:
                N, facs = make_instance(spec)
                t0 = time.time()
                rec = {"block": block["name"], "spec": spec, "N": N, "bits": N.bit_length()}
                try:
                    res = dw_pipeline.factor_dw(N, pcfg, spec["seed"])
                    rec["ok"] = (res["factors"] == facs)
                    rec["out_factors"] = res["factors"]
                    log = res["log"]
                    for L in log:
                        if L["step"] == "split":
                            n = L["n"]
                            L["omega"] = len({p for p in facs if n % p == 0})
                            L["sum_inv_p"] = sum(1.0 / p for p in set(facs) if n % p == 0)
                    rec["log"] = log
                except Exception as ex:
                    rec["ok"] = False
                    rec["error"] = repr(ex)
                rec["true_factors"] = facs
                rec["seconds"] = round(time.time() - t0, 2)
                fh.write(json.dumps(rec) + "\n"); fh.flush()
                print(block["name"], spec, rec["ok"], rec["seconds"], flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
