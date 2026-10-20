"""ptesplit pre-registered runs (notes/ptesplit_prereg.md).  Usage: python -m src.ptesplit.run_experiment E1|E2"""
from __future__ import annotations

import json
import resource
import sys
import time

from flint import fmpz_poly

from . import census, nccheck

CFG = json.load(open("experiments/configs/ptesplit.json"))


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20  # macOS: bytes


def E1():
    with open("experiments/raw/ptesplit_E1.jsonl", "w") as fh:
        for M, R in CFG["E1_cells"]:
            t = time.time()
            for row in census.census_cell(M, R, ks=tuple(CFG["E1_ks"])):
                row["secs"] = round(time.time() - t, 1)
                row["rss_mb"] = round(rss_mb(), 1)
                fh.write(json.dumps(row) + "\n")
                fh.flush()
                print({k: row[k] for k in ("M", "R", "k", "pairs", "box_pairs", "max_fibre", "injection_bound")},
                      row.get("kleiman_violations"), row.get("min_constant"), row["rss_mb"], flush=True)


def E2():
    with open("experiments/raw/ptesplit_E2.jsonl", "w") as fh:
        for Fc in CFG["E2_F"]:
            F = fmpz_poly(Fc)
            for M, B in CFG["E2_cells"]:
                for h in CFG["E2_h"]:
                    row = nccheck.nc_row(F, M, B, h)
                    fh.write(json.dumps(row) + "\n")
                    fh.flush()
                    print(row, flush=True)


if __name__ == "__main__":
    {"E1": E1, "E2": E2}[sys.argv[1]]()
