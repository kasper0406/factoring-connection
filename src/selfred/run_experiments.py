"""Run E1 / E2 of notes/selfred_prereg.md.  Usage:
   .venv/bin/python -m src.selfred.run_experiments e1 14 16 18
   .venv/bin/python -m src.selfred.run_experiments e2
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy

from . import census as C

RAW = Path(__file__).resolve().parents[2] / "experiments" / "raw"


def run_e1(bs):
    out = RAW / "selfred_e1.jsonl"
    with out.open("w") as f:
        for b in bs:
            t = time.time()
            for r in C.e1_census(b):
                f.write(json.dumps(r) + "\n")
            print(f"E1 b={b} done in {time.time() - t:.1f}s", flush=True)


def run_e2():
    out = RAW / "selfred_e2.jsonl"
    ps = []
    for start in (10 ** 4, 10 ** 5):
        p = start
        for _ in range(40):
            p = int(sympy.nextprime(p))
            ps.append(p)
    with out.open("w") as f:
        for D in (1, 2, 3, 5, 8):
            for k in range(8):
                h = C.e2_polys(D, k)
                U = C.unit_fixed_divisor(h)
                for p in ps:
                    r = C.e2_case(h, p)
                    r.update(D0=D, k=k, U=U, P2_ok=(r["Gr"] == U // p ** C.vl(U, p)))
                    f.write(json.dumps(r) + "\n")
            print(f"E2 D={D} done", flush=True)


def run_e1_posthoc(bs):
    out = RAW / "selfred_e1_posthoc.jsonl"
    with out.open("w") as f:
        for b in bs:
            for r in C.e1_posthoc_dedup(b):
                f.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "e1posthoc":
        run_e1_posthoc([int(x) for x in sys.argv[2:]])
    elif sys.argv[1] == "e1":
        run_e1([int(x) for x in sys.argv[2:]])
    else:
        run_e2()
