"""systeq experiment runner (pre-registered in notes/systeq_prereg.md).

    .venv/bin/python -m src.systeq.run_experiment [E1|E2|E3|E4|all]

Writes experiments/raw/systeq_<E>.jsonl.  Configuration: experiments/configs/systeq.json.
Fractions are stored exactly as "num/den" strings plus a float for display only.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CFG = json.loads((ROOT / "experiments/configs/systeq.json").read_text())


def _ser(v):
    if isinstance(v, Fraction):
        return {"exact": f"{v.numerator}/{v.denominator}", "float": float(v)}
    if isinstance(v, dict):
        return {k: _ser(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_ser(x) for x in v]
    return v


def _write(name, rows):
    out = ROOT / f"experiments/raw/systeq_{name}.jsonl"
    with out.open("w") as fh:
        for r in rows:
            fh.write(json.dumps(_ser(r)) + "\n")
    print(f"wrote {out}")


def e1():
    from .bgjt import control_summary
    rows = []
    for p, e0 in CFG["E1"]["fields"]:
        t = time.time()
        s = control_summary(p, e0, seed=CFG["seed"])
        s["seconds"] = round(time.time() - t, 2)
        rows.append(s)
        print({k: (float(v) if isinstance(v, Fraction) else v) for k, v in s.items()})
    _write("E1", rows)


def e2():
    from .intanalog import random_baseline, vfamily_bits, vfamily_yield
    c = CFG["E2"]
    B = c["B"]
    rows = []
    for n, m in c["instances"]:
        bits = vfamily_bits(n, m, c["H"], c["count"], CFG["seed"])
        base = random_baseline(bits, B, c["baseline_samples"], CFG["seed"])
        r = vfamily_yield(n, m, c["H"], c["count"], B, base, CFG["seed"])
        rows.append(r)
        print({k: (float(v) if isinstance(v, Fraction) else v) for k, v in r.items() if k != "v_bits"})
    _write("E2", rows)


def e3():
    from .intanalog import descartes_count
    c = CFG["E3"]
    rows = []
    for n in c["n"]:
        r = descartes_count(n, c["H"])
        rows.append(r)
        print(r)
    _write("E3", rows)


def e4():
    from .intanalog import cyclo_bits, cyclo_calibration, random_baseline
    c = CFG["E4"]
    rows = []
    for n, lo, hi in c["instances"]:
        base = random_baseline(cyclo_bits(n, lo, hi), c["B"], c["baseline_samples"], CFG["seed"])
        r = cyclo_calibration(n, lo, hi, c["B"], base)
        rows.append(r)
        print({k: (float(v) if isinstance(v, Fraction) else v) for k, v in r.items()})
    _write("E4", rows)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for name, fn in [("E1", e1), ("E2", e2), ("E3", e3), ("E4", e4)]:
        if which in (name, "all"):
            fn()
