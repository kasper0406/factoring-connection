"""HARNESS runner for Addendum A (E1b, E2, E3).  Usage:
.venv/bin/python -m src.nfsfruit.run_add OUT.jsonl [workers]"""
import json
import sys
from multiprocessing import Pool

from .harness2 import run_e1b, run_e2, run_e3

E1CFG = json.load(open("experiments/configs/nfsfruit_e1.json"))
LIN = {"40": dict(B=1500, A=600, Bb=200, excess=60), "44": dict(B=2000, A=800, Bb=250, excess=60)}


def _job(j):
    try:
        kind = j[0]
        if kind == "e1b":
            return run_e1b(*j[1:], orderings=5)
        if kind == "e2":
            return run_e2(*j[1:])
        if kind == "e3":
            return run_e3(*j[1:])
    except Exception as e:
        return dict(exp=j[0], args=str(j[1:4]), status="error", error=repr(e))


def jobs():
    out = []
    for fam_, per in E1CFG["families"].items():
        for bits in (40, 44):
            for cls in ("blum", "mixed", "both1", "tri"):
                for sd in range(20):
                    out.append(("e1b", fam_, cls, bits, sd, per[str(bits)]))
    for bits in (40, 44):
        for sd in range(40):
            out.append(("e2", "qs", bits, sd, E1CFG["families"]["qs"][str(bits)]))
            for f in ("lin_base", "lin_rand"):
                c = dict(E1CFG["families"][f][str(bits)])
                c["excess"] = 700
                out.append(("e2", f, bits, sd, c))
            for cls in ("blum", "mixed", "both1", "tri"):
                out.append(("e3", cls, bits, sd, LIN[str(bits)]))
    out.sort(key=lambda j: j[0] == "e1b" and j[1].startswith("nfs3") and j[3] == 44, reverse=True)
    return out


if __name__ == "__main__":
    outp = sys.argv[1]
    w = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    with open(outp, "w") as fo, Pool(w) as pool:
        for rec in pool.imap_unordered(_job, jobs()):
            fo.write(json.dumps(rec) + "\n")
            fo.flush()
