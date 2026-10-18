"""HARNESS runner for pre-registered experiment E1 (notes/nfsfruit_prereg.md).
Usage: .venv/bin/python -m src.nfsfruit.run_e1 experiments/configs/nfsfruit_e1.json OUT.jsonl [workers] [--pilot SEED]"""
import json
import sys
from multiprocessing import Pool

from .harness import run_instance


def _job(args):
    fam, cls, bits, seed, cfg, st, lc = args
    try:
        return run_instance(fam, cls, bits, seed, cfg, subset_trials=st, lin_checks=lc)
    except Exception as e:  # record, never hide
        return dict(family=fam, cls=cls, bits=bits, seed=seed, status="error", error=repr(e))


def main():
    cfgp, outp = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else 8
    C = json.load(open(cfgp))
    seeds = list(range(*C["seeds"]))
    if "--pilot" in sys.argv:
        seeds = [int(sys.argv[sys.argv.index("--pilot") + 1])]
    jobs = []
    for fam, per in C["families"].items():
        for bits in C["bits"]:
            for cls in C["classes"]:
                for sd in seeds:
                    jobs.append((fam, cls, bits, sd, per[str(bits)], C["subset_trials"], C["lin_checks"]))
    jobs.sort(key=lambda j: (j[0].startswith("nfs3") and j[2] == 44), reverse=True)  # slow first
    with open(outp, "w") as fo, Pool(workers) as pool:
        for rec in pool.imap_unordered(_job, jobs):
            fo.write(json.dumps(rec) + "\n")
            fo.flush()


if __name__ == "__main__":
    main()
