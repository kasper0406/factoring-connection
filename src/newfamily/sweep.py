"""HARNESS: precision sweep for the Prop SAL decoder (Salie oracle at N).  Writes
experiments/analysis/newfamily_salie_sweep.json.  Run: PYTHONPATH=. .venv/bin/python -m src.newfamily.sweep"""
import json

from src.common.instances import balanced_semiprime
from src.newfamily import diag


def main():
    res = []
    for n in (64, 128, 256):
        for delta in (1e-2, 3e-3, 1e-3, 1e-4, 1e-6):
            ok, T, qs = 0, 20, 0
            for seed in range(T):
                inst = balanced_semiprime(n, 24_500_000_000_000 + 7_000_000 + seed)
                p, q = inst.labels["factors"]
                orc = diag.SalieOracle(p, q, delta=delta, seed=seed)
                g = diag.salie_factor(p * q, orc, c=7 + seed, delta=delta)
                ok += g in (p, q)
                qs = orc.queries
            r = {"n": n, "delta": delta, "success": ok, "trials": T, "queries": qs}
            print(json.dumps(r), flush=True)
            res.append(r)
    with open("experiments/analysis/newfamily_salie_sweep.json", "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
