"""HARNESS for experiment thmA E1/E1b (notes/thmA_prereg.md).  Knows the factorisations (labels); the
operational pipeline (src/thmA/pipeline.py) receives only N, the public config and a seed.

Usage: .venv/bin/python -m src.thmA.harness_e1 <config.json> <out.jsonl>
"""
from __future__ import annotations

import json
import random
import sys
import time
from math import prod

import gmpy2

from src.common import instances as I
from . import pipeline


def _rand_prime(rng: random.Random, bits: int) -> int:
    while True:
        x = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
        if gmpy2.is_prime(x):
            return x


def make_instance(spec: dict):
    fam, bits, seed = spec["family"], spec["bits"], spec["seed"]
    if fam == "semi":
        inst = I.balanced_semiprime(bits, seed); return inst.N, inst.factors
    if fam == "tri":
        inst = I.multi_prime(bits, 3, seed); return inst.N, inst.factors
    if fam == "unbal":
        inst = I.unbalanced_semiprime(bits, spec["small_bits"], seed); return inst.N, inst.factors
    if fam == "ppow":
        inst = I.prime_power(bits, spec["k"], seed); return inst.N, inst.factors
    if fam == "prime":
        inst = I.prime_instance(bits, seed); return inst.N, inst.factors
    rng = random.Random(f"thmA_harness|{fam}|{bits}|{seed}")
    if fam == "small_x_semi":
        s = _rand_prime(rng, spec["small_bits"])
        inst = I.balanced_semiprime(bits - spec["small_bits"], seed)
        return s * inst.N, sorted([s] + inst.factors)
    if fam == "p2q":
        while True:
            p, q = _rand_prime(rng, bits // 3), _rand_prime(rng, bits - 2 * (bits // 3))
            if p != q:
                return p * p * q, sorted([p, p, q])
    if fam == "semi_sq":
        inst = I.balanced_semiprime(bits // 2, seed)
        return inst.N ** 2, sorted(inst.factors * 2)
    raise ValueError(fam)


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
                    res = pipeline.factor(N, pcfg, spec["seed"])
                    rec["ok"] = (res["factors"] == facs)
                    rec["out_factors"] = res["factors"]
                    log = res["log"]
                    for L in log:                       # harness-side annotation: omega of the split input
                        if L["step"] == "split":
                            n = L["n"]
                            L["omega"] = len({p for p in facs if n % p == 0})
                            L["sum_inv_p"] = sum(1.0 / p for p in set(facs) if n % p == 0)
                    rec["log"] = log
                except Exception as ex:                  # cap reached or bug: recorded as failure
                    rec["ok"] = False
                    rec["error"] = repr(ex)
                rec["true_factors"] = facs
                rec["seconds"] = round(time.time() - t0, 2)
                fh.write(json.dumps(rec) + "\n"); fh.flush()
                print(block["name"], spec, rec["ok"], rec["seconds"], flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
