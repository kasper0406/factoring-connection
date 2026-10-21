"""E0b extension (MODEL; disclosed): larger n for 800/1024 bits, where the E0b grid (n <= 1200) was too short.
Usage: .venv/bin/python -m src.schnorrcvp.run_e0b_ext -> experiments/raw/schnorrcvp_E0b_ext.jsonl"""
import json

import mpmath as mp

from .core import first_primes
from .psi import ideal_enumeration_cost
from .run_e0b import Lcrude

if __name__ == "__main__":
    out = open("experiments/raw/schnorrcvp_E0b_ext.jsonl", "w")
    for bits, ns in [(800, [1600, 2000, 2400, 3200]), (1024, [1600, 2400, 3200, 4000])]:
        for n in ns:
            r = ideal_enumeration_cost(first_primes(n), bits * mp.log(2))
            lv = r.get("log2V")
            T = None if lv is None else 0.292 * n + max(0.0, lv - 0.2075 * n)
            rec = {"bits": bits, "n": n, "log2V": lv, "log2T_sieve": T, "lnX": r.get("lnX"),
                   "log2_QS": Lcrude(bits, 0.5, 1.0), "log2_GNFS": Lcrude(bits, 1 / 3, (64 / 9) ** (1 / 3))}
            out.write(json.dumps(rec) + "\n"); out.flush(); print(rec, flush=True)
