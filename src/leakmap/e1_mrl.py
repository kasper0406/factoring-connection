"""E1 addendum (v2, disclosed deviation): P1 for the registered MR ladder (exact exponents (N-1)/2^j, j <= v2(N-1)),
which v1 replaced by the fixed floor battery MR2.  DIAGNOSTIC."""
import json

from . import netrees as nt
from .e1 import population

if __name__ == "__main__":
    res = []
    for (b, kind, sample) in [(12, "iid", 0), (14, "iid", 0), (16, "iid", 200000), (12, "fermat", 0),
                              (14, "fermat", 0), (16, "fermat", 0)]:
        pairs = population(b, kind, sample, seed=f"leakmapE1|{2 * b}")
        mism = splits = 0
        for p, q in pairs:
            N = p * q
            tree = nt.tree_mr_ladder(N)
            out = nt.run_on_N(tree, N)
            tp, tq = nt.diag_transcript(tree, N, p), nt.diag_transcript(tree, N, q)
            if out[0] == "split":
                splits += 1
                mism += int(tp == tq or out[1] not in (p, q))
            else:
                mism += int(not (tp == tq == out[1]))
        r = dict(b=b, population=kind, pairs=len(pairs), p1_mismatch=mism, split_rate=splits / len(pairs))
        print(r, flush=True)
        res.append(r)
    json.dump(res, open("experiments/analysis/leakmap_e1_mrladder_p1.json", "w"), indent=1)
