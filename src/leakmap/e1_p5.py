"""E1 P5 post-processing (DIAGNOSTIC): UNIVERSAL shares by population, and the exact ARCH pin count.

ARCH node k deviates UNIVERSALLY at the smaller prime p iff (p-1) | N - k, i.e. q = p - 1 + k (balanced range),
so the number of pairs with a UNIVERSAL ARCH deviation at some node must equal #{pairs: 1 <= q - p + 1 <= 32}.
"""
import json
from collections import Counter

from . import netrees as nt
from .e1 import population


def arch_pin_check(b, kind, sample=0):
    pairs = population(b, kind, sample, seed=f"leakmapE1|{2 * b}")
    tree = nt.tree_arch()
    nodes = tree[1]
    pins = univ_pairs = 0
    for p, q in pairs:
        N = p * q
        if 1 <= q - p + 1 <= 32:
            pins += 1
        hit = False
        for nd in nodes:
            for x in (p, q):
                if pow(2, N, x) == pow(2, nd.t_target, x) and (N - nd.t_target) % (x - 1) == 0:
                    hit = True
        univ_pairs += hit
    return dict(b=b, population=kind, pairs=len(pairs), pins=pins, universal_pairs=univ_pairs, ok=pins == univ_pairs)


def shares(path="experiments/raw/leakmap_e1.jsonl"):
    out = []
    for l in open(path):
        r = json.loads(l)
        c = Counter()
        for k, v in r["dev_classes"].items():
            c[k.split("|")[0]] += v
        tot = sum(c.values())
        out.append(dict(n=r["n"], population=r["population"], tree=r["tree"], deviations=tot,
                        universal_share=(c["UNIVERSAL"] / tot if tot else None)))
    return out


if __name__ == "__main__":
    res = {"arch_pins": [arch_pin_check(b, k, s) for (b, k, s) in
                         [(12, "iid", 0), (14, "iid", 0), (16, "iid", 200000), (12, "fermat", 0), (14, "fermat", 0),
                          (16, "fermat", 0)]],
           "shares": shares()}
    json.dump(res, open("experiments/analysis/leakmap_e1_p5.json", "w"), indent=1)
    for r in res["arch_pins"]:
        print(r)
    for r in res["shares"]:
        print(r)
