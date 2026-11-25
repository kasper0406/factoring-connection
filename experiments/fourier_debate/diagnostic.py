"""Post hoc, truth-using diagnostic, separate from all timed search results.

On all forms with exactly one matched witness, check its eight ancestors. Compare
the two child scores, then replace that witness's pair count by the unconditional
primitive-point expectation, retaining both child populations and every other
point's contribution. Score ties at 1/2. Pool first and last four split accuracies
and bootstrap whole forms within bit size, holding the singleton counts fixed.
"""
import json
from collections import Counter

import numpy as np

from common import HERE, K, SIZES, exclusion, hatch, primerange
from guide import box_sum, children, pair_tables


def preference(delta):
    return 1 if delta > 0 else 0 if delta < 0 else 0.5


def scores(row):
    form = row["form"]
    primes = list(primerange(row["B"]+1))
    wt, pop, num, den, _ = pair_tables(form, primes, False, Counter())
    (u, v), = row["matched_points"]
    own = box_sum(wt, (u-1,v-1,u,v))
    value = hatch.norm(form,u,v)
    ex = exclusion(form)
    good = sum(value % p == 0 and ex % p != 0 for p in primes)
    assert own == good*(good-1)//2
    box = (0,0,K,K)
    res = []
    for _ in range(8):
        left, right = children(box)
        x,y,X,Y = left
        target, other = (left,right) if x < u <= X and y < v <= Y else (right,left)
        delta = (box_sum(wt,target)-box_sum(wt,other))*den
        delta -= (box_sum(pop,target)-box_sum(pop,other))*num
        res.append([preference(delta), preference(delta-own*den+num)])
        box = target
    return np.array(res)


def main():
    rows = [json.loads(line) for line in (HERE/"heldout.jsonl").read_text().splitlines()]
    rng = np.random.default_rng(202611253)
    obs = np.zeros((8,2))
    boots = np.zeros((10000,8,2))
    n = 0
    for bits in SIZES:
        vals = np.array([scores(row) for row in rows if row["bits"]==bits and len(row["matched_points"])==1])
        m = len(vals)
        n += m
        obs += vals.sum(axis=0)
        W = rng.multinomial(m, np.full(m, 1/m), size=len(boots))
        boots += (W @ vals.reshape(m,16)).reshape(-1,8,2)
    obs /= n
    boots /= n
    result = dict(singleton_forms=n)
    for key, group in (("upper_four", slice(0,4)), ("lower_four",slice(4,8))):
        p = obs[group].mean(axis=0)
        q = np.quantile(boots[:,group].mean(axis=1), [.025,.975], axis=0)
        result[key] = {mode:[float(p[i]), float(q[0,i]), float(q[1,i])] for i,mode in enumerate(["original","replaced"])}
    print(json.dumps(result))
    (HERE/"diagnostic.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
