"""Experiment G1 (pre-registered in notes/genus2rh_prereg.md): exact check of
Proposition A (separation) and Lemma 3 (reduction) on toy N = pq.

HARNESS code: it knows p, q. The operational routine `one_curve` sees only N, v, w, rng.
Run:  .venv/bin/python -m src.genus2rh.experiment_g1
"""
import json
import os
import random
from collections import Counter

from sympy import primerange

from .cantor import Jac, NonUnit
from .diagnostic import jac_order, largest_prime, order_of, squarefree_mod
from .hecm import exponent_schedule, monic_model, one_curve

SEED, NPAIRS, NCURVES, V, W, LO, HI = 20261018, 40, 30, 251, 50, 101, 251
OUT = "experiments/raw/genus2rh_g1.jsonl"


def smooth(n, w):
    from sympy import factorint
    return max(factorint(n)) <= w if n > 1 else True


def reduce_D(D, r):
    u = [c % r for c in D[0]]
    v = [c % r for c in D[1]]
    while v and v[-1] == 0:
        v.pop()
    return (tuple(u), tuple(v))


def main():
    rng = random.Random(SEED)
    primes = list(primerange(LO, HI + 1))
    sched = exponent_schedule(V, W)
    k = 1
    for r in sched:
        k *= r
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    tally = Counter()
    with open(OUT, "w") as fh:
        for pair in range(NPAIRS):
            p, q = sorted(rng.sample(primes, 2))
            N = p * q
            crng = random.Random(SEED * 1000 + pair)
            for c in range(NCURVES):
                rec = {}
                d, stage = one_curve(N, V, W, crng, sched, record=rec)
                f, a, b = rec["f"], rec["a"], rec["b"]
                info = {"p": p, "q": q, "curve": c, "found": d, "stage": stage}
                data = {}
                for r in (p, q):
                    fr = [x % r for x in f]
                    ok = fr[5] != 0
                    if ok:
                        F, A, B = monic_model(r, fr, a % r, b % r)
                        ok = squarefree_mod(r, F)
                    if not ok:
                        data[r] = None
                        continue
                    n = jac_order(r, F)
                    D = ((-A % r, 1), (B,) if B else ())
                    o = order_of(r, F, D, n)
                    data[r] = {"F": F, "D": D, "n": n, "ord": o,
                               "smooth": smooth(n, W), "ell": largest_prime(o),
                               "ord_divides_k": k % o == 0}
                covered = False
                for r, s in ((p, q), (q, p)):
                    if data[r] and data[s] and data[r]["smooth"] and \
                            data[s]["n"] % data[r]["ell"] != 0:
                        covered = True
                info["covered"] = covered
                info["J"] = {str(r): (data[r]["n"] if data[r] else None) for r in (p, q)}
                info["ord_divides_k"] = {str(r): (data[r]["ord_divides_k"] if data[r] else None)
                                         for r in (p, q)}
                tally["runs"] += 1
                tally["covered"] += covered
                if d is not None:
                    tally["found_" + stage] += 1
                    if covered:
                        tally["covered_found"] += 1
                    elif stage == "main":
                        tally["uncovered_found_main"] += 1
                    if stage == "main" and data[p] and data[q] and \
                            not data[p]["ord_divides_k"] and not data[q]["ord_divides_k"]:
                        tally["main_found_neither_killed"] += 1
                else:
                    if covered:
                        tally["P1_violation"] += 1
                    if stage is None:
                        # P2: reduction lemma -- recompute over F_p and F_q
                        mism = False
                        for r in (p, q):
                            if data[r] is None:
                                continue
                            Jr = Jac(r, data[r]["F"])
                            Dr = data[r]["D"]
                            for t in sched:
                                Dr = Jr.smul(t, Dr)
                            if reduce_D(rec["final"], r) != Dr:
                                mism = True
                        tally["P2_checked"] += 1
                        tally["P2_mismatch"] += mism
                        info["P2_mismatch"] = mism
                fh.write(json.dumps(info) + "\n")
    print(json.dumps(dict(tally), indent=1, sort_keys=True))
    return tally


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# G1b: POST-HOC supplement (disclosed in notes/genus2rh_results.md). In G1 every run
# found a divisor (accidental intermediate coincidences dominate at p ~ 200), so P2 was
# never exercised. G1b checks only the reduction lemma (P2) at p, q ~ 2^30, where no
# point counting is needed: rerun the same schedule over F_p and F_q and compare.
OUTB = "experiments/raw/genus2rh_g1b.jsonl"


def main_b(npairs=20, ncurves=10, seed=SEED + 1, v=2 ** 30, w=50):
    from sympy import randprime
    rng = random.Random(seed)
    sched = exponent_schedule(v, w)
    tally = Counter()
    with open(OUTB, "w") as fh:
        for pair in range(npairs):
            p = randprime(2 ** 29, 2 ** 30)
            q = randprime(2 ** 29, 2 ** 30)
            while q == p:
                q = randprime(2 ** 29, 2 ** 30)
            N = p * q
            crng = random.Random(seed * 1000 + pair)
            for c in range(ncurves):
                rec = {}
                d, stage = one_curve(N, v, w, crng, sched, record=rec)
                tally["runs"] += 1
                info = {"p": p, "q": q, "curve": c, "found": d, "stage": stage}
                if d is not None:
                    tally["found_" + stage] += 1
                elif stage is None:
                    mism = False
                    for r in (p, q):
                        F = [x % r for x in rec["F"]]
                        A, B = rec["A"] % r, rec["B"] % r
                        Jr = Jac(r, F)
                        Dr = ((-A % r, 1), (B,) if B else ())
                        for t in sched:
                            Dr = Jr.smul(t, Dr)
                        if reduce_D(rec["final"], r) != Dr:
                            mism = True
                    tally["P2_checked"] += 1
                    tally["P2_mismatch"] += mism
                    info["P2_mismatch"] = mism
                    info["final_identity"] = Jac.is_identity(rec["final"])
                    tally["final_identity"] += info["final_identity"]
                fh.write(json.dumps(info) + "\n")
    print(json.dumps(dict(tally), indent=1, sort_keys=True))
    return tally
