"""Run archive validation, then the locked, interleaved held-out comparison."""
import json
import math
import random
import sys
from collections import Counter
from time import perf_counter_ns, process_time_ns

from common import HERE, K, REPS, SIZES, archived, exclusion, hatch, primerange, truth
from guide import matched, pair_tables, prefix, spatial


def audit_archive(bounds):
    for bits, _, _, form in archived():
        primes = list(primerange(bounds[str(bits)] + 1))
        reference = truth(form, primes)[1]
        _, original_hits, _ = hatch.forced_search(form, primes, False)
        assert original_hits == set(reference)
        pair_data = pair_tables(form, primes, True, Counter())
        inc_data = pair_tables(form, primes, False, Counter())
        assert pair_data == inc_data
        ex = exclusion(form)
        # A good prime dividing a primitive norm has a projective root
        # necessarily. This direct calculation does not call line_roots.
        possible = [p for p in primes if ex % p]
        direct = []
        for v in range(1, K+1):
            for u in range(1, K+1):
                if math.gcd(u, v) > 1:
                    direct.append(0)
                    continue
                value = hatch.norm(form, u, v)
                count = sum(value % p == 0 for p in possible)
                direct.append(math.comb(count, 2))
        assert prefix(direct) == pair_data[0]
    print("18 archive cases validated", flush=True)


def run():
    bounds = json.loads((HERE / "calibration.json").read_text())["bounds"]
    audit_archive(bounds)
    used = {n for _, _, n, _ in archived()}
    dest = HERE / "heldout.jsonl"
    with dest.open("x") as out:
        for rep in range(REPS):
            for bits in SIZES:
                seed = 202611252000 + 10000*bits + rep
                rng = random.Random(seed)
                bound = bounds[str(bits)]
                primes = list(primerange(bound+1))
                try:
                    n = hatch.fixtures.fixture(bits, rng)
                    form = hatch.fixtures.form(n, primes, rng)
                except (ValueError, AssertionError) as exc:
                    print(json.dumps(dict(bits=bits, rep=rep, seed=seed, generation_error=repr(exc))),
                          file=out, flush=True)
                    continue
                assert n not in used
                used.add(n)
                public_rng = random.Random(n ^ 0x124AC)
                tie_coins = public_rng.getrandbits(256)
                order = ["pairs", "incidence", "blind", "matched"]
                public_rng.shuffle(order)
                arms = {}
                audit = {}
                for arm in order:
                    then = perf_counter_ns()
                    cpu_then = process_time_ns()
                    if arm == "matched":
                        hit, counters = matched(form, primes)
                    else:
                        hit, counters = spatial(form, primes, tie_coins, arm, audit)
                    cpu_ns = process_time_ns() - cpu_then
                    wall_ns = perf_counter_ns() - then
                    arms[arm] = dict(wall_ns=wall_ns, cpu_ns=cpu_ns, hit=hit, counts=counters)
                assert audit["pairs"] == audit["incidence"]
                whole, matches = truth(form, primes)
                for arm in arms.values():
                    hit = arm["hit"]
                    assert (hit is not None) == bool(matches)
                    if hit is not None:
                        assert tuple(hit[:2]) in matches
                assert arms["pairs"]["hit"] == arms["incidence"]["hit"]
                assert arms["pairs"]["counts"]["nodes"] == arms["incidence"]["counts"]["nodes"]
                print(json.dumps(dict(bits=bits, rep=rep, seed=seed, N=n, form=form, B=bound,
                                      unrestricted=whole, matched_points=sorted(matches), order=order,
                                      arms=arms)), file=out, flush=True)
            if (rep+1) % 32 == 0:
                print(f"{rep+1} of {REPS} forms per size", flush=True)


if __name__ == "__main__":
    run()
