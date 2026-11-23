"""Export all survivor measurements and the descriptive follow-up for the strongest hit."""

import gzip, itertools, math
from fractions import Fraction
from collections import defaultdict
import gmpy2 as g
from grammar import gate_specs, PROJ
from statistics_run import *


def describe():
    with open(ROOT / "gates.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["index", "gate", "depth"])
        for i, (name, d) in enumerate(gate_specs()):
            w.writerow([i, name, d])
    with open(ROOT / "targets.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["target_id", "name", "n", "definition"])
        for n in SIZES:
            defs = [f"bit{j}(p)" for j in (n // 2 - 2, n // 4, 1, 2, 3)]
            defs.extend(f"parity(p & {mask:#x})" for mask in masks_gl(n))
            defs.extend(f"Legendre({ell}|p) < 0" for ell in (3, 5, 7, 11))
            defs.extend(["p mod 5 == 1", "p mod 7 in {1,2,3}"])
            defs.extend(f"bit{j}(floor(8p/q))" for j in range(3))
            for j, (label, form) in enumerate(zip(TARGETS, defs)):
                w.writerow([j, label, n, form])


def export():
    surv = list(csv.DictReader(open(R / "survivors.tsv"), delimiter="\t"))
    with gzip.open(R / "survivor_measurements.tsv.gz", "wt") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "n",
                "split",
                "view",
                "test",
                "target",
                "discovery_direction",
                "m",
                "n00",
                "n01",
                "n10",
                "n11",
            ]
            + FIELDS
            + ["oriented_adv", "oriented_lo", "oriented_hi"]
        )
        for n in SIZES:
            for s in SPLITS:
                data = np.load(R / f"{n}_{s}_counts.npz")
                cs = data["counts"]
                ss = np.load(R / f"{n}_{s}_stats.npz")["stats"]
                for i, row in enumerate(surv):
                    t = int(row["test"])
                    j = int(row["target_id"])
                    k = t * T + j if s == "D" else i
                    direction = int(row["direction"])
                    for v in range(4):
                        a = ss[k, v]
                        c = cs[k, v]
                        av, lo, hi = a[:3]
                        oriented = [av, lo, hi] if direction == 1 else [-av, -hi, -lo]
                        w.writerow(
                            [
                                n,
                                s,
                                VIEWS[v],
                                t,
                                TARGETS[j],
                                direction,
                                c.sum(),
                                *c,
                                *a,
                                *oriented,
                            ]
                        )


def exact24():
    # Descriptive population calculation for the first strongest hit. There is no selection on this calculation.
    ps = [x for x in range(2048, 4096) if g.is_prime(x)]
    strata = defaultdict(lambda: np.zeros((2, 2), dtype="i8"))
    for p, q in itertools.combinations(ps, 2):
        N = p * q
        if N.bit_length() != 24:
            continue
        x = 1 - ((N % 16).bit_count() % 2)
        y = (p >> 3) & 1
        strata[N >> 16][x, y] += 1
    cells = [[Fraction(0) for y in range(2)] for x in range(2)]
    for c in strata.values():
        for x in range(2):
            for y in range(2):
                cells[x][y] += Fraction(
                    int(c[x, y]) * int(c[:, 1 - y].sum()), int(c.sum())
                )
    total = sum(sum(v) for v in cells)
    prob = np.array([[float(k / total) for k in rr] for rr in cells])
    adv = (cells[0][0] + cells[1][1]) / total - Fraction(1, 2)
    mi = np.sum(prob * np.log2(prob / (prob.sum(1)[:, None] * prob.sum(0))))
    output = dict(
        test=386,
        expression=CAT[386]["expression"],
        direction=-1,
        target="low3",
        population=sum(int(c.sum()) for c in strata.values()),
        p_discordant=float(total / sum(int(c.sum()) for c in strata.values())),
        matched_adv=float(adv),
        matched_mi_bits=float(mi),
    )
    json.dump(output, open(R / "exact24.json", "w"), indent=2)
    print(output)


if __name__ == "__main__":
    describe()
    export()
    exact24()
