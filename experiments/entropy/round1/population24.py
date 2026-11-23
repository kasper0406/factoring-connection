"""Enumerate the finite 24-bit population at the registered grammar and eight-bit matching rule."""

import itertools, math, subprocess
import gmpy2 as g
from statistics_run import *


def main():
    ps = [p for p in range(2048, 4096) if g.is_prime(p)]
    factors = np.array(
        [(p, q) for p, q in itertools.combinations(ps, 2) if p * q >= 2**23],
        dtype=np.uint64,
    )
    p, q = factors.T
    N = p * q
    size = len(N)
    # All members of the complete population split within the inclusive Fermat window.
    assert all(
        (int(a) + int(b)) // 2 - (1 + math.isqrt(int(a) * int(b) - 1)) < 256
        for a, b in factors
    )
    gate_input = np.column_stack((factors, factors)).astype("<u8")
    st = str(ROOT / "data" / "population24")
    gate_input.tofile(st + ".rows")
    subprocess.run(
        [str(ROOT / "gates"), "24", "data", st + ".rows", st + ".g", st + ".f"],
        check=True,
    )
    assert np.all(np.fromfile(st + ".f", dtype="u1") & 8)
    a = np.fromfile(st + ".g", dtype="<u2").reshape(size, 50)
    unary = np.array([((a[:, i] >> j) & 1) != 0 for op, i, j in SPECS if op == 0])
    tasks = all_tasks()
    targ = tasks[:, 1]
    ys = labels(p, q, 24)
    raw = np.zeros((len(tasks), 4), dtype="i4")
    matched = np.zeros((len(tasks), 4))
    pref = N >> 16
    for b in np.unique(pref):
        ix = np.flatnonzero(pref == b)
        y = ys[:, ix]
        c = tally(tasks, pack(unary[:, ix]), y, np.ones((T, 1, len(ix)), dtype=bool))[
            :, 0, :
        ]
        raw += c.astype("i4")
        # Exact population integration: a primary in this bin with label 0 or 1
        # survives matching with probability #opposite labels / #rows in the bin.
        freq = y.sum(1)[targ] / len(ix)
        for cell in range(4):
            matched[:, cell] += c[:, cell] * (1 - freq if cell % 2 else freq)
    joint = np.divide(
        matched,
        matched.sum(1)[:, None],
        out=np.full_like(matched, np.nan),
        where=matched.sum(1)[:, None] != 0,
    )
    adv = joint[:, 0] + joint[:, 3] - 0.5
    mi = np.maximum(
        0,
        hbinary(joint[:, 0] + joint[:, 1])
        + hbinary(joint[:, 0] + joint[:, 2])
        - hpart(joint).sum(1),
    )
    np.savez_compressed(
        R / "population24.npz",
        pairs=tasks,
        raw_counts=raw,
        matched_joint=joint,
        matched_adv=adv,
        matched_mi_bits=mi,
        discordance=matched.sum(1) / size,
    )
    surv = list(csv.DictReader(open(R / "survivors.tsv"), delimiter="\t"))
    oriented = np.array(
        [
            int(d["direction"]) * adv[int(d["test"]) * T + int(d["target_id"])]
            for d in surv
        ]
    )
    example = 386 * T + 4
    assert np.isclose(
        -adv[example],
        json.load(open(R / "exact24.json"))["matched_adv"],
        rtol=0,
        atol=1e-14,
    )
    result = dict(
        population=size,
        grammar_test_target_pairs=len(tasks),
        discovery_survivors=len(surv),
        positive_in_frozen_direction=int((oriented > 0).sum()),
        min_oriented_adv=float(oriented.min()),
        max_oriented_adv=float(oriented.max()),
        fermat_splits=size,
    )
    json.dump(result, open(R / "population24_summary.json", "w"), indent=2)
    print(result)


if __name__ == "__main__":
    main()
