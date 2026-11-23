from common import *
import csv, math, random, subprocess, json
from collections import Counter, defaultdict


def gate_specs():
    gates = [("N", 0)]
    gates += [
        (f"pow({x},{e},N)", 1)
        for x in (2, 3, 5, 7)
        for e in ("N", "(N-1)/2", "N-1", "N+1")
    ]
    gates += [(f"LucasV(N,{p},{q}) mod N", 1) for p, q in ((3, 1), (5, 1), (3, -1))]
    gates += [(f"pow(pow({x},N,N){c:+d},N,N)", 2) for x in (2, 3) for c in (-1, 1)]
    for k in (3, 5, 7, 11, 13, 16, 31):
        gates += [(f"floor(N/{k})", 1), (f"N mod {k}", 1)]
    gates += [
        ("isqrt(N)", 1),
        ("frac16(sqrt(N))", 1),
        ("frac16(N sqrt(2))", 1),
        ("frac16(N sqrt(3))", 1),
    ]
    gates += [(f"cf({a})[{k}]", 1) for a in ("sqrt(N)", "N/2^h") for k in range(1, 5)]
    return gates


PROJ = [f"bit{i}" for i in range(8)] + [
    "parity8",
    "parity16",
    "parity",
    "jacobi_neg",
    "jacobi_zero",
]
OPS = ["unary", "xor", "and", "or"]


def catalog():
    return list(csv.DictReader(open(ROOT / "catalog.tsv"), delimiter="\t"))


def witnesses():
    blocks = []
    for n in SIZES:
        rng = random.Random(seed("witness", n))
        raw = []
        while len(raw) < 64:
            x = rng.getrandbits(n) | (1 << (n - 1)) | 1
            if math.isqrt(x) ** 2 != x:
                raw.append([x % (1 << 64), x >> 64, 0, 0])
        st = str(ROOT / "data" / f"witness{n}")
        np.array(raw, dtype="<u8").tofile(st + ".rows")
        subprocess.run(
            [
                str(ROOT / "gates"),
                str(n),
                "witness",
                st + ".rows",
                st + ".g",
                st + ".f",
            ],
            check=True,
        )
        blocks.append(np.fromfile(st + ".g", dtype="<u2").reshape(64, -1))
    return np.concatenate(blocks)


def build():
    a = witnesses()
    assert a.shape == (512, 50)
    rng = random.Random(seed("grammar"))
    full = (1 << len(a)) - 1
    seen = {0, full}
    unary = []
    items = []
    for g, (name, dep) in enumerate(gate_specs()):
        for b, proj in enumerate(PROJ):
            key = int.from_bytes(
                np.packbits((a[:, g] >> b) & 1, bitorder="little").tobytes(), "little"
            )
            if key in seen:
                continue
            seen.add(key)
            idx = len(unary)
            unary.append((dep + 1, key, f"{proj}({name})"))
            items.append((idx, dep + 1, "unary", g, b, unary[-1][2]))
    buckets = defaultdict(list)
    for i, (d, _, _) in enumerate(unary):
        for j in range(i):
            dep = 1 + max(d, unary[j][0])
            for op in (1, 2, 3):
                buckets[(dep, op)].append((j, i))
    for pairs in buckets.values():
        rng.shuffle(pairs)
    while len(items) < 12000:
        before = len(items)
        for (dep, op), pairs in sorted(buckets.items()):
            while pairs:
                i, j = pairs.pop()
                x = unary[i][1]
                y = unary[j][1]
                key = (x ^ y) if op == 1 else ((x & y) if op == 2 else (x | y))
                if key in seen:
                    continue
                seen.add(key)
                items.append((len(items), dep, OPS[op], i, j, f"{OPS[op]}(t{i},t{j})"))
                break
            if len(items) == 12000:
                break
        if len(items) == before:
            raise RuntimeError("exhausted")
    with open(ROOT / "catalog.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "depth", "op", "a", "b", "expression"])
        w.writerows(items)
    print(
        json.dumps(
            {
                "total": len(items),
                "unary": len(unary),
                "by_depth": dict(Counter(i[1] for i in items)),
                "distinct_witness_vectors": len(seen) - 2,
            }
        )
    )


if __name__ == "__main__":
    build()
