import csv, json, random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from config import *

NAMES = (R / "variants.txt").read_text().splitlines()
READS = [f"bit{j}" for j in range(6)] + [
    "parity8",
    "parity",
    "jac_neg",
    "jac_zero",
    "frac_bit3",
    "frac_bit2",
    "frac_bit1",
    "weight_gt_half",
    "negative",
    "dense_parity",
]


def engine(mode, vs, src, dst, **kw):
    return run(["./engine", mode, ",".join(map(str, vs)), src, dst], **kw)


def tab(file, header, rows):
    with open(file, "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def build(frozen=False):
    ids = list(range(len(NAMES)))
    passes = []
    if frozen:
        with open(RES / "timings.tsv") as f:
            times = np.array(
                [float(row["median_us"]) for row in csv.DictReader(f, delimiter="\t")]
            )
        keep = [i for i in ids if times[i] + 0.05 <= 50]
    else:
        for k in range(3):
            outs = engine(
                "time",
                ids,
                DATA / "cost.rows",
                "/dev/null",
                text=True,
                capture_output=True,
            )
            passes.append([float(x.split()[1]) for x in outs.stdout.splitlines()])
        times = np.median(passes, axis=0)
        keep = [i for i in ids if times[i] + 0.05 <= 50]
        tab(
            RES / "timings.tsv",
            ["variant", "name", "us1", "us2", "us3", "median_us", "eligible"],
            [
                (i, NAMES[i], *(t[i] for t in passes), times[i], int(i in keep))
                for i in ids
            ],
        )
    panels = []
    for n in (40, 48, 56, 64, 80, 96):
        fn = DATA / f"wit{n}.g"
        engine("features", keep, DATA / f"wit{n}.rows", fn)
        panels.append(np.fromfile(fn, "<u2").reshape(len(keep), 86, 5))
    truth = np.concatenate(panels, axis=1)[:, :512]
    full = (1 << 512) - 1
    seen = {0}
    un = []
    spec = []
    words = []
    depcost = []

    def include(x):
        key = min(x, full ^ x)
        if key in seen:
            return False
        seen.add(key)
        return True

    for i, v in enumerate(keep):
        for slot in range(4):
            for b in range(16):
                x = int.from_bytes(
                    np.packbits(
                        (truth[i, :, slot] >> b) & 1, bitorder="little"
                    ).tobytes(),
                    "little",
                )
                if include(x):
                    un.append((v, slot, b))
                    words.append(x)
                    depcost.append(times[v] + 0.05)
                    spec.append(
                        [
                            len(spec),
                            "unary",
                            v,
                            slot,
                            b,
                            times[v] + 0.05,
                            f"{NAMES[v]}[{slot}].{READS[b]}",
                        ]
                    )
    rng = random.Random(seed("grammar"))
    nu = len(un)
    while len(spec) < 20000:
        op = "xor" if len(spec) % 2 else "maj3"
        pick = sorted(rng.sample(range(nu), 2 if op == "xor" else 3))
        parents = {un[j][0] for j in pick}
        fam = {NAMES[v].split("/")[0] for v in parents}
        if len(fam) < 2:
            continue
        cost = sum(times[i] for i in parents) + 0.15
        if cost > 50:
            continue
        a, b = pick[:2]
        c = pick[-1] if op == "maj3" else -1
        x = (
            words[a] ^ words[b]
            if c == -1
            else ((words[a] & words[b]) | (words[a] & words[c]) | (words[b] & words[c]))
        )
        if include(x):
            spec.append(
                [
                    len(spec),
                    op,
                    a,
                    b,
                    c,
                    cost,
                    f"{op}(" + ",".join(f"t{j}" for j in pick) + ")",
                ]
            )
    tab(RES / "catalog.tsv", ["id", "op", "a", "b", "c", "cost_us", "expression"], spec)
    rec = dict(
        eligible=keep,
        total=len(spec),
        unary=nu,
        count_by_op=dict(Counter(x[1] for x in spec)),
        witnesses=512,
    )
    (RES / "catalog.json").write_text(json.dumps(rec, indent=2))
    print(rec, flush=True)


if __name__ == "__main__":
    import sys

    build("--frozen" in sys.argv)
