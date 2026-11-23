import ctypes as ct
import csv, random, sys, json, time
import numpy as np
from common import *
from grammar import catalog, OPS
from scipy.stats import binom, beta
from scipy.special import xlogy

TARGETS = (
    ["upper", "middle", "low1_chi4", "low2", "low3"]
    + [f"gl{i}" for i in range(4)]
    + [f"leg{e}_neg" for e in (3, 5, 7, 11)]
    + ["mod5_eq1", "mod7_123"]
    + [f"ratio{i}" for i in range(3)]
)
VIEWS = ["raw", "raw_clean", "matched", "matched_clean"]
T = len(TARGETS)
CAT = catalog()
SPECS = np.array(
    [(OPS.index(i["op"]), int(i["a"]), int(i["b"])) for i in CAT], dtype="i4"
)
FAMILY = len(CAT) * T * len(SIZES) * 2
THRESH = 0.05 / FAMILY


def masks_gl(n):
    rng = random.Random(seed("GL", n))
    h = n // 2
    ms = []
    while len(ms) < 4:
        m = rng.getrandbits(h - 2) << 1
        if m.bit_count() >= 2 and m not in ms:
            ms.append(m)
    return ms


def labels(p, q, n):
    bits = [(p >> j) & 1 for j in (n // 2 - 2, n // 4, 1, 2, 3)]
    bits += [np.bitwise_count(p & m) & 1 for m in masks_gl(n)]
    for ell in (3, 5, 7, 11):
        neg = np.array([pow(j, (ell - 1) // 2, ell) == ell - 1 for j in range(ell)])
        if ell % 4 == 3:  # reciprocity for (ell|p)
            bits.append(neg[p % ell] ^ ((p % 4) == 3))
        else:
            bits.append(neg[p % ell])
    bits += [p % 5 == 1, np.isin(p % 7, [1, 2, 3])]
    ratio = 8 * p // q
    bits += [(ratio >> i) & 1 for i in range(3)]
    return np.array(bits, dtype=bool)


def prefix(p, q, n):
    return np.fromiter(
        (int(x) * int(y) >> (n - 8) for x, y in zip(p, q)), dtype=np.uint8, count=len(p)
    )


def pack(x):
    raw = np.packbits(x, axis=-1, bitorder="little")
    pad = (-raw.shape[-1]) % 8
    if pad:
        raw = np.pad(raw, [(0, 0)] * (raw.ndim - 1) + [(0, pad)])
    return np.ascontiguousarray(raw).view(np.uint64)


def unaries(n, s):
    a = np.memmap(str(stem(n, s)) + ".g", dtype="<u2", mode="r", shape=(SPLITS[s], 50))
    return np.ascontiguousarray(
        [pack((a[:, i] >> j) & 1) for op, i, j in SPECS if op == 0]
    )


def context(n, s):
    p, q, r, t = rows(n, s).T
    y = labels(p, q, n)
    yr = labels(r, t, n)
    pref = prefix(p, q, n)
    assert np.array_equal(pref, prefix(r, t, n))
    flags = np.fromfile(str(stem(n, s)) + ".f", dtype="u1").reshape(-1, 2)
    ok = flags[:, 0] == 0
    pair = ok & (flags[:, 1] == 0)
    keep = np.stack(
        [np.ones_like(y), np.broadcast_to(ok, y.shape), y != yr, (y != yr) & pair],
        axis=1,
    )
    return y, keep, pref


def all_tasks():
    return np.column_stack(
        (np.repeat(np.arange(len(CAT)), T), np.tile(np.arange(T), len(CAT)))
    ).astype("i4")


def tally(tasks, unary, y, keep):
    # Axes: target, view, primary label, packed row. The two counts are X=1 at Y=0,1.
    if not len(tasks):
        return np.zeros((0, keep.shape[1], 4), dtype="u4")
    masks = np.stack((keep & ~y[:, None, :], keep & y[:, None, :]), axis=2)
    margins = masks.sum(axis=-1, dtype=np.uint32)
    packed = pack(masks)
    tasks = np.ascontiguousarray(tasks, dtype="i4")
    out = np.empty((len(tasks), keep.shape[1], 2), dtype="u4")
    fun = ct.CDLL(str(ROOT / "score.dylib")).measure
    fun.restype = None
    ptr = lambda a: ct.c_void_p(a.ctypes.data)
    fun(
        unary.shape[-1],
        len(y),
        keep.shape[1],
        len(tasks),
        ptr(tasks),
        ptr(SPECS),
        ptr(unary),
        ptr(packed),
        ptr(out),
    )
    return np.concatenate([margins[tasks[:, 1]] - out, out], axis=-1)


def make_controls(n, s, unary, y, keep, pref):
    ids = [
        next(int(z["id"]) for z in CAT if z["expression"] == f"bit{i}(N)")
        for i in (1, 2, 3)
    ]
    pairs = np.array([(i, j) for i in ids for j in (2, 3, 4)], dtype="i4")
    counts = tally(
        pairs,
        unary,
        y,
        keep[:, 3:4],
    )
    # Permutations are within-prefix and use a separate stream on each split and size.
    ys = y[[1]].copy()
    clean = keep[1, 3]
    rng = np.random.Generator(np.random.PCG64(seed("shuffle", n, s)))
    for b in np.unique(pref[clean]):
        ix = np.flatnonzero(clean & (pref == b))
        ys[0, ix] = rng.permutation(ys[0, ix])
    ts = np.arange(len(CAT)) if s == "D" else np.array(ids[:1])
    shuffle = tally(
        np.column_stack((ts, np.zeros(len(ts), dtype="i4"))),
        unary,
        ys,
        clean[None, None, :],
    )
    first = ids[0]
    raw = np.unpackbits(unary[first].view("u1"), bitorder="little")[: len(pref)].astype(
        bool
    )
    margins = []
    for b in np.unique(pref[clean]):
        ix = clean & (pref == b)
        margins.append(
            [
                int(ix.sum()),
                int(raw[ix].sum()),
                int(ys[0, ix].sum()),
                int((raw[ix] & ys[0, ix]).sum()),
            ]
        )
    return pairs, counts[:, 0], ts, shuffle[:, 0], margins


def run(n, s, tasks):
    t = time.monotonic()
    u = unaries(n, s)
    y, keep, pref = context(n, s)
    counts = tally(tasks, u, y, keep)
    pairs, low, ids, shuf, marg = make_controls(n, s, u, y, keep, pref)
    out = ROOT / "results" / f"{n}_{s}"
    np.savez_compressed(str(out) + "_counts.npz", pairs=tasks, counts=counts)
    np.savez_compressed(
        str(out) + "_controls.npz",
        low_pairs=pairs,
        low=low,
        shuffle_ids=ids,
        shuffle=shuf,
        shuffle_margins=np.array(marg, dtype="i4").reshape(-1, 4),
    )
    print(n, s, len(tasks), round(time.monotonic() - t, 2), flush=True)


def agreement_p(c):
    total = c.sum(axis=-1)
    hit = c[..., 0] + c[..., 3]
    return np.minimum(1, 2 * binom.cdf(np.minimum(hit, total - hit), total, 0.5))


def select():
    tasks = all_tasks()
    best = np.ones(len(tasks))
    bn = np.zeros(len(tasks), "i4")
    bv = bn.copy()
    sign = bn.copy()
    for n in SIZES:
        d = np.load(ROOT / "results" / f"{n}_D_counts.npz")
        assert np.array_equal(tasks, d["pairs"])
        c = d["counts"]
        for v in (2, 3):
            a = c[:, v]
            pv = agreement_p(a)
            sel = pv < best
            # All observed p-values are representable; prevent ambiguous direction ties on extensions.
            if np.any(pv == 0):
                raise ArithmeticError("exact tail underflow in direction selection")
            best[sel] = pv[sel]
            bn[sel] = n
            bv[sel] = v
            sign[sel] = np.where(
                2 * (a[sel, 0] + a[sel, 3]) >= a[sel].sum(axis=1), 1, -1
            )
    survive = np.flatnonzero(best <= THRESH)
    with open(ROOT / "results" / "survivors.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["test", "target_id", "target", "direction", "n", "view", "p_discovery"]
        )
        for i in survive:
            t, j = tasks[i]
            w.writerow([t, j, TARGETS[j], sign[i], bn[i], VIEWS[bv[i]], best[i]])
    print(
        "Discovery survivors",
        len(survive),
        "family",
        FAMILY,
        "cutoff",
        THRESH,
        flush=True,
    )
    return tasks[survive]


if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor

    def batch(s, tasks):
        with ThreadPoolExecutor(3) as pool:
            list(pool.map(lambda n: run(n, s, tasks), SIZES))

    batch("D", all_tasks())
    tasks = select()
    for s in ("V", "C"):
        batch(s, tasks)
