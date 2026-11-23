from config import *
from labels import TARGETS, BASE
from grammar import tab
import ctypes as ct, csv, json
from scipy.stats import binom, norm
from collections import Counter

CAT = list(csv.DictReader(open(RES / "catalog.tsv"), delimiter="\t"))
VARS = json.loads((RES / "catalog.json").read_text())["eligible"]
OPS = np.array(
    [
        ({"unary": 0, "xor": 1, "maj3": 2}[r["op"]], *[int(r[k]) for k in "abc"])
        for r in CAT
    ],
    dtype="i4",
)
NU = int(np.sum(OPS[:, 0] == 0))
T = len(TARGETS)
FAMILY = len(CAT) * T * 4
TH = 0.05 / FAMILY


def pack(a):
    x = np.packbits(a, axis=-1, bitorder="little")
    pad = (-x.shape[-1]) % 8
    if pad:
        x = np.pad(x, [(0, 0)] * (x.ndim - 1) + [(0, pad)])
    return np.ascontiguousarray(x).view("<u8")


def strata(a, n, mask, fine=False):
    p, q = a.T
    N = p * q
    if (
        n > 64
    ):  # exact high part through Python integers in the small locked-candidate path
        pp = np.fromiter(
            (
                (int(x) * int(y) - (9 << (n - 4))) >> (n - 4 - (12 if fine else 10))
                for x, y in a
            ),
            dtype="u8",
        )
    else:
        pp = (N - (9 << (n - 4))) >> (n - 4 - (12 if fine else 10))
    top = pp // 7
    mod = 430080 if fine else 3840
    residue = p % mod * (q % mod) % mod
    full = (residue * (4096 if fine else 1024) + top).astype("i8")
    coarse = (p % 256 * (q % 256) % 256 * 64 + (top >> (6 if fine else 4))).astype("i8")
    ix = np.flatnonzero(mask)
    keys, inv, counts = np.unique(full[ix], return_inverse=True, return_counts=True)
    compact = np.where(counts[inv] >= 8, full[ix], -1 - coarse[ix])
    return ix, compact


def pair_bins(a, n, s, mask, fine=False):
    ix, key = strata(a, n, mask, fine)
    rng = np.random.default_rng(seed("matched", n, s, int(fine)))
    order = rng.permutation(len(ix))
    order = order[np.argsort(key[order], kind="stable")]
    sortedkey = key[order]
    if len(order) < 2:
        return np.empty(0, "i4"), np.empty(0, "i4")
    starts = np.r_[0, np.flatnonzero(np.diff(sortedkey)) + 1]
    sizes = np.diff(np.r_[starts, len(order)])
    pos = np.arange(len(order)) - np.repeat(starts, sizes)
    left = np.flatnonzero((pos % 2 == 0) & (pos + 1 < np.repeat(sizes, sizes)))
    return ix[order[left]], ix[order[left + 1]]


class Dataset:
    def __init__(self, n, s):
        self.n = n
        self.s = s
        self.a = load(n, s)
        self.m = len(self.a)
        self.f = {}
        for j in range(0, len(VARS), 3):
            group = VARS[j : j + 3]
            block = np.memmap(
                DATA / f"{n}_{s}_v{j}.g",
                mode="r",
                dtype="<u2",
                shape=(len(group), self.m, 5),
            )
            for k, v in enumerate(group):
                self.f[v] = block[k]
        self.budget = np.fromfile(str(path(n, s)) + ".budget", dtype="u1")
        self.flags = np.zeros(self.m, "u1")
        self.reason = {}
        for v, m in self.f.items():
            flags = m[:, 4] & 1
            self.flags |= flags.astype("u1")
            self.reason[v] = int(flags.sum())
        self.clean = (self.flags == 0) & (self.budget == 0)
        with np.load(DATA / f"{n}_{s}.labels.npz") as f:
            self.y = np.unpackbits(f["y"], axis=1, count=self.m).astype(bool)
            self.ok = np.unpackbits(f["ok"], axis=1, count=self.m).astype(bool)

    def x(self, t, ix=slice(None)):
        op, a, b, c = OPS[t]
        if op == 0:
            return ((self.f[a][ix, b] >> c) & 1).astype(bool)
        l = self.x(a, ix)
        r = self.x(b, ix)
        if op == 1:
            return l ^ r
        third = self.x(c, ix)
        return (l & r) | (l & third) | (r & third)

    def parents(self, t):
        op, a, b, c = OPS[t]
        return (
            {a}
            if op == 0
            else set.union(
                *(self.parents(j) for j in ([a, b] if op == 1 else [a, b, c]))
            )
        )

    def matched_screen(self):
        l, r = pair_bins(self.a, self.n, self.s, self.clean)
        if not len(l):
            return (
                np.zeros((len(CAT), T + 1, 2), "u4"),
                np.zeros(len(CAT), "u4"),
                np.zeros(T + 1, "u4"),
                0,
            )
        X = np.ascontiguousarray(
            [[pack(self.x(t, ix)) for t in range(NU)] for ix in (l, r)]
        )
        ok = self.ok[:, l] & self.ok[:, r]
        yl = self.y[:, l]
        yr = self.y[:, r]
        yp = yl & ~yr & ok
        ym = yr & ~yl & ok
        # Fixed shuffled sentinel uses a separately randomized within-pair label exchange.
        mid = TARGETS.index("min/midhalf")
        swap = (
            np.random.default_rng(seed("shuffle", self.n))
            .integers(0, 2, len(l), dtype="u1")
            .astype(bool)
        )
        yp = np.vstack((yp, np.where(swap, ym[mid], yp[mid])))
        ym = np.vstack((ym, np.where(swap, (yl[mid] & ~yr[mid] & ok[mid]), ym[mid])))
        Y = pack(np.stack((yp, ym), axis=1))
        c = np.empty((len(CAT), T + 1, 2), "u4")
        ones = np.empty(len(CAT), "u4")
        m = np.r_[ok.sum(1), ok[mid].sum()]
        fn = ct.CDLL(str(R / "matched.dylib")).matched
        fn.restype = None
        ptr = lambda a: ct.c_void_p(a.ctypes.data)
        fn(
            X.shape[-1],
            len(CAT),
            T + 1,
            NU,
            ptr(OPS),
            ptr(X),
            ptr(Y),
            ptr(c),
            ptr(ones),
        )
        return c, ones, m, len(l)


def pv(c):
    m = c.sum(axis=-1, dtype="i8")
    return np.minimum(1.0, 2 * binom.cdf(c.min(axis=-1), m, 0.5))


def marginal_pass(ones, pairs):
    return (
        (ones >= 0.0001 * 2 * pairs) & (ones <= (1 - 0.0001) * 2 * pairs) & (pairs > 0)
    )


def one(n):
    d = Dataset(n, "D")
    c, o, m, p = d.matched_screen()
    np.savez_compressed(
        RES / f"D{n}_screen.npz", counts=c, ones=o, eligible_pairs=m, matched_pairs=p
    )
    rec = dict(
        n=n,
        rows=d.m,
        clean=int(d.clean.sum()),
        matched_pairs=p,
        budget_clean=int((d.budget == 0).sum()),
        gate_hits=int(d.flags.sum()),
        by_variant=d.reason,
    )
    (RES / f"D{n}_flags.json").write_text(json.dumps(rec))
    print(n, rec["clean"], p, flush=True)


def centered(x, y, arr, stratum):
    if not len(x):
        return dict(m=0, adv=np.nan, Z=0, U=0, var=0, mh=0, dyad=0)
    _, g, c = np.unique(stratum, return_inverse=True, return_counts=True)
    xx = np.bincount(g, weights=x, minlength=len(c))
    yy = np.bincount(g, weights=y, minlength=len(c))
    e = (x - xx[g] / c[g]) * (y - yy[g] / c[g])
    U = e.sum()
    V = np.sum((xx * (1 - xx / c)) * (yy * (1 - yy / c)) / np.maximum(c - 1, 1))
    _, f = np.unique(arr.reshape(-1), return_inverse=True)
    res = e - e.mean()
    cl = np.bincount(f.ravel(), weights=np.repeat(res, 2))
    dyad = max(0.0, cl @ cl - res @ res)
    v = max(V, dyad) if arr.max() < (1 << 24) else V
    return dict(
        m=len(x),
        adv=2 * U / len(x),
        Z=U / np.sqrt(v) if v > 0 else 0,
        U=U,
        var=v,
        mh=V,
        dyad=dyad,
    )


def fullstat(d, t, j, clean, fine=False):
    keep = clean & d.ok[j]
    ix, h = strata(d.a, d.n, keep, fine)
    return centered(d.x(t, ix), d.y[j, ix], d.a[ix], h)


def discover():
    tables = []
    passes = []
    data = []
    per = []
    small = []
    for n in (40, 48, 56):
        z = np.load(RES / f"D{n}_screen.npz")
        c = z["counts"][:, :T, :]
        tables.append(c)
        keep = marginal_pass(z["ones"], int(z["matched_pairs"]))
        passes.append(keep)
        p = pv(c)
        p[~keep] = 1
        per.append(p)
        crop = np.flatnonzero(p.ravel() <= TH)
        small.extend((n, int(i // T), int(i % T), float(p.ravel()[i])) for i in crop)
    pooled = np.sum(tables, axis=0)
    pp = pv(pooled)
    pp[~np.any(passes, axis=0)] = 1
    per.append(pp)
    small.extend(
        ("pool", int(i // T), int(i % T), float(pp.ravel()[i]))
        for i in np.flatnonzero(pp.ravel() <= TH)
    )
    candidate = sorted({(i, j) for _, i, j, _ in small})
    print("nominal candidates", len(candidate), "cells", len(small), flush=True)
    cache = {}
    promote = []
    details = []
    if candidate:
        data = [Dataset(n, "D") for n in (40, 48, 56)]
    for i, j in candidate:
        stats = [fullstat(d, i, j, d.clean) for d in data]
        cells = []
        for z in range(4):
            p = float(per[z][i, j])
            tc = (tables[z] if z < 3 else pooled)[i, j]
            sign = 1 if tc[0] >= tc[1] else -1
            if z == 3:
                v = sum(st["var"] for st in stats)
                Z = sum(st["U"] for st in stats) / np.sqrt(v) if v else 0.0
            else:
                Z = stats[z]["Z"]
            prob = max(p, 2 * norm.sf(sign * Z))
            cells.append((prob, z, sign, p, Z))
            details.append([i, j, z, p, sign, Z, prob])
        best = min(cells)
        if (
            best[0] <= TH
            and not BASE[j % len(BASE)].startswith("low")
            and not (OPS[i, 0] == 0 and OPS[i, 1] == 42)
        ):
            promote.append((best[0], i, j, best[2], best[1]))
    promote = sorted(promote)[:64]
    tab(
        RES / "discovery_sensitivity.tsv",
        ["test", "target", "screen", "p_sign", "sign", "z_max_variance", "p_gate"],
        details,
    )
    tab(
        RES / "promoted.tsv",
        ["p_gate", "test", "target", "direction", "screen"],
        promote,
    )
    rec = {
        "family": FAMILY,
        "cutoff": TH,
        "nominal_pairs": len(candidate),
        "nominal_cells": len(small),
        "promoted": len(promote),
        "min_p": min(float(p.min()) for p in per),
        "shuffled_discoveries": {
            n: int(
                (
                    pv(np.load(RES / f"D{n}_screen.npz")["counts"][:, -1])
                    <= 0.05 / (len(CAT) * 3)
                ).sum()
            )
            for n in (40, 48, 56)
        },
    }
    (RES / "discovery.json").write_text(json.dumps(rec, indent=2))
    print(rec, flush=True)
    return promote


def validation(promote):
    if not promote:
        tab(RES / "locks.tsv", ["rank", "test", "target", "direction", "holm_p"], [])
        tab(
            RES / "validation.tsv",
            ["test", "target", "direction", "clean_sign_p", "holm_p", "fine_direction"],
            [],
        )
        return
    data = [Dataset(n, "V") for n in (40, 48, 56)]
    acc = []
    for _, t, j, direction, _ in promote:
        pos = neg = 0
        fineU = 0
        for d in data:
            l, r = pair_bins(d.a, d.n, d.s, d.clean)
            ok = d.ok[j, l] & d.ok[j, r]
            X = d.x(t).astype("i1")
            Y = d.y[j].astype("i1")
            z = (X[l] - X[r]) * (Y[l] - Y[r])
            pos += np.sum(ok & (z == direction))
            neg += np.sum(ok & (z == -direction))
            fineU += fullstat(d, t, j, d.clean, True)["U"]
        p = binom.sf(pos - 1, pos + neg, 0.5) if pos + neg else 1.0
        acc.append((p, t, j, direction, fineU * direction))
    acc.sort()
    adj = np.maximum.accumulate(
        [min(1.0, x[0] * (len(acc) - i)) for i, x in enumerate(acc)]
    )
    tab(
        RES / "validation.tsv",
        ["test", "target", "direction", "clean_sign_p", "holm_p", "fine_direction"],
        [(t, j, d, p, a, u) for (p, t, j, d, u), a in zip(acc, adj)],
    )
    locked = [
        (t, j, d, a) for (p, t, j, d, u), a in zip(acc, adj) if a <= 0.01 and u > 0
    ][:20]
    tab(
        RES / "locks.tsv",
        ["rank", "test", "target", "direction", "holm_p"],
        [(i + 1, *v) for i, v in enumerate(locked)],
    )


if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(2) as ex:
        list(ex.map(one, (40, 48, 56)))
    validation(discover())
