from inference import *
import ctypes as ct


def nodes(n):
    while (
        not (RES / f"identities{n}.json").is_file()
        or (RES / f"identities{n}.json").stat().st_size < 50
    ):
        time.sleep(2)
    run(["./nodes", n, DATA / f"identities{n}.list", DATA / f"id{n}"])
    print(n, "nodes ready", flush=True)


def estimate(n, cohort):
    while not (RES / f"{cohort}{n}_measure.tsv").exists():
        time.sleep(2)
    uniq, _ = plans(n)
    pcount = (DATA / f"id{n}.repeated_primes").stat().st_size // 8
    pairs = np.fromfile(DATA / f"id{n}.repeated_pairs", "<u8").reshape(-1, 2)
    # Store the exact repeated-edge correction only for rows with recurrent pair identities.
    duplicate = dict(((int(p), int(q)), i) for i, (p, q) in enumerate(pairs))
    fingerprint = lambda a: (a[:, 0] * np.uint64(0x9E3779B185EBCA87)) ^ a[:, 1]
    dup_h = fingerprint(pairs)
    fun = ct.CDLL(str(R / "dyad.dylib")).accumulate
    fun.restype = ct.c_double
    ptr = lambda a: ct.c_void_p(a.ctypes.data)
    keys = [k for k in uniq if not k.startswith("plant") and k != "negative"]
    out = []
    for start in range(0, len(keys), 4):
        group = keys[start : start + 4]
        cache = {}
        for k, mats, fiber in tables(n, cohort):
            if k not in group:
                continue
            counts = mats[1]
            cs = counts.sum(1)
            large = cs >= 8
            extra = reduce_bins(counts, large, COARSE)[int(large.sum()) :]
            m = counts.sum()

            def moment(c):
                count = c.sum(1)
                div = np.maximum(1, count)
                px = (c[:, 2] + c[:, 3]) / div
                py = (c[:, 1] + c[:, 3]) / div
                return np.column_stack(
                    (px * py, -px * (1 - py), -(1 - px) * py, (1 - px) * (1 - py))
                )

            lut = moment(counts)
            lut[~large] = moment(extra)[COARSE[~large]]
            avg = np.sum(counts * lut) / m
            lut -= avg
            cache[k] = dict(
                lut=lut,
                m=m,
                s2=np.sum(counts * lut**2),
                delta=0.0,
                edge=0.0,
                sums=np.zeros(pcount + 1),
                edges=np.zeros(len(pairs)),
            )
        for phase in ("C", "S") if cohort == "augmented" and n in (64, 96) else ("C",):
            fs = (
                jobs(n)
                if phase == "C"
                else [
                    DATA / f"S{n}_{i:03d}.pairs" for i in range(128 if n == 64 else 170)
                ]
            )
            sub = {
                k: uniq[k]
                for k in group
                if not (phase == "S" and n == 96 and uniq[k][0] >= 4)
            }
            if not sub:
                continue
            for fn in fs:
                a, h, b, cols = coded(n, fn, sub)
                node = np.fromfile(str(fn) + ".nodes", "<u4").reshape(-1, 2)
                sparse = n > 64
                if len(pairs):
                    di = np.flatnonzero(np.isin(fingerprint(a), dup_h))
                    di = [(i, duplicate.get((int(a[i, 0]), int(a[i, 1])))) for i in di]
                else:
                    di = []
                for k, x, y, good, clean, _ in cols:
                    c = cache[k]
                    e = c["lut"][h, 2 * x.astype("u1") + y]
                    e[~clean] = 0
                    c["delta"] += fun(
                        ct.c_int64(len(node)),
                        ptr(node),
                        ptr(e),
                        ptr(c["sums"]),
                        ct.c_int(int(sparse)),
                    )
                    for i, j in di:
                        if j is None:
                            continue
                        c["edge"] += 2 * c["edges"][j] * e[i]
                        c["edges"][j] += e[i]
        for k, c in cache.items():
            v = max(0, c["s2"] + c["delta"] - c["edge"])
            out.append(
                dict(
                    cohort=cohort,
                    n=n,
                    key=k,
                    m=c["m"],
                    dyad_se=2 * np.sqrt(v) / c["m"],
                    iid_sandwich_se=2 * np.sqrt(c["s2"]) / c["m"],
                    design_effect=v / c["s2"] if c["s2"] > 0 else np.nan,
                )
            )
    write_records(f"{cohort}{n}_dyad.tsv", out)
    print(cohort, n, "dyad", flush=True)


if __name__ == "__main__":
    for n in (40, 48, 56, 80, 64, 96):
        nodes(n)
        estimate(n, "original")
        if n in (64, 96):
            estimate(n, "augmented")
