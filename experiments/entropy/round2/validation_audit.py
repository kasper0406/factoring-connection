from decay import rowsfile
from inference import *
from scipy.stats import rankdata


def floorids(arr, n):
    N = arr[:, 0] * arr[:, 1]
    off = N - (9 << (n - 4))
    top = (off >> (n - 16)) // 7
    return top.astype("i4")


def loglift(train, test, i, j, dk, vk):
    # Discovery-only Jeffreys base rates; conditional tables shrink toward the stratum base at strength 8.
    def ids(d, k):
        a = d.a[k]
        p, q = a.T
        top = floorids(a, d.n) // 4
        co = p % 256 * (q % 256) % 256 * 64 + top // 16
        fine = p % 3840 * (q % 3840) % 3840 * 1024 + top
        return fine.astype("i4"), co.astype("i4")

    dk = dk & train.ok[j]
    vk = vk & test.ok[j]
    f, c = ids(train, dk)
    fv, cv = ids(test, vk)
    x = train.x(i)[dk]
    y = train.y[j, dk]
    xx = test.x(i)[vk]
    yy = test.y[j, vk]

    def estimate(f, x, y, fv, xx):
        keys, inv, cnt = np.unique(f, return_inverse=True, return_counts=True)
        where = np.searchsorted(keys, fv).clip(0, len(keys) - 1)
        seen = keys[where] == fv
        totals = np.bincount(inv, weights=y, minlength=len(keys))
        floor = (totals + 0.5) / (cnt + 1)
        c1 = np.bincount(inv * 2 + x, minlength=2 * len(keys))
        t1 = np.bincount(inv * 2 + x, weights=y, minlength=2 * len(keys))
        pred = (t1 + 8 * np.repeat(floor, 2)) / (c1 + 8)
        return seen, cnt[where], floor[where], pred[2 * where + xx]

    if not len(fv) or not len(f):
        return np.nan
    one, n1, b1, m1 = estimate(f, x, y, fv, xx)
    two, n2, b2, m2 = estimate(c, x, y, cv, xx)
    globalb = (y.sum() + 0.5) / (len(y) + 1)
    base = np.where(one & (n1 >= 8), b1, np.where(two, b2, globalb))
    model = np.where(one & (n1 >= 8), m1, np.where(two, m2, globalb))
    return np.mean(np.where(yy, np.log(model / base), np.log((1 - model) / (1 - base))))


def main():
    promoted = rowsfile(RES / "promoted.tsv")
    checks = []
    cal = []
    usage = []
    ll = []
    grids = []
    sanity = []
    for n in (40, 48, 56):
        d = Dataset(n, "D")
        v = Dataset(n, "V")
        both = []
        for s, dset in [("D", d), ("V", v)]:
            a = dset.a
            uni, c = np.unique(a, return_counts=True)
            both.append(uni)
            assert len(np.unique(np.ascontiguousarray(a).view("V16"))) == len(a)
            usage.append(
                dict(
                    n=n,
                    split=s,
                    rows=len(a),
                    distinct_products=len(a),
                    distinct_primes=len(uni),
                    mean_endpoint_uses=c.mean(),
                    max_endpoint_uses=c.max(),
                )
            )
            pos = dset.ok[0] & ((dset.budget & 48) == 0)
            vals = np.load(DATA / f"{n}_{s}.labels.npz")["pit"][pos] / 65536
            for k in (1, 2, 3, 4):
                cal.append(
                    dict(
                        n=n,
                        split=s,
                        metric=f"PIT_bit{k}_mean",
                        value=dset.y[k - 1, pos].mean(),
                        m=int(pos.sum()),
                    )
                )
            raw = (a[pos, 0] >> (n // 2 - 3)) & 1
            cal.append(
                dict(
                    n=n,
                    split=s,
                    metric="raw_sorted_first_unfixed_bit_mean",
                    value=raw.mean(),
                    m=len(raw),
                )
            )
            pub = (a[:, 0] * a[:, 1]) % 256
            for k in (1, 2, 3, 4):
                x = (pub >> k) & 1
                j = BASE.index(f"low{k}")
                ix, st = strata(a, n, dset.clean & dset.ok[j])
                s0 = centered(x[ix], dset.y[j, ix], a[ix], st)
                sanity.append(
                    dict(
                        n=n,
                        split=s,
                        floor_bit=k,
                        prime_bit=k,
                        m=s0["m"],
                        MH_variance=s0["mh"],
                        centered=s0["adv"],
                    )
                )
                assert s0["mh"] == 0 and s0["U"] == 0
        assert not len(np.intersect1d(*both))
        trainmask = d.ok[0] & ((d.budget & 48) == 0)
        valmask = v.ok[0] & ((v.budget & 48) == 0)
        f = floorids(d.a[trainmask], n)
        fv = floorids(v.a[valmask], n)
        p = d.a[trainmask, 0]
        pv = v.a[valmask, 0]
        up = np.load(DATA / f"{n}_V.labels.npz")["pit"][valmask] / 65536
        emp = np.empty(len(pv))
        coarse = [np.sort(p[f // 64 == k]) for k in range(64)]
        train_ind = np.argsort(f, kind="stable")
        starts = np.searchsorted(f[train_ind], np.arange(4097))
        vsort = np.argsort(fv, kind="stable")
        vs = np.searchsorted(fv[vsort], np.arange(4097))
        for k in range(4096):
            at = vsort[vs[k] : vs[k + 1]]
            ref = np.sort(p[train_ind[starts[k] : starts[k + 1]]])
            ref = ref if len(ref) >= 8 else coarse[k // 64]
            if not len(at):
                continue
            emp[at] = (np.searchsorted(ref, pv[at], side="left") + 0.5) / (len(ref) + 1)
            grids.append(
                dict(
                    n=n,
                    bin=k,
                    D_count=starts[k + 1] - starts[k],
                    V_count=len(at),
                    PIT_mean=up[at].mean(),
                    empirical_mean=emp[at].mean(),
                )
            )
        for name, val in [
            ("PIT_empirical_rms", np.sqrt(np.mean((emp - up) ** 2))),
            (
                "PIT_empirical_midpoint_disagreement",
                np.mean((emp >= 0.5) != (up >= 0.5)),
            ),
        ]:
            cal.append(dict(n=n, split="V", metric=name, value=val, m=len(emp)))
        for row in promoted:
            t = int(row["test"])
            j = int(row["target"])
            sg = int(row["direction"])
            for kind in ("inclusive", "full_clean", "path_clean"):

                def keep(ds):
                    if kind == "inclusive":
                        return (ds.budget & 48) == 0
                    if kind == "full_clean":
                        return ds.clean
                    ok = ds.budget == 0
                    for a in ds.parents(t):
                        ok &= (ds.f[a][:, 4] & 1) == 0
                    return ok

                kd = keep(d)
                kv = keep(v)
                for fine in (False, True):
                    z = fullstat(v, t, j, kv, fine)
                    checks.append(
                        dict(
                            n=n,
                            test=t,
                            target=TARGETS[j],
                            view=kind,
                            fine=int(fine),
                            sign=sg,
                            **z,
                        )
                    )
                if kind == "path_clean":
                    ll.append(
                        dict(
                            n=n,
                            test=t,
                            target=TARGETS[j],
                            view=kind,
                            validation_loglift_nats_per_row=loglift(d, v, t, j, kd, kv),
                        )
                    )
    for name, table in [
        ("validation_views", checks),
        ("validation_loglift", ll),
        ("pit_calibration", cal),
        ("pit_empirical_by_bin", grids),
        ("discovery_prime_reuse", usage),
        ("floor_sanity", sanity),
    ]:
        write_records(name + ".tsv", table)


if __name__ == "__main__":
    main()
