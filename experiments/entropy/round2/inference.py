"""Confirmation estimates. Candidate identities, signs and the two cohorts stay explicit."""

from aggregate import *
from scipy.special import betaincinv, xlogy
from scipy.stats import chi2
import warnings


def write_records(name, rows):
    with open(RES / name, "w") as f:
        if rows:
            w = csv.DictWriter(
                f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            w.writeheader()
            w.writerows(rows)


def cp(k, m, alpha=0.05):
    k = np.asarray(k, dtype=float)
    lo = np.zeros(k.shape)
    hi = np.ones(k.shape)
    z = k > 0
    lo[z] = betaincinv(k[z], m - k[z] + 1, alpha / 2)
    z = k < m
    hi[z] = 1 - betaincinv(m - k[z], k[z] + 1, alpha / 2)
    return lo, hi


def mass_support(key):
    if "_y" not in key:
        return np.ones(H, bool)
    name = BASE[int(key.split("_y")[1]) % len(BASE)]
    if name.startswith("gen"):
        d = int(name[3:])
        if d == -1:
            keep = RESIDUES % 4 == 1
        elif d == 2:
            keep = RESIDUES % 8 == 1
        elif 3840 % abs(d) == 0:
            l = abs(d)
            keep = np.array([pow(int(v % l), (l - 1) // 2, l) == 1 for v in RESIDUES])
        else:
            keep = np.ones(1024, bool)
        return np.repeat(keep, 1024)
    return np.ones(H, bool)


def reduce_bins(mat, large, coarse):
    sparse = ~large
    extra = np.column_stack(
        [
            np.bincount(coarse[sparse], weights=mat[sparse, j], minlength=8192)
            for j in range(4)
        ]
    )
    return np.vstack((mat[large], extra))


def miller(table):
    t = np.atleast_2d(table).astype("f8")
    n = t.sum(1)
    m = n.sum()
    if not m:
        return np.nan, np.nan
    x = np.column_stack((t[:, 0] + t[:, 1], t[:, 2] + t[:, 3]))
    y = np.column_stack((t[:, 0] + t[:, 2], t[:, 1] + t[:, 3]))
    # I = (sum c ln c + n ln n - sum x ln x - sum y ln y)/(m ln 2).
    info = (
        xlogy(t, t).sum() + xlogy(n, n).sum() - xlogy(x, x).sum() - xlogy(y, y).sum()
    ) / (m * np.log(2))
    df = (t > 0).sum() - (x > 0).sum() - (y > 0).sum() + (n > 0).sum()
    return info, info - df / (2 * m * np.log(2))


def exact_cov_box(table, support, large, coarse):
    # Full atomic cells are fixed before sampling. Simultaneous coverage survives adaptive pooling.
    m = int(table.sum())
    h = len(table)
    k = np.arange(int(table.max()) + 1)
    lo, hi = cp(k, m, 0.025 / (4 * h))
    a = reduce_bins(lo[table], large, coarse)
    b = reduce_bins(hi[table], large, coarse)
    nlo = a.sum(1)
    nhi = b.sum(1)
    low = 2 * (a[:, 0] * a[:, 3] - b[:, 1] * b[:, 2])
    high = 2 * (b[:, 0] * b[:, 3] - a[:, 1] * a[:, 2])
    with np.errstate(divide="ignore", invalid="ignore"):
        c0 = np.where(low < 0, low / nlo, low / nhi)
        c1 = np.where(high > 0, high / nlo, high / nhi)
    c0 = np.where(np.isnan(c0), -0.5 * nhi, c0)
    c1 = np.where(np.isnan(c1), 0.5 * nhi, c1)
    l = np.maximum(-0.5 * nhi, c0).sum()
    u = np.minimum(0.5 * nhi, c1).sum()
    # Uniform over every grouping: |2*E Cov| <= 2*min(P(X=1),P(X=0)).
    xl, xu = cp(np.array([table[:, 2:].sum()]), m, 0.025)
    cap = min(0.5, 2 * xu[0], 2 * (1 - xl[0]))
    return max(-cap, l), min(cap, u)


def estimated(key, raw, suppfiber=False):
    support = mass_support(key)
    if key == "negative":
        support &= np.repeat(RESIDUES % 4 == 3, 1024)
    if suppfiber:
        support &= np.repeat(RESIDUES % 8 == 1, 1024)
    assert not np.any(raw[~support])
    t = raw[support]
    co = COARSE[support]
    total = t.sum(0, dtype="i8")
    m = int(total.sum())
    if not m:
        return {}
    large = t.sum(1) >= 8
    g = reduce_bins(t, large, co)
    g = g[g.sum(1) > 0]
    c = g.sum(1)
    xx = g[:, 2] + g[:, 3]
    yy = g[:, 1] + g[:, 3]
    px = xx / c
    py = yy / c
    cov = g[:, 3] / c - px * py
    U = (c * cov).sum()
    mh = np.sum(xx * (1 - px) * yy * (1 - py) / np.maximum(c - 1, 1))
    residual = (
        np.column_stack((px * py, -px * (1 - py), -(1 - px) * py, (1 - px) * (1 - py)))
        - cov[:, None]
    )
    se = 2 * np.sqrt((g * residual**2).sum()) / m
    hit = int(total[0] + total[3])
    al, ah = cp(np.array([hit]), m)
    lo, hi = exact_cov_box(t, support, large, co)
    mi, mi_bc = miller(total)
    cmi, cmi_bc = miller(g)
    return dict(
        m=m,
        n00=total[0],
        n01=total[1],
        n10=total[2],
        n11=total[3],
        adv=hit / m - 0.5,
        lo=al[0] - 0.5,
        hi=ah[0] - 0.5,
        centered=2 * U / m,
        centered_exact_lo=lo,
        centered_exact_hi=hi,
        centered_se=se,
        mh_z=U / np.sqrt(mh) if mh > 0 else 0.0,
        mh_se=2 * np.sqrt(mh) / m,
        mi_bits=mi,
        mi_miller_bits=mi_bc,
        conditional_mi_bits=cmi,
        conditional_mi_miller_bits=cmi_bc,
        coarse_fraction=t[~large].sum() / m,
        x_rate=(total[2] + total[3]) / m,
        y_rate=(total[1] + total[3]) / m,
    )


def tables(n, cohort):
    km, names = plans(n)
    arr = np.load(RES / f"C{n}_strata.npy", mmap_mode="r")
    if cohort == "original" or n not in (64, 96):
        for i, k in enumerate(km):
            yield k, arr[i], False
        return
    ext = np.load(RES / f"S{n}_strata.npy", mmap_mode="r")
    sk = list(np.load(RES / f"S{n}_masks.npz")["keys"])
    for i, (k, (rank, *_)) in enumerate(km.items()):
        if k not in sk:
            yield k, arr[i], False
            continue
        if n == 96 and rank < 0:
            if k == "negative":
                yield k, arr[i], False
                continue
            a = arr[i].copy()
            a[:, np.repeat(RESIDUES % 8 != 1, 1024)] = 0
            yield k, a + ext[sk.index(k)], True
        else:
            yield k, arr[i] + ext[sk.index(k)], False


def safe_power(m, delta, alpha):
    k = int(binom.isf(alpha, m, 0.5)) + 1
    return binom.sf(k - 1, m, 0.5 + delta)


def size(n, cohort):
    res = {}
    km, allr = plans(n)
    norg = len(locks(n))
    for k, stats, fiber in tables(n, cohort):
        for v in (0, 1):
            res[k, v] = estimated(k, stats[v], fiber)
    rows = []
    power = []
    for idx, (k, sign, target) in enumerate(allr):
        for v in (0, 1):
            d = res[k, v].copy()
            if not d:
                continue
            m = d["m"]
            o = (d["lo"], d["hi"]) if sign > 0 else (-d["hi"], -d["lo"])
            cl = (
                (d["centered_exact_lo"], d["centered_exact_hi"])
                if sign > 0
                else (-d["centered_exact_hi"], -d["centered_exact_lo"])
            )
            stat = d["mh_z"] * sign if idx < norg else d["adv"] * 2 * np.sqrt(m)
            p = (
                norm.sf(stat)
                if idx < norg
                else binom.sf(int(d["n00"] + d["n11"]) - 1, m, 0.5)
            )
            row = dict(
                cohort=cohort,
                n=n,
                view="clean" if v else "all",
                rank=idx + 1 if idx < norg else 0,
                key=k,
                target=target,
                sign=sign,
                **d,
                signed_agreement=sign * d["adv"],
                signed_lo=o[0],
                signed_hi=o[1],
                signed_centered=sign * d["centered"],
                signed_centered_lo=cl[0],
                signed_centered_hi=cl[1],
                z=stat,
                one_sided_p=p,
                holm_p=np.nan,
            )
            rows.append(row)
            if v:
                delta = {
                    "plant1": 1 / n,
                    "plant16sq": 16 / n**2,
                    "plant1sq": 1 / n**2,
                }.get(k, 1 / n**2)
                power.append(
                    dict(
                        cohort=cohort,
                        n=n,
                        rank=row["rank"],
                        key=k,
                        target=target,
                        m=m,
                        reference_adv=delta,
                        power_single=safe_power(m, delta, 0.01),
                        power_first_holm=safe_power(m, delta, 0.01 / max(1, norg)),
                        observed_agreement=row["signed_agreement"],
                        exact_lo=row["signed_lo"],
                        exact_hi=row["signed_hi"],
                        one_sided_p=p,
                    )
                )
    for v in ("all", "clean"):
        sub = sorted(
            [r for r in rows if r["rank"] and r["view"] == v],
            key=lambda r: r["one_sided_p"],
        )
        cur = 0
        for i, r in enumerate(sub):
            cur = max(cur, min(1.0, r["one_sided_p"] * (len(sub) - i)))
            r["holm_p"] = cur
    write_records(f"{cohort}{n}_measure.tsv", rows)
    write_records(f"{cohort}{n}_power.tsv", power)
    print(cohort, n, len(rows), flush=True)


if __name__ == "__main__":
    import sys, time

    cohort = sys.argv[1]
    for n in (40, 48, 56, 64, 96, 80):
        prefix = "S" if cohort == "augmented" and n in (64, 96) else "C"
        while not (RES / f"{prefix}{n}_masks.npz").exists():
            time.sleep(3)
        size(n, cohort)
