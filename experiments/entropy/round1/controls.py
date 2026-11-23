import csv
from scipy.stats import binom, hypergeom
from statistics_run import *


def write_rows(path, rows):
    if not rows:
        return
    with open(R / path, "w") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def exact_power(m, p, upper_alpha):
    if not m:
        return np.nan
    cutoff = int(binom.isf(upper_alpha, m, 0.5)) + 1
    return float(binom.sf(cutoff - 1, m, p))


def sentinel(marg):
    if not len(marg):
        return np.nan
    total = np.array([1.0])
    start = 0
    for m, x, y, xy in marg:
        lo = max(0, int(x + y) - int(m))
        hi = min(x, y)
        mass = hypergeom.pmf(np.arange(lo, hi + 1), m, x, y)
        mass /= mass.sum()
        total = np.convolve(total, mass)
        start += lo
    i = marg[:, 3].sum() - start
    return min(1, 2 * min(total[: i + 1].sum(), total[i:].sum()))


def control_tables():
    planted = []
    low = []
    shuffle = []
    for n in SIZES:
        for s in SPLITS:
            y, keep, pref = context(n, s)
            for k in (1, 2):
                rng = np.random.Generator(np.random.PCG64(seed("plant", k, n, s)))
                p = 0.5 + 1 / n**k
                # Construct the actual synthetic bit. The target is the primary middle bit.
                same = rng.random(SPLITS[s]) < p
                test = y[1] ^ ~same
                for v in (0, 2, 3):
                    ix = keep[1, v]
                    m = int(ix.sum())
                    c = np.bincount(
                        (2 * test[ix].astype("i4") + y[1, ix]), minlength=4
                    ).reshape(1, 4)
                    st = for_size(c, m)[0]
                    hit = int(c[0, 0] + c[0, 3])
                    planted.append(
                        dict(
                            n=n,
                            split=s,
                            view=VIEWS[v],
                            k=k,
                            m=m,
                            hits=hit,
                            nominal_adv=p - 0.5,
                            **dict(zip(FIELDS, st)),
                            p_one=binom.sf(hit - 1, m, 0.5) if m else np.nan,
                            p_two=agreement_p(c)[0] if m else np.nan,
                            power_one_05=exact_power(m, p, 0.05),
                            power_discovery_positive=exact_power(m, p, THRESH / 2),
                        )
                    )
            d = np.load(R / f"{n}_{s}_controls.npz")
            for (i, j), c in zip(d["low_pairs"], d["low"]):
                m = int(c.sum())
                st = for_size(c[None, :], m)[0]
                low.append(
                    dict(
                        n=n,
                        split=s,
                        test=CAT[i]["expression"],
                        target=TARGETS[j],
                        m=m,
                        **dict(zip(FIELDS, st)),
                        p=float(agreement_p(c)),
                    )
                )
            cs = d["shuffle"]
            pv = agreement_p(cs)
            m = int(cs[0].sum())
            first = for_size(cs[:1], m)[0]
            shuffle.append(
                dict(
                    n=n,
                    split=s,
                    m=m,
                    tests=len(cs),
                    reject=sum(pv <= 0.05 / (len(CAT) * len(SIZES))) if m else 0,
                    min_binom_p=pv.min() if m else np.nan,
                    sentinel="bit1(N)",
                    **dict(zip(FIELDS, first)),
                    sentinel_conditional_p=sentinel(d["shuffle_margins"]),
                )
            )
    order = np.argsort([c["p"] for c in low])
    adj = np.maximum.accumulate(
        np.array([low[i]["p"] for i in order]) * (len(order) - np.arange(len(order)))
    )
    for i, a in zip(order, adj):
        low[i]["p_holm"] = min(1.0, a)
    write_rows("planted.tsv", planted)
    write_rows("low_bit_controls.tsv", low)
    write_rows("shuffle.tsv", shuffle)


if __name__ == "__main__":
    control_tables()
