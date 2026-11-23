from scipy.stats import chi2
from controls import *


def fit(n, adv, m, poly):
    grid = np.arange(801) / 100 if poly else np.arange(501) / 1000
    A = (n / 24.0) ** (-grid[:, None]) if poly else np.exp2(-grid[:, None] * (n - 24))
    w = 4 * m
    amp = np.maximum(0, ((A * w) * adv).sum(1) / (A * A * w).sum(1))
    rss = ((A * amp[:, None] - adv) ** 2 * w).sum(1)
    i = int(rss.argmin())
    return amp[i], grid[i], rss[i]


def main():
    surv = list(csv.DictReader(open(R / "survivors.tsv"), delimiter="\t"))
    records = []
    bysize = {}
    cutoff = 0.05 / (4 * max(1, len(surv)))
    for n in SIZES:
        for s in ("V", "C"):
            x = np.load(R / f"{n}_{s}_counts.npz")
            c = x["counts"]
            m = c.sum(-1)
            bysize[(n, s)] = (c, m)
    for i, row in enumerate(surv):
        d = int(row["direction"])
        tid = int(row["test"])
        j = int(row["target_id"])
        ns = []
        ad = []
        ms = []
        passes = []
        endpoints = {}
        for n in SIZES:
            for s in ("V", "C"):
                c, m = bysize[(n, s)]
                count = int(m[i, 3])
                a = c[i, 3]
                hit = int(a[0] + a[3])
                hit = hit if d == 1 else count - hit
                av = hit / count - 0.5 if count else np.nan
                if s == "C" and count >= 100:
                    ns.append(n)
                    ad.append(av)
                    ms.append(count)
                if n >= 80:
                    p = binom.sf(hit - 1, count, 0.5) if count else 1
                    passes.append(av > 0 and p <= cutoff)
                    endpoints[f"adv_{n}_{s}"] = av
                    endpoints[f"p_{n}_{s}"] = p
        ns, ad, ms = map(np.asarray, (ns, ad, ms))
        ap, k, rp = fit(ns, ad, ms, True)
        ae, b, re = fit(ns, ad, ms, False)
        fits = bool(len(ns) >= 5 and rp <= chi2.ppf(0.95, len(ns) - 2) and re - rp >= 2)
        records.append(
            dict(
                test=tid,
                target=TARGETS[j],
                direction=d,
                n_fit=len(ns),
                amp_poly_at24=ap,
                k=k,
                rss_poly=rp,
                amp_exp_at24=ae,
                beta=b,
                rss_exp=re,
                delta=re - rp,
                gof_cutoff=chi2.ppf(0.95, len(ns) - 2),
                scale_pass=int(fits),
                endpoints_pass=int(all(passes)),
                lead=int(all(passes) and fits),
                **endpoints,
            )
        )
    write_rows("scaling.tsv", records)
    print(
        "survivors",
        len(surv),
        "cutoff",
        cutoff,
        "scale passes",
        sum(x["scale_pass"] for x in records),
        "endpoint passes",
        sum(x["endpoints_pass"] for x in records),
        "leads",
        sum(x["lead"] for x in records),
    )


if __name__ == "__main__":
    main()
