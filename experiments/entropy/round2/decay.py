from inference import *
from scipy.stats import norm, chi2

MODELS = {
    "poly": (0.0, 6.0),
    "barrier": (1 / 6, 2.0),
    "free": (0.02, 2.0),
    "sub": (0.0, 12.0),
    "wide_poly": (0.0, 32.0),
}


def design(ns, kind):
    if kind in ("poly", "wide_poly"):
        return np.log(ns / ns[0])
    if kind == "sub":
        return np.sqrt(ns * np.log(ns)) - np.sqrt(ns[0] * np.log(ns[0]))
    return (ns - ns[0]) * np.log(2)


def profile(Y, ns, se, kind):
    # Row-wise, profiled nonnegative amplitude with scalar slope refinement at the global grid optimum.
    Y = np.atleast_2d(Y)
    t = design(ns, kind)
    l, u = MODELS[kind]
    weight = 1 / se**2
    grid = np.linspace(l, u, 257)
    F = np.exp(-grid[:, None] * t)
    dot = Y @ (F * weight).T
    den = (F * F) @ weight
    ix = (np.maximum(dot, 0) ** 2 / den).argmax(1)
    a = grid[np.maximum(ix - 1, 0)]
    b = grid[np.minimum(ix + 1, len(grid) - 1)]
    yw = Y * weight

    def value(s):
        f = np.exp(-s[:, None] * t)
        amp = np.maximum(0, (f * yw).sum(1) / (f * f * weight).sum(1))
        return ((Y - amp[:, None] * f) ** 2 * weight).sum(1), amp

    for _ in range(32):
        v = (2 * a + b) / 3
        w = (a + 2 * b) / 3
        left = value(v)[0] <= value(w)[0]
        a = np.where(left, a, v)
        b = np.where(left, w, b)
    slope = (a + b) / 2
    rss, amp = value(slope)
    for k in (l, u):
        v = np.full(len(Y), k)
        r, A = value(v)
        ok = r < rss
        slope = np.where(ok, v, slope)
        amp = np.where(ok, A, amp)
        rss = np.minimum(rss, r)
    return rss, amp, slope


def local(ns, y, se):
    if len(ns) < 2:
        return np.nan, np.inf
    _, amp, k = profile(y, ns, se, "wide_poly")
    w = 1 / se**2
    t = design(ns, "poly")
    f = np.exp(-k[0] * t)
    J = np.column_stack((f, -amp[0] * f * t))
    info = (J.T * w) @ J
    return k[0], np.sqrt(np.linalg.inv(info)[1, 1]) if np.linalg.det(
        info
    ) > 1e-14 else np.inf


def rowsfile(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def process(cohort):
    packed = {}
    for n in (40, 48, 56, 64, 80, 96):
        dy = {}
        if n <= 48:
            dy = {
                d["key"]: float(d["dyad_se"])
                for d in rowsfile(RES / f"original{n}_dyad.tsv")
            }
        for r in rowsfile(RES / f"{cohort}{n}_measure.tsv"):
            ident = (r["rank"], r["target"])
            packed.setdefault(ident, {}).setdefault(n, {})[r["view"]] = r
            if n <= 48 and r["key"] in dy:
                r["dyad_se"] = dy[r["key"]]
    results = []
    cache = {}
    for ident, data in packed.items():
        rank, target = ident
        organic = int(rank) > 0
        key = next(iter(data.values()))["clean"]["key"]
        sig = int(next(iter(data.values()))["clean"]["sign"])
        cachekey = (key, sig)
        if cachekey in cache:
            results.append(dict(cache[cachekey], rank=rank, target=target))
            continue
        data = {n: v for n, v in data.items() if "clean" in v and "all" in v}
        ns = np.array(sorted(data))
        y = []
        se = []
        z = []
        changes = []
        thresholds = []
        for n in ns:
            rc = data[n]["clean"]
            ra = data[n]["all"]
            if organic:
                effect = float(rc["signed_centered"])
                inc = float(ra["signed_centered"])
                s = float(rc["centered_se"])
                if n <= 48:
                    s = max(s, float(rc["dyad_se"]), float(rc["mh_se"]))
            else:
                effect = float(rc["signed_agreement"])
                inc = float(ra["signed_agreement"])
                s = np.sqrt((0.25 - effect**2) / int(rc["m"]))
            y.append(effect)
            se.append(s)
            z.append(float(rc["z"]))
            thresholds.append(
                float(rc["holm_p"]) if organic else float(rc["one_sided_p"])
            )
            if n >= 64:
                changes.append(abs(inc - effect) / max(abs(effect), 1e-30))
        y = np.array(y)
        se = np.array(se)
        ok = se > 0
        ns = ns[ok]
        y = y[ok]
        se = se[ok]
        fit = {k: tuple(float(a[0]) for a in profile(y, ns, se, k)) for k in MODELS}
        expP = {}
        rng = np.random.default_rng(seed("bootstrap", cohort, key))
        for k in ("barrier", "free", "sub"):
            rr, am, rate = fit[k]
            mu = am * np.exp(-rate * design(ns, k))
            sim = mu + rng.normal(size=(2000, len(ns))) * se
            d = profile(sim, ns, se, k)[0] - profile(sim, ns, se, "poly")[0]
            obs = rr - fit["poly"][0]
            expP[k] = (int((d >= obs - 1e-10).sum()) + 1) / 2001
        rp, A, k = fit["poly"]
        rb, _, b = fit["barrier"]
        rf, _, f = fit["free"]
        rs, _, g = fit["sub"]
        low = ns <= 56
        high = ns >= 64
        kl, sl = local(ns[low], y[low], se[low])
        kh, sh = local(ns[high], y[high], se[high])
        c1 = all(
            n in data
            and float(data[n]["clean"]["z"]) >= 3
            and (not organic or float(data[n]["clean"]["holm_p"]) <= 0.01)
            for n in [80, 96]
        )
        c2 = rb - rp >= 9.21 and expP["barrier"] <= 0.01
        c3 = (
            rp <= chi2.ppf(0.99, len(ns) - 2)
            and np.isfinite(sh)
            and kh - kl <= 2 * np.hypot(sh, sl)
        )
        c4 = max(changes) < 0.2 if changes else False
        o = dict(
            cohort=cohort,
            rank=rank,
            key=key,
            target=target,
            n_sizes=len(ns),
            sizes=",".join(map(str, ns)),
            poly_rss=rp,
            poly_amplitude_at_first=A,
            k=k,
            exp_barrier_rss=rb,
            beta_barrier=b,
            LR_barrier_vs_poly=rb - rp,
            bootstrap_p_barrier=expP["barrier"],
            exp_free_rss=rf,
            beta_free=f,
            LR_free_vs_poly=rf - rp,
            bootstrap_p_free=expP["free"],
            sub_rss=rs,
            gamma=g,
            LR_sub_vs_poly=rs - rp,
            bootstrap_p_sub=expP["sub"],
            wide_poly_rss=fit["wide_poly"][0],
            wide_k=fit["wide_poly"][2],
            k_low=kl,
            se_k_low=sl,
            k_high=kh,
            se_k_high=sh,
            max_relative_split_shift=max(changes) if changes else np.nan,
            C1=int(c1),
            C2=int(c2),
            C3=int(c3),
            C4=int(c4),
            lead=int(organic and c1 and c2 and c3 and c4),
        )
        cache[cachekey] = o
        results.append(o)
    write_records(f"{cohort}_decay.tsv", results)
    print(cohort, len(results), "leads", sum(r["lead"] for r in results), flush=True)


if __name__ == "__main__":
    for cohort in ("original", "augmented"):
        while (
            not all(
                (RES / f"{cohort}{n}_measure.tsv").exists()
                for n in (40, 48, 56, 64, 80, 96)
            )
            or not (RES / "original48_dyad.tsv").exists()
        ):
            time.sleep(3)
        process(cohort)
