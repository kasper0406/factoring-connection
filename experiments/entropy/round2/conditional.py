"""Finite conditional distribution at the two largest sizes and a common-odds sensitivity interval."""

import heapq
from scipy.special import logsumexp, gammaln
from scipy.stats import hypergeom
from scipy.optimize import brentq
from scipy.signal import fftconvolve
from inference import *


def multiply(a, b):
    c = (
        np.convolve(a, b)
        if len(a) * len(b) < 120000
        else np.maximum(0, fftconvolve(a, b))
    )
    return c / c.sum()


def finite_test(mat, direction, key):
    support = mass_support(key)
    c = mat[support]
    cc = reduce_bins(c, c.sum(1) >= 8, COARSE[support]).astype("i8")
    m = cc.sum(1)
    x = cc[:, 2] + cc[:, 3]
    y = cc[:, 1] + cc[:, 3]
    tot = m.sum()
    lo = np.maximum(0, x + y - m)
    hi = np.minimum(x, y)
    active = hi > lo
    obs = int(cc[:, 3].sum() - lo.sum())
    rows, rep = np.unique(
        np.column_stack((m[active], x[active], y[active])), axis=0, return_counts=True
    )
    m, x, y = rows.T
    lo = np.maximum(0, x + y - m)
    width = np.minimum(x, y) - lo
    grid = np.arange(int(width.max()) + 1)
    a = grid + lo[:, None]
    lp = hypergeom.logpmf(a, m[:, None], x[:, None], y[:, None])
    lp = np.where(grid <= width[:, None], lp, -np.inf)
    lp -= logsumexp(lp, axis=1)[:, None]

    def gen(theta):
        w = lp + theta * grid
        z = logsumexp(w, axis=1)
        p = np.exp(w - z[:, None])
        return float(rep @ z), float(rep @ (p @ grid)), p

    _, e0, _ = gen(0)
    theta = brentq(lambda u: gen(u)[1] - obs, -35, 35)
    zstar, mean, p = gen(theta)
    heap = []
    idx = 0
    for row, w, n in zip(p, width, rep):
        a = row[: w + 1]
        if w == 1 and n > 1:
            b = binom.pmf(np.arange(n + 1), n, a[1])
        else:
            b = np.array([1.0])
            k = n
            while k:
                if k & 1:
                    b = multiply(b, a)
                k //= 2
                if k:
                    a = multiply(a, a)
        heapq.heappush(heap, (len(b), idx, b))
        idx += 1
    while len(heap) > 1:
        a = heapq.heappop(heap)[2]
        b = heapq.heappop(heap)[2]
        c = multiply(a, b)
        heapq.heappush(heap, (len(c), idx, c))
        idx += 1
    pmf = heap[0][2]

    def logtail(u, high):
        d = np.arange(len(pmf) - obs) if high else -np.arange(obs + 1)
        prob = pmf[obs:] if high else pmf[obs::-1]
        z, _, _ = gen(u)
        # For u on the appropriate side of the mean-matching tilt all weights are <=1.
        keep = prob > 0
        return (
            zstar
            - z
            + (u - theta) * obs
            + logsumexp(np.log(prob[keep]) + (u - theta) * d[keep])
        )

    consistent = direction * (obs - e0) >= 0
    if consistent:
        lnp = min(0.0, logtail(0, direction > 0))
    else:
        # Our locked upper/lower tail contains the mean. Use its complementary strict tail.
        d = np.arange(1, obs + 1) if direction > 0 else np.arange(1, len(pmf) - obs)
        p0 = pmf[obs - 1 :: -1] if direction > 0 else pmf[obs + 1 :]
        w = (theta if direction > 0 else -theta) * d
        mask = p0 > 0
        lt = zstar - theta * obs + logsumexp(np.log(p0[mask]) + w[mask])
        lnp = np.log1p(-np.exp(min(lt, -1e-15)))
    l = brentq(lambda u: logtail(u, True) - np.log(0.025), theta - 15, theta)
    h = brentq(lambda u: logtail(u, False) - np.log(0.025), theta, theta + 15)
    ci = sorted([direction * 2 * (gen(z)[1] - e0) / tot for z in (l, h)])
    return dict(
        m=tot,
        informative_strata=int(rep.sum()),
        conditional_log_p=lnp,
        conditional_p=np.exp(lnp),
        signed_centered=direction * 2 * (obs - e0) / tot,
        common_log_odds=direction * theta,
        common_odds_conditional_adv_lo=ci[0],
        common_odds_conditional_adv_hi=ci[1],
    )


def main():
    for cohort in ("original", "augmented"):
        out = []
        for n in (80, 96):
            while not (RES / f"{cohort}{n}_measure.tsv").exists():
                time.sleep(3)
            km, ar = plans(n)
            cache = {}
            for k, mat, _ in tables(n, cohort):
                if k.startswith("t"):
                    sign = next(x[1] for x in ar if x[0] == k)
                    for v in (0, 1):
                        cache[k, v] = finite_test(mat[v], sign, k)
            for i, (k, sign, target) in enumerate(ar):
                if (k, 0) in cache:
                    for v in (0, 1):
                        out.append(
                            dict(
                                cohort=cohort,
                                n=n,
                                rank=i + 1,
                                key=k,
                                target=target,
                                view="clean" if v else "all",
                                **cache[k, v],
                            )
                        )
        write_records(f"{cohort}_conditional.tsv", out)
        print(cohort, "complete", flush=True)


if __name__ == "__main__":
    main()
