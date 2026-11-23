"""Finite-sample intervals and estimates; field order is saved inside every archive."""

from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from scipy.special import betaincinv, entr
from scoring import *

FIELDS = ["adv", "adv_lo", "adv_hi", "mi_bits", "mi_lo", "mi_hi"]
R = ROOT / "results"


def cp(k, m, alpha=0.05):
    k = np.asarray(k)
    a = alpha / 2
    lo = np.zeros(k.shape)
    hi = np.ones(k.shape)
    ok = k > 0
    lo[ok] = betaincinv(k[ok], m - k[ok] + 1, a)
    ok = k < m
    hi[ok] = 1 - betaincinv(m - k[ok], k[ok] + 1, a)
    return lo, hi


def hpart(x):
    return entr(x) / np.log(2)


def hbinary(x):
    return hpart(x) + hpart(1 - x)


def entropy_box(lo, hi, binary):
    h = hbinary if binary else hpart
    left = h(lo)
    right = h(hi)
    top_at = 0.5 if binary else np.exp(-1)
    return np.minimum(left, right), h(np.clip(top_at, lo, hi))


def for_size(c, m):
    if not m:
        return np.tile([np.nan, -0.5, 0.5, np.nan, 0, 1], (len(c), 1))
    c = c.astype("i4")
    hit = c[:, 0] + c[:, 3]
    u, inv = np.unique(hit, return_inverse=True)
    a, b = cp(u, m)
    # Invert each distinct integer once. Six simultaneous intervals bound H(X)+H(Y)-H(X,Y).
    k = np.column_stack((c[:, 2] + c[:, 3], c[:, 1] + c[:, 3], c))
    unique, idx = np.unique(k, return_inverse=True)
    idx = idx.reshape(k.shape)
    l, h = cp(unique, m, 0.05 / 6)
    ml, mh = entropy_box(l[idx[:, :2]], h[idx[:, :2]], True)
    jl, jh = entropy_box(l[idx[:, 2:]], h[idx[:, 2:]], False)
    lower = np.maximum(0, ml.sum(1) - jh.sum(1))
    upper = np.minimum.reduce([np.ones(len(c)), mh.min(1), mh.sum(1) - jl.sum(1)])
    est = np.maximum(0, hbinary(k[:, :2] / m).sum(1) - hpart(c / m).sum(1))
    return np.column_stack(
        (hit / m - 0.5, a[inv] - 0.5, b[inv] - 0.5, est, lower, np.maximum(0, upper))
    )


def summarize_archive(job):
    n, s = job
    data = np.load(R / f"{n}_{s}_counts.npz")
    c = data["counts"]
    out = np.empty(c.shape[:-1] + (6,), dtype="f8")
    for v in range(c.shape[1]):
        ms = c[:, v].sum(-1)
        for m in np.unique(ms):
            ids = np.flatnonzero(ms == m)
            out[ids, v] = for_size(c[ids, v], int(m))
    good = np.isfinite(out[..., 0])
    assert np.all(out[..., 3][good] + 1e-12 >= out[..., 4][good])
    assert np.all(out[..., 3][good] - 1e-12 <= out[..., 5][good])
    np.savez_compressed(R / f"{n}_{s}_stats.npz", fields=FIELDS, views=VIEWS, stats=out)
    print(n, s, "statistics", flush=True)


if __name__ == "__main__":
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(summarize_archive, [(n, s) for s in SPLITS for n in SIZES]))
