from confirm_eval import *
from labels import ybits
import math

H = 1024 * 1024
RESIDUES = np.array([i for i in range(3840) if math.gcd(i, 3840) == 1])
assert len(RESIDUES) == 1024
COARSE = np.repeat(RESIDUES % 256 // 2, 1024) * 64 + np.tile(
    np.arange(1024) // 16, 1024
)


def plans(n):
    unique = {}
    allrows = []
    for rank, row in enumerate(locks(n)):
        j = int(row["target"])
        t = int(row["test"])
        j0 = j % len(BASE)
        if BASE[j0].startswith("gen"):
            j = j0
        key = f"t{t}_y{j}"
        if key not in unique:
            unique[key] = (rank, t, j)
        allrows.append((key, int(row["direction"]), TARGETS[int(row["target"])]))
    for name in ("plant1", "plant16sq", "plant1sq", "negative"):
        unique[name] = (-1, -1, -1)
        allrows.append((name, 1, name))
    return unique, allrows


def subset_target(p, q, coin, n, j):
    name = BASE[j % len(BASE)]
    f = np.where(coin, q, p) if j >= len(BASE) else p
    if name.startswith("gen"):
        d = int(name[3:])

        def ch(x):
            if d == -1:
                return x % 4 == 3
            if d == 2:
                return np.isin(x % 8, [3, 5])
            ell = abs(d)
            lut = np.array([pow(v, (ell - 1) // 2, ell) == ell - 1 for v in range(ell)])
            return lut[x % ell]

        elig = (p % 8 * (q % 8) % 8 == 1) if d == 2 else ch(p) == ch(q)
        return ch(f), elig
    if name == "mod3half":
        return f % 3 == 1, np.ones(len(f), bool)
    # The general path includes 192-bit PIT and all original readouts.
    y, ok, _ = ybits(np.column_stack((p, q)), n, coin, only=[j])
    return y[j], ok[j]


def coded(n, fn, kmap):
    a = np.memmap(fn, dtype="<u8", mode="r").reshape(-1, 2)
    p, q = a.T
    h, xb, fb, ctrl = (
        np.memmap(str(fn) + ".locked", dtype="<u4", mode="r").reshape(-1, 4).T
    )
    b = np.memmap(str(fn) + ".budget", dtype="u1", mode="r")
    coin = (ctrl & 1) != 0
    mid = (ctrl >> 1 & 1).astype(bool)
    target = {}

    def gen():
        for name, (rank, t, j) in kmap.items():
            if rank >= 0:
                if j not in target:
                    target[j] = subset_target(p, q, coin, n, j)
                y, eligible = target[j]
                x = (xb >> rank) & 1
                fl = (fb >> rank) & 1
            elif name == "negative":
                x = (ctrl >> 5) & 1
                y = (ctrl >> 6) & 1
                eligible = (ctrl >> 7 & 1) != 0
                fl = 0
            else:
                bit = {"plant1": 2, "plant16sq": 3, "plant1sq": 4}[name]
                x = mid ^ ((ctrl >> bit & 1) == 0)
                y = mid
                eligible = np.ones(len(a), bool)
                fl = 0
            good = eligible & ((b & 48) == 0)
            clean = good & ((b & 15) == 0) & (fl == 0)
            yield name, x, y, good, clean, fl

    return a, h, b, gen()


def aggregate(n, kind="C"):
    km, ar = plans(n)
    if kind == "S" and n == 96:
        km = {k: v for k, v in km.items() if v[0] < 4}
    files = (
        jobs(n)
        if kind == "C"
        else [DATA / f"S{n}_{i:03d}.pairs" for i in range(128 if n == 64 else 170)]
    )
    hist = np.lib.format.open_memmap(
        RES / f"{kind}{n}_strata.npy", "w+", dtype="u4", shape=(len(km), 2, H, 4)
    )
    hist[:] = 0
    reasons = np.zeros((len(km), 32), "i8")
    budget = np.zeros(64, "i8")
    for fn in files:
        while not Path(str(fn) + ".locked").exists():
            time.sleep(3)
        a, h, b, columns = coded(n, fn, km)
        budget += np.bincount(b, minlength=64)
        h4 = h * 4
        for i, (k, x, y, ok, clean, fl) in enumerate(columns):
            ix = h4 + 2 * x.astype("u4") + y
            for v, m in enumerate((ok, clean)):
                if m.any():
                    c = np.bincount(ix[m], minlength=H * 4)
                    hist[i, v].reshape(-1)[:] += c.astype("u4")
            reasons[i] += np.bincount(((b & 15) + 16 * fl)[ok], minlength=32)
    hist.flush()
    np.savez(
        RES / f"{kind}{n}_masks.npz",
        keys=list(km),
        by_measurement=reasons,
        budget=budget,
    )
    print(n, budget.sum(), list(km), flush=True)


if __name__ == "__main__":
    import time

    for n in (40, 48, 56, 64, 96, 80):
        while any(not Path(str(f) + ".locked").exists() for f in jobs(n)):
            time.sleep(2)
        aggregate(n)
