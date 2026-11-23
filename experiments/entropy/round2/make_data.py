from concurrent.futures import ThreadPoolExecutor
from config import *
import sys


def forbid(n):
    blocks = []
    for phase in "DVC":
        a = np.fromfile(R.parent / "data" / f"{n}_{phase}.rows", dtype="<u8").reshape(
            -1, 2
        )
        blocks.append(a)
    return np.unique(np.vstack(blocks).view("V16").ravel())


def small(n):
    run(["./sampler", "table", n])
    prev = forbid(n)
    for si, (s, m) in enumerate([("D", 240000), ("V", 120000), ("C", 250000)]):
        fn = path(n, s)
        tmp = Path(str(fn) + ".tmp")
        run(["./sampler", n, si, int(m * 1.12) + 2000, "low0", tmp])
        arr = np.fromfile(tmp, "<u8").reshape(-1, 2)
        v = arr.view("V16").ravel()
        eligible = ~np.isin(v, prev)
        if s != "C":
            _, idx = np.unique(v, return_index=True)
            keep = np.zeros(len(v), bool)
            keep[idx] = True
            eligible &= keep
        xn = arr[eligible][:m]
        assert len(xn) == m
        xn.tofile(fn)
        tmp.unlink()
        print(n, s, len(xn), flush=True)


if __name__ == "__main__":
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(small, [40, 48, 56]))
    run(["./sampler", "table", 64])
    (DATA / "table64.ready").touch()
    print("64 prime tables ready", flush=True)
