from config import *
from concurrent.futures import ThreadPoolExecutor


def prepare(n):
    fn = DATA / f"round1_{n}.u64"
    if not fn.exists():
        with open(fn, "wb") as out:
            for s in "DVC":
                arr = np.memmap(
                    R.parent / "data" / f"{n}_{s}.rows", dtype="<u8", mode="r"
                ).reshape(-1, 2)
                out.write(arr)
    return fn


def go(job):
    n, f = job
    run(
        [
            "./budget",
            n,
            DATA / f"round1_{n}.u64",
            f,
            str(f) + ".budget",
            str(f) + ".hash",
        ]
    )


if __name__ == "__main__":
    jobs = []
    for n in (40, 48, 56, 64, 80, 96):
        prepare(n)
        if n <= 56:
            jobs.extend((n, path(n, s)) for s in "DVC")
        else:
            jobs.extend((n, f) for f in sorted(DATA.glob(f"C{n}_*.pairs")))
    assert len(jobs) == 309
    with ThreadPoolExecutor(5) as pool:
        list(pool.map(go, jobs))
    print(len(jobs), "budget files")
