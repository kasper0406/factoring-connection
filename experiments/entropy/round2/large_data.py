from concurrent.futures import ThreadPoolExecutor
from config import *


def chunk(job):
    n, i = job
    fn = DATA / f"C{n}_{i:03d}.pairs"
    if fn.is_file() and fn.stat().st_size == 16000000:
        return
    if n == 64:
        while not (DATA / "table64.ready").exists():
            time.sleep(1)
    tmp = Path(str(fn) + ".partial")
    t = time.monotonic()
    run(["./sampler", n, 2, 1000000, f"chunk{i}", tmp])
    assert tmp.stat().st_size == 16000000
    tmp.replace(fn)
    print(n, i, round(time.monotonic() - t, 2), flush=True)


if __name__ == "__main__":
    jobs = (
        [(96, i) for i in range(4)]
        + [(64, i) for i in range(40)]
        + [(96, i) for i in range(4, 220)]
        + [(80, i) for i in range(40)]
    )
    with ThreadPoolExecutor(6) as pool:
        list(pool.map(chunk, jobs))
