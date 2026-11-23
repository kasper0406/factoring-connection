from concurrent.futures import ThreadPoolExecutor
import subprocess, time, sys
from common import *


def compute(job):
    n, s = job
    path = str(stem(n, s))
    ts = time.monotonic()
    subprocess.run(
        [str(ROOT / "gates"), str(n), "data", path + ".rows", path + ".g", path + ".f"],
        check=True,
    )
    assert Path(path + ".g").stat().st_size == SPLITS[s] * 100
    assert Path(path + ".f").stat().st_size == SPLITS[s] * 2
    print(n, s, round(time.monotonic() - ts, 2), flush=True)


if __name__ == "__main__":
    with ThreadPoolExecutor(4) as pool:
        list(pool.map(compute, [(n, s) for s in SPLITS for n in SIZES]))
