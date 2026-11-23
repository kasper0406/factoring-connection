import subprocess, time, sys
from concurrent.futures import ThreadPoolExecutor
from common import *


def make(job):
    n, s = job
    out = str(stem(n, s)) + ".rows"
    t = time.monotonic()
    subprocess.run(
        [
            str(ROOT / "sample"),
            str(n),
            str(SPLITS[s]),
            str(seed("primary", n, s)),
            str(seed("reference", n, s)),
            out,
        ],
        check=True,
    )
    print(n, s, "sampled", round(time.monotonic() - t, 1), flush=True)


if __name__ == "__main__":
    with ThreadPoolExecutor(4) as pool:
        list(pool.map(make, [(n, s) for s in SPLITS for n in SIZES]))
