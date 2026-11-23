from pathlib import Path
import hashlib, os, subprocess, time
import numpy as np

R = Path(__file__).resolve().parent
DATA = R / "data"
RES = R / "results"
FLINT = "/Volumes/git/factoring-connection/.venv/lib/python3.12/site-packages/flint/.dylibs/libflint.24.0.dylib"
ENV = {**os.environ, "FLINTLIB": FLINT}


def seed(*key):
    return int.from_bytes(
        hashlib.sha256(("|".join(map(str, ("ent-c-R2-v1", *key)))).encode()).digest()[
            :8
        ],
        "little",
    )


def run(args, **kw):
    return subprocess.run(list(map(str, args)), cwd=R, env=ENV, check=True, **kw)


def path(n, s):
    return DATA / f"{n}_{s}.u64"


def load(n, s):
    return np.memmap(path(n, s), "<u8", mode="r").reshape(-1, 2)
