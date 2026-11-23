from pathlib import Path
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parent
SIZES = (24, 32, 40, 48, 56, 64, 80, 96)
SPLITS = {"D": 100000, "V": 50000, "C": 50000}


def seed(*parts):
    return int.from_bytes(
        hashlib.sha256(
            "|".join(map(str, ("ent-c-20261123-v1", *parts))).encode()
        ).digest()[:8],
        "big",
    )


def stem(n, s):
    return ROOT / "data" / f"{n}_{s}"


def rows(n, s):
    return np.memmap(
        str(stem(n, s)) + ".rows", dtype="<u8", mode="r", shape=(SPLITS[s], 4)
    )
