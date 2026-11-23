from config import *
from functools import lru_cache
import gmpy2 as g
import json, sys

BASE = (
    [f"pit{k}" for k in range(1, 5)]
    + [f"gl{k}" for k in range(8)]
    + [f"pitgl{k}" for k in range(4)]
    + ["midhalf", "midhalf1", "midthird"]
    + [f"low{k}" for k in range(1, 5)]
    + [f"gen{d}" for d in (-1, 2, -3, 5, -7, -11)]
    + ["cubic7", "quartic5", "mod3half", "mod5half"]
)
TARGETS = [f"{fac}/{name}" for fac in ("min", "random") for name in BASE]


def masklist(n):
    h = n // 2
    top = max(4, (n - 1).bit_length())
    low = (1 << (h // 2 + 1)) - 2
    out = []
    i = 0
    while len(out) < 8:
        m = int.from_bytes(
            hashlib.shake_128(f"GL|{n}|{i}".encode()).digest(8), "little"
        ) & ((1 << (h - top)) - 2)
        i += 1
        if (m & low).bit_count() >= 4 and m not in out:
            out.append(m)
    return out


PITMASK = [0x9B75, 0xE71A, 0xADB3, 0x766D]


def pit_words(a, n, coin):
    with g.context(g.get_context(), precision=192):
        h = n // 2
        lo = 3 << (h - 2)
        hi = 1 << h
        ll = g.log(lo)
        lh = g.log(hi)
        one = g.mpfr(1)
        log = lru_cache(maxsize=65536)(g.log)
        small = np.empty(len(a), "<u2")
        rand = small.copy()
        ok = np.empty(len(a), bool)
        for i, (u, v) in enumerate(a):
            p, q = int(u), int(v)
            N = p * q
            lp, lq = log(p), log(q)
            L = lp + lq
            bound_is_lo = N <= lo * hi
            sm = ll if bound_is_lo else L - lh
            lower = lo if bound_is_lo else g.mpfr(N) / hi
            ok[i] = g.sqrt(N) - lower >= 16 * L
            U = one - g.log(lq / lp) / g.log((L - sm) / sm)
            assert 0 <= U <= 1
            word = min(65535, int(g.floor(65536 * U)))
            small[i] = word
            star = one - U / 2 if coin[i] else U / 2
            rand[i] = min(65535, int(g.floor(65536 * star)))
        return small, rand, ok


def ybits(a, n, coin, only=None):
    # only controls the optional high-scale PIT work on externally locked organic targets.
    p, q = a.T
    h = n // 2
    uses_pit = only is None or any(BASE[j % len(BASE)].startswith("pit") for j in only)
    if uses_pit:
        w1, w2, pok = pit_words(a, n, coin)
    else:
        w1 = w2 = np.zeros(len(a), dtype="<u2")
        pok = np.ones(len(a), bool)
    out = []
    elig = []
    masks = masklist(n)
    allok = np.ones(len(a), bool)
    for f, w in [(p, w1), (np.where(coin, q, p), w2)]:
        ys = [(w >> (16 - k)) & 1 for k in (1, 2, 3, 4)]
        es = [pok] * 4
        ys.extend(np.bitwise_count(f & m) & 1 for m in masks)
        es.extend([allok] * 8)
        ys.extend(np.bitwise_count(w & m) & 1 for m in PITMASK)
        es.extend([pok] * 4)
        ys.extend((f >> k) & 1 for k in (h // 2, h // 2 + 1, h // 3, 1, 2, 3, 4))
        es.extend([allok] * 7)
        for d in (-1, 2, -3, 5, -7, -11):

            def chi(z):
                if d == -1:
                    return (z % 4) == 3
                if d == 2:
                    return np.isin(z % 8, [3, 5])
                l = abs(d)
                return np.array([pow(j, (l - 1) // 2, l) == l - 1 for j in range(l)])[
                    z % l
                ]

            ys.append(chi(f))
            es.append(((p % 8) * (q % 8) % 8 == 1) if d == 2 else (chi(p) == chi(q)))
        ys.extend(
            [np.isin(f % 7, [1, 6]), f % 5 == 1, f % 3 == 1, np.isin(f % 5, [1, 2])]
        )
        es.extend(
            [
                p % 7 * (q % 7) % 7 == 1,
                (p % 5 * (q % 5) % 5 == 1)
                & np.isin(p % 5, [1, 4])
                & np.isin(q % 5, [1, 4]),
                allok,
                allok,
            ]
        )
        out.extend(ys)
        elig.extend(es)
    return np.array(out, bool), np.array(elig, bool), w1


def prepare(n, s):
    a = load(n, s)
    rng = np.random.default_rng(seed("orientation", n, s))
    coin = rng.integers(0, 2, len(a), dtype="u1").astype(bool)
    y, ok, w = ybits(a, n, coin)
    fn = DATA / f"{n}_{s}.labels.npz"
    np.savez_compressed(
        fn,
        y=np.packbits(y, axis=1),
        ok=np.packbits(ok, axis=1),
        coin=np.packbits(coin),
        pit=w,
    )
    print(n, s, len(a), flush=True)


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(4) as ex:
        jobs = [ex.submit(prepare, n, s) for s in "DV" for n in (40, 48, 56)]
        for j in jobs:
            j.result()
    (RES / "targets.json").write_text(
        json.dumps(
            {
                "names": TARGETS,
                "GL": {n: masklist(n) for n in [40, 48, 56, 64, 80, 96]},
                "PIT_GL": PITMASK,
            },
            indent=2,
        )
    )
