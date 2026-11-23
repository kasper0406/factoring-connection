from grammar import *
import sys


def low(job):
    n, s = job
    vs = json.loads((RES / "catalog.json").read_text())["eligible"]
    src = path(n, s)

    # Concurrent independent variant blocks; each output is variant-major.
    def block(tup):
        start, group = tup
        dst = DATA / f"{n}_{s}_v{start}.g"
        engine("features", group, src, dst)
        assert dst.stat().st_size == len(group) * src.stat().st_size // 16 * 10

    with ThreadPoolExecutor(5) as pool:
        list(pool.map(block, [(j, vs[j : j + 3]) for j in range(0, len(vs), 3)]))
    print(n, s, "features", flush=True)


if __name__ == "__main__":
    for s in sys.argv[1] if len(sys.argv) > 1 else "D":
        for n in [40, 48, 56]:
            low((n, s))
