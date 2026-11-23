from inference import *


def main():
    pos = [hashlib.sha256() for _ in CAT]
    neg = [hashlib.sha256() for _ in CAT]
    rare = []
    for n in (40, 48, 56):
        d = Dataset(n, "D")
        u = np.array([pack(d.x(t)) for t in range(NU)])
        assert d.m % 64 == 0
        for t, (op, a, b, c) in enumerate(OPS):
            v = (
                u[t]
                if not op
                else (
                    u[a] ^ u[b]
                    if op == 1
                    else (u[a] & u[b]) | (u[b] & u[c]) | (u[a] & u[c])
                )
            )
            pos[t].update(v)
            neg[t].update(~v)
        z = np.load(RES / f"D{n}_screen.npz")
        m = int(z["matched_pairs"]) * 2
        r = z["ones"] / m
        freq = np.minimum(r, 1 - r)
        rare.append(
            dict(
                n=n,
                paired_rows=m,
                degenerate=int((freq == 0).sum()),
                below_cutoff=int((freq < 1e-4).sum()),
                rare=int(((freq >= 1e-4) & (freq < 1e-2)).sum()),
                common=int((freq >= 1e-2).sum()),
            )
        )
    groups = {}
    for i, (a, b) in enumerate(zip(pos, neg)):
        key = min(a.digest(), b.digest())
        groups.setdefault(key, []).append(i)
    write_records("discovery_frequency.tsv", rare)
    write_records(
        "empirical_catalog.tsv",
        [
            dict(test=t, first_on_D=v[0], class_size=len(v))
            for v in groups.values()
            for t in v
        ],
    )
    print(
        len(groups), "distinct up to complements on 720000 discovery rows", flush=True
    )


if __name__ == "__main__":
    main()
