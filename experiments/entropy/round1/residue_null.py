"""Post hoc low-bit calibration with an exactly uniform, product-compatible residue."""

import csv
from controls import write_rows
from statistics_run import *


def main():
    report = []
    inverse = np.array([pow(u, -1, 16) if u % 2 else 0 for u in range(16)], dtype="u1")
    for n in SIZES:
        for s, m in SPLITS.items():
            rng = np.random.default_rng(seed("posthoc-residue-null", n, s))
            a = rows(n, s)
            residues = ((a[:, [0, 2]] % 16) * (a[:, [1, 3]] % 16) % 16).astype("u1")
            # For each product independently, U is uniform odd and U*V = N modulo 16.
            u = 2 * rng.integers(0, 8, size=residues.shape, dtype="u1") + 1
            v = residues * inverse[u] % 16
            assert np.array_equal(u * v % 16, residues)
            flags = np.fromfile(str(stem(n, s)) + ".f", "u1").reshape(-1, 2)
            for b in (1, 2, 3):
                x = (residues[:, 0] >> b) & 1
                for k in (1, 2, 3):
                    y = (u[:, 0] >> k) & 1
                    yr = (u[:, 1] >> k) & 1
                    for clean in (False, True):
                        keep = y != yr
                        if clean:
                            keep &= (flags == 0).all(1)
                        c = np.bincount((2 * x + y)[keep], minlength=4)
                        count = int(c.sum())
                        stats = for_size(c[None, :], count)[0]
                        report.append(
                            dict(
                                n=n,
                                split=s,
                                view="matched_clean" if clean else "matched",
                                product_bit=b,
                                synthetic_residue_bit=k,
                                m=count,
                                **dict(zip(FIELDS, stats)),
                                p=float(agreement_p(c)),
                            )
                        )
    order = np.argsort([d["p"] for d in report])
    size = len(order)
    p_adjusted = np.maximum.accumulate(
        np.array([report[i]["p"] for i in order]) * (size - np.arange(size))
    )
    for i, p in zip(order, p_adjusted):
        report[i]["p_holm"] = min(1.0, p)
    write_rows("uniform_residue_null.tsv", report)
    print(
        "comparisons", size, "Holm rejections", sum(d["p_holm"] <= 0.05 for d in report)
    )


if __name__ == "__main__":
    main()
