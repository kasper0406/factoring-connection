"""Conditional MI above its finite-sample null expectation in the locked CF tests."""

from decay import rowsfile
from inference import COARSE, RES, reduce_bins, tables, write_records
import numpy as np
from scipy.special import gammaln, xlogy


def excess(cells):
    cells = np.asarray(cells, dtype="i8")
    n = cells.sum(axis=1)
    m = n.sum()
    x = np.minimum(cells[:, 2:].sum(axis=1), cells[:, :2].sum(axis=1))
    y = np.minimum(cells[:, 1::2].sum(axis=1), cells[:, ::2].sum(axis=1))
    use = (x > 0) & (y > 0)
    observed = xlogy(cells[use], cells[use]).sum()
    n, x, y = n[use], x[use], y[use]
    if not n.size:
        return 0.0
    # Swapping levels preserves MI. Both margins are now at most n/2.
    total = np.zeros_like(n, dtype=float)
    norm = np.zeros_like(total)
    constant = (
        gammaln(y + 1)
        + gammaln(n - y + 1)
        + gammaln(x + 1)
        + gammaln(n - x + 1)
        - gammaln(n + 1)
    )
    for k in range(int(np.minimum(x, y).max()) + 1):
        active = (x >= k) & (y >= k)
        remaining = [
            np.full(active.sum(), k),
            x[active] - k,
            y[active] - k,
            n[active] - x[active] - y[active] + k,
        ]
        probability = np.exp(constant[active] - sum(gammaln(z + 1) for z in remaining))
        total[active] += probability * sum(xlogy(z, z) for z in remaining)
        norm[active] += probability
    assert np.allclose(norm, 1, rtol=0, atol=1e-8)
    return (observed - (total / norm).sum()) / (m * np.log(2))


def main():
    # The two 2x2 tables below occur with equal probability under their fixed margins.
    assert abs(excess([[0, 1, 1, 0]]) + excess([[1, 0, 0, 1]])) < 1e-14
    assert np.isclose(excess([[2, 0, 0, 2]]), 2 / 3)
    assert np.isclose(excess([[1, 1, 1, 1]]), -1 / 3)
    outputs = []
    for cohort in ("original", "augmented"):
        for n in (40, 48, 56, 64, 80, 96):
            if cohort == "augmented" and n not in (64, 96):
                outputs.extend(
                    dict(r, cohort=cohort)
                    for r in list(outputs)
                    if r["cohort"] == "original" and r["n"] == n
                )
                continue
            measured = rowsfile(RES / f"{cohort}{n}_measure.tsv")
            for key, views, _ in tables(n, cohort):
                locks = [r for r in measured if r["key"] == key and int(r["rank"])]
                if not locks:
                    continue
                for v, view in enumerate(("all", "clean")):
                    raw = views[v]
                    counts = reduce_bins(raw, raw.sum(1) >= 8, COARSE)
                    adj = excess(counts)
                    for r in locks:
                        if r["view"] == view:
                            observed = float(r["conditional_mi_bits"])
                            outputs.append(
                                dict(
                                    cohort=cohort,
                                    n=n,
                                    rank=r["rank"],
                                    key=key,
                                    target=r["target"],
                                    view=view,
                                    conditional_plugin_bits=observed,
                                    fixed_margin_null_mean_bits=observed - adj,
                                    conditional_excess_bits=adj,
                                )
                            )
    write_records("mi_null.tsv", outputs)


if __name__ == "__main__":
    main()
