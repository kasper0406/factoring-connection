"""Measure synthetic agreement on the exact rows retained by each organic lock."""

from concurrent.futures import ThreadPoolExecutor
from decay import rowsfile
from inference import (
    DATA,
    RES,
    binom,
    coded,
    cp,
    jobs,
    np,
    plans,
    safe_power,
    write_records,
)


def per_size(n):
    km, manifest = plans(n)
    km = {key: spec for key, spec in km.items() if spec[0] >= 0}
    counts = {key: np.zeros(3, "i8") for key in km}
    output = []
    for cohort in ("original", "augmented"):
        files = jobs(n) if cohort == "original" else sorted(DATA.glob(f"S{n}_*.pairs"))
        subset = {
            key: spec
            for key, spec in km.items()
            if cohort == "original" or n != 96 or spec[0] < 4
        }
        for fn in files:
            _, _, _, columns = coded(n, fn, subset)
            ctrl = np.memmap(str(fn) + ".locked", "<u4", "r").reshape(-1, 4)[:, 3]
            weak = (ctrl & 16) != 0
            strong = (ctrl & 4) != 0
            for key, _, _, _, clean, _ in columns:
                counts[key] += [
                    clean.sum(),
                    (clean & strong).sum(),
                    (clean & weak).sum(),
                ]
        measured = rowsfile(RES / f"{cohort}{n}_measure.tsv")
        nlock = sum(k in km for k, _, _ in manifest)
        for rank, (key, _, target) in enumerate(manifest, 1):
            if key not in km:
                continue
            m, h1, h2 = map(int, counts[key])
            reference = next(
                r for r in measured if int(r["rank"]) == rank and r["view"] == "clean"
            )
            assert m == int(reference["m"])
            for label, delta, hits in [("1/n", 1 / n, h1), ("1/n^2", 1 / n**2, h2)]:
                lo, hi = cp(hits, m)
                output.append(
                    dict(
                        cohort=cohort,
                        n=n,
                        rank=rank,
                        key=key,
                        target=target,
                        leak=label,
                        clean=m,
                        hits=hits,
                        observed_advantage=hits / m - 0.5,
                        exact_lo=lo - 0.5,
                        exact_hi=hi - 0.5,
                        exact_one_sided_p=binom.sf(hits - 1, m, 0.5),
                        power_single=safe_power(m, delta, 0.01),
                        power_first_holm=safe_power(m, delta, 0.01 / nlock),
                    )
                )
    return output


if __name__ == "__main__":
    with ThreadPoolExecutor(3) as workers:
        data = list(workers.map(per_size, (40, 48, 56, 64, 80, 96)))
    write_records("matched_plants.tsv", [row for batch in data for row in batch])
