from decay import rowsfile
from inference import *
from scipy.stats import binom


def export():
    ns = (40, 48, 56, 64, 80, 96)
    maskrows = []
    allrows = []
    powrows = []
    fits = []
    cond = []
    cluster = []
    for cohort in ("original", "augmented"):
        for n in ns:
            allrows.extend(rowsfile(RES / f"{cohort}{n}_measure.tsv"))
            powrows.extend(rowsfile(RES / f"{cohort}{n}_power.tsv"))
            cc = cohort if n in (64, 96) else "original"
            cluster.extend(
                [dict(x, cohort=cohort) for x in rowsfile(RES / f"{cc}{n}_dyad.tsv")]
            )
            base = np.load(RES / f"C{n}_masks.npz")
            km, ar = plans(n)
            masks = {k: t.copy() for k, t in zip(base["keys"], base["by_measurement"])}
            if cohort == "augmented" and n in (64, 96):
                more = np.load(RES / f"S{n}_masks.npz")
                fiber = None
                if n == 96:
                    fiber = np.zeros(32, "i8")
                    for fn in jobs(n):
                        code = np.memmap(str(fn) + ".locked", "<u4", "r").reshape(
                            -1, 4
                        )[:, 0]
                        b = np.memmap(str(fn) + ".budget", "u1", "r")
                        ok = ((b & 48) == 0) & (RESIDUES[code // 1024] % 8 == 1)
                        fiber[:16] += np.bincount(b[ok] & 15, minlength=16)
                for k, p in zip(more["keys"], more["by_measurement"]):
                    if n == 96 and k == "negative":
                        continue
                    if n == 96 and k.startswith("plant"):
                        masks[k] = fiber.copy()
                    masks[k] += p
            for i, (key, sign, target) in enumerate(ar):
                m = masks[key]
                idx = np.arange(32)
                row = dict(
                    cohort=cohort,
                    n=n,
                    rank=i + 1 if key.startswith("t") else 0,
                    key=key,
                    target=target,
                    eligible=int(m.sum()),
                    clean=int(m[0]),
                    any_split=int(m.sum() - m[0]),
                    fraction_split=1 - m[0] / m.sum(),
                )
                for j, name in enumerate(
                    [
                        "close",
                        "Lehman",
                        "minus_64_smooth",
                        "plus_64_smooth",
                        "gate_proper_gcd",
                    ]
                ):
                    row[name] = int(m[(idx & (1 << j)) > 0].sum())
                maskrows.append(row)
                for v, z in [("clean", row["clean"]), ("all", row["eligible"])]:
                    ms = [
                        d
                        for d in allrows
                        if d["cohort"] == cohort
                        and int(d["n"]) == n
                        and d["target"] == target
                        and d["key"] == key
                        and d["view"] == v
                    ]
                    assert len(ms) == 1 and int(ms[0]["m"]) == z, (
                        cohort,
                        n,
                        key,
                        v,
                        z,
                        ms,
                    )
        fits.extend(rowsfile(RES / f"{cohort}_decay.tsv"))
        cond.extend(rowsfile(RES / f"{cohort}_conditional.tsv"))
    raw_tails = {}
    for row in allrows:
        m = int(row["m"])
        hits = int(row["n00"]) + int(row["n11"])
        if int(row["sign"]) < 0:
            hits = m - hits
        raw_p = float(binom.sf(hits - 1, m, 0.5))
        row["agreement_p_one_sided_exact"] = raw_p
        row["primary_p_estimand"] = (
            "centered_MH" if int(row["rank"]) else "raw_agreement"
        )
        if row["view"] == "clean":
            raw_tails[row["cohort"], row["n"], row["key"], row["target"]] = raw_p
    for row in powrows:
        row["agreement_p_one_sided_exact"] = raw_tails[
            row["cohort"], row["n"], row["key"], row["target"]
        ]
        row["primary_p_estimand"] = (
            "centered_MH" if int(row["rank"]) else "raw_agreement"
        )
    for name, data in [
        ("confirmation.tsv", allrows),
        ("power.tsv", powrows),
        ("scaling.tsv", fits),
        ("conditional_exact.tsv", cond),
        ("cluster_sensitivity.tsv", cluster),
        ("split_flags.tsv", maskrows),
    ]:
        write_records(name, data)


if __name__ == "__main__":
    export()
