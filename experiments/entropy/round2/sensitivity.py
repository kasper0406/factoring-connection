from decay import *


def main():
    for cohort in ("original", "augmented"):
        output = []
        fiber = []
        for n in (64, 80, 96):
            dy_file = RES / f"{cohort if n in (64, 96) else 'original'}{n}_dyad.tsv"
            while not dy_file.exists():
                time.sleep(3)
        for rank in (1, 3, 5):
            for drop in (False, True):
                ns = []
                Y = []
                se = []
                for n in (40, 48, 56, 64, 80, 96):
                    if n == 40 and drop:
                        continue
                    r = next(
                        r
                        for r in rowsfile(RES / f"{cohort}{n}_measure.tsv")
                        if int(r["rank"]) == rank and r["view"] == "clean"
                    )
                    cc = cohort if n in (64, 96) else "original"
                    dy = next(
                        z
                        for z in rowsfile(RES / f"{cc}{n}_dyad.tsv")
                        if z["key"] == r["key"]
                    )
                    se.append(
                        max(
                            float(r["mh_se"]),
                            float(r["centered_se"]),
                            float(dy["dyad_se"]),
                        )
                    )
                    ns.append(n)
                    Y.append(float(r["signed_centered"]))
                ns, Y, se = map(np.array, (ns, Y, se))
                fit = {k: profile(Y, ns, se, k) for k in ["poly", "wide_poly", "free"]}
                output.append(
                    dict(
                        cohort=cohort,
                        rank=rank,
                        exclude40=int(drop),
                        poly_rss=fit["poly"][0][0],
                        exp_rss=fit["free"][0][0],
                        beta=fit["free"][2][0],
                        wide_k=fit["wide_poly"][2][0],
                        wide_poly_rss=fit["wide_poly"][0][0],
                    )
                )
        for n in (64, 80, 96):
            for k, mats, _ in tables(n, cohort):
                if k != "t829_y31":
                    continue
                for residue in (1, 2):
                    sel = np.repeat(RESIDUES % 3 == residue, 1024)
                    c = mats[1][sel].astype("i8")
                    s = c.sum(1)
                    ok = s > 0
                    c = c[ok]
                    s = s[ok]
                    x = c[:, 2] + c[:, 3]
                    y = c[:, 1] + c[:, 3]
                    U = np.sum(c[:, 3] - x * y / s)
                    m = s.sum()
                    v = (x * (1 - x / s) * y * (1 - y / s) / np.maximum(s - 1, 1)).sum()
                    fiber.append(
                        dict(
                            cohort=cohort,
                            n=n,
                            N_mod3=residue,
                            m=m,
                            centered_full_strata=2 * U / m,
                            mh_z=U / np.sqrt(v) if v else 0,
                        )
                    )
        write_records(f"{cohort}_fit_sensitivity.tsv", output)
        write_records(f"{cohort}_mod3_fibers.tsv", fiber)
        print(cohort, "complete", flush=True)


if __name__ == "__main__":
    main()
