from decay import rowsfile
from inference import *
from labels import BASE, PITMASK
from scipy.stats import binomtest
import tempfile, hmac, sys

sys.path.append("/Volumes/git/factoring-connection/.venv/lib/python3.12/site-packages")
import mpmath as mp


def main():
    checks = dict()
    total = 0
    for n in (40, 48, 56):
        prime_sets = []
        for si, s in enumerate("DVC"):
            a = load(n, s)
            ps = np.unique(a)
            prime_sets.append(ps)
            assert all(
                (lambda t: 0 if t < 19661 else 1 if t < 39322 else 2)(
                    int.from_bytes(
                        hmac.digest(
                            b"ent-c-R2-v1/prime-pool",
                            int(p).to_bytes(8, "little"),
                            "sha256",
                        )[:2],
                        "big",
                    )
                )
                == si
                for p in ps
            )
            if s != "C":
                v = np.ascontiguousarray(a).view("V16").ravel()
                assert len(np.unique(v)) == len(v)
                archive = np.load(DATA / f"{n}_{s}.labels.npz")
                pit = archive["pit"]
                bits = np.unpackbits(archive["y"], axis=1)
                coin = np.unpackbits(archive["coin"])
                h = n // 2
                lo = 3 << (h - 2)
                hi = 1 << h
                with mp.workdps(75):
                    for i in np.linspace(0, len(a) - 1, 8).astype(int):
                        p, q = map(int, a[i])
                        N = p * q
                        L = mp.log(N)
                        lower = mp.log(max(mp.mpf(lo), mp.mpf(N) / hi))
                        upper = L / 2

                        def kernel(z):
                            return 1 / (z * (L - z))

                        u = mp.quad(kernel, [lower, mp.log(p)]) / mp.quad(
                            kernel, [lower, upper]
                        )
                        assert min(65535, int(65536 * u)) == pit[i]
                        for offset, prob in [
                            (0, u),
                            (len(BASE), 1 - u / 2 if coin[i] else u / 2),
                        ]:
                            w = min(65535, int(65536 * prob))
                            for j in range(4):
                                assert bits[offset + j, i] == (w >> (15 - j) & 1)
                                assert bits[offset + 12 + j, i] == (
                                    (w & PITMASK[j]).bit_count() & 1
                                )
                        total += 1
        assert all(
            not len(np.intersect1d(prime_sets[i], prime_sets[j]))
            for i, j in ((0, 1), (1, 2), (0, 2))
        )
    checks["independent_PIT_quadratures"] = total
    excerpts = 0
    for n in (64, 80, 96):
        for supp in [False, True] if n != 80 else [False]:
            fn = DATA / (f"S{n}_000.pairs" if supp else f"C{n}_000.pairs")
            a = np.memmap(fn, "<u8", "r").reshape(-1, 2)
            with tempfile.TemporaryDirectory(dir=DATA) as tmp:
                sample = Path(tmp) / "sample"
                out = Path(tmp) / "output"
                seedname = "power-addendum-0" if supp else "chunk0"
                run(
                    ["./sampler", n, 2, 2048, seedname, sample]
                    + (["equal8"] if supp and n == 96 else [])
                )
                assert np.array_equal(
                    np.fromfile(sample, "<u8").reshape(-1, 2), a[:2048]
                )
                cfgfile = DATA / ("sup96.cfg" if supp and n == 96 else f"lock{n}.cfg")
                run(["./locked", n, fn.name, cfgfile, sample, out, "1"])
                assert (
                    out.read_bytes()
                    == np.memmap(str(fn) + ".locked", "u1", "r")[: 2048 * 16].tobytes()
                )
                excerpts += 2048
    checks["regenerated_high_size_inputs_outcomes"] = excerpts
    cis = 0
    for cohort in ("original", "augmented"):
        obs = {}
        for n in (40, 48, 56, 64, 80, 96):
            for d in rowsfile(RES / f"{cohort}{n}_measure.tsv"):
                m = int(d["m"])
                hit = int(d["n00"]) + int(d["n11"])
                i = binomtest(hit, m).proportion_ci()
                assert np.allclose(
                    [i.low - 0.5, i.high - 0.5],
                    [float(d["lo"]), float(d["hi"])],
                    atol=2e-10,
                    rtol=0,
                )
                assert (
                    float(d["centered_exact_lo"]) - 1e-12
                    <= float(d["centered"])
                    <= float(d["centered_exact_hi"]) + 1e-12
                )
                obs[(d["n"], d["rank"], d["view"])] = float(d["signed_centered"])
                cis += 1
        for row in rowsfile(RES / f"{cohort}_conditional.tsv"):
            av = obs[row["n"], row["rank"], row["view"]]
            assert abs(av - float(row["signed_centered"])) < 1e-12
            assert (
                float(row["common_odds_conditional_adv_lo"])
                < av
                < float(row["common_odds_conditional_adv_hi"])
            )
    checks["exact_CI_comparisons"] = cis
    (RES / "end_to_end_checks.json").write_text(json.dumps(checks, indent=2))
    print(checks)


if __name__ == "__main__":
    main()
