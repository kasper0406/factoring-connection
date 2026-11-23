"""Post-screen descriptive control diagnostics and power planning. No lead decisions."""

from scipy.stats import norm
from controls import *


def main():
    diagnostic = []
    powers = []
    for n in SIZES:
        for s in SPLITS:
            p, q, r, t = rows(n, s).T
            y, keep, pref = context(n, s)
            x = (((p % 4) * (q % 4) % 4) >> 1).astype(bool)
            target = y[2]
            assert np.array_equal(x == target, q % 4 == 1)
            for j, name in enumerate(VIEWS):
                ix = keep[2, j]
                m = int(ix.sum())
                tab = np.bincount(x[ix].astype(int) * 2 + target[ix], minlength=4)
                st = for_size(tab[None, :], m)[0]
                diagnostic.append(
                    dict(
                        n=n,
                        split=s,
                        view=name,
                        m=m,
                        q_mod4_eq1=int(((q % 4 == 1) & ix).sum()),
                        **dict(zip(FIELDS, st)),
                    )
                )
        for k in (1, 2):
            p = 0.5 + 1 / n**k
            for label, alpha in [
                ("single_one_sided", 0.05),
                ("discovery_positive_tail", THRESH / 2),
            ]:
                # Initial planning approximation, rounded upward; verify the achieved probability exactly.
                m = max(
                    1000,
                    int(
                        np.ceil(
                            (
                                0.5 * norm.isf(alpha)
                                + np.sqrt(p * (1 - p)) * norm.ppf(0.8)
                            )
                            ** 2
                            / (p - 0.5) ** 2
                            / 1000
                        )
                    )
                    * 1000,
                )
                while exact_power(m, p, alpha) < 0.8:
                    m += 1000
                cutoff = int(binom.isf(alpha, m, 0.5)) + 1
                rng = np.random.default_rng(seed("extra-plant-power", n, k, label))
                reps = 20000
                recovered = int((rng.binomial(m, p, size=reps) >= cutoff).sum())
                lo, hi = cp(np.array([recovered]), reps)
                powers.append(
                    dict(
                        n=n,
                        k=k,
                        procedure=label,
                        planned_m=m,
                        exact_power=exact_power(m, p, alpha),
                        synthetic_replicates=reps,
                        recovered=recovered,
                        simulated_power=recovered / reps,
                        sim_lo=lo[0],
                        sim_hi=hi[0],
                    )
                )
    write_rows("low1_diagnostic.tsv", diagnostic)
    write_rows("power_planning.tsv", powers)


if __name__ == "__main__":
    main()
