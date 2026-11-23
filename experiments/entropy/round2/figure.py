import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(__file__).resolve().parent.parent / ".mplcache")
)
os.environ.setdefault("XDG_CACHE_HOME", os.environ["MPLCONFIGDIR"])
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from decay import rowsfile
from config import *


def main():
    table = rowsfile(RES / "confirmation.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0), layout="constrained")
    plt.rcParams.update({"font.size": 10})
    a, b = axes
    for rank, off, col, mark, label in [
        ("1", -0.6, "#125a70", "o", "Even-index count parity → (2|p)"),
        ("3", 0, "#be683d", "s", "Low-byte parity → (2|p)"),
        ("5", 0.6, "#725487", "^", "Jacobi readout → p mod 3"),
    ]:
        rr = sorted(
            [
                r
                for r in table
                if r["cohort"] == "augmented"
                and r["rank"] == rank
                and r["view"] == "clean"
                and int(r["n"]) >= 48
            ],
            key=lambda r: int(r["n"]),
        )
        n = np.array([int(r["n"]) for r in rr])
        factor = 2.0 ** (n / 4)
        y = np.array([float(r["signed_centered"]) for r in rr]) * factor
        se = np.array([float(r["mh_se"]) for r in rr]) * factor
        a.errorbar(
            n + off,
            y,
            yerr=1.96 * se,
            color=col,
            marker=mark,
            ms=4,
            capsize=2,
            ls="none",
            label=label,
        )
    a.set(
        title="Three distinct locked comparisons",
        xlabel="Product bit length",
        ylabel="Signed centered advantage × $2^{n/4}$",
        xticks=[48, 56, 64, 80, 96],
    )
    a.axhline(0, color="#777", lw=0.8)
    a.legend(frameon=False, fontsize=8, loc="center right", bbox_to_anchor=(0.99, 0.43))
    for cohort, off, col, marker, label in [
        ("original", -0.65, "#72808b", "o", "Original cohort"),
        ("augmented", 0.65, "#14665d", "D", "After eligible-row top-up"),
    ]:
        rr = sorted(
            [
                r
                for r in table
                if r["cohort"] == cohort
                and r["key"] == "plant1sq"
                and r["view"] == "clean"
                and int(r["n"]) >= 64
                and not (cohort == "augmented" and int(r["n"]) == 80)
            ],
            key=lambda r: int(r["n"]),
        )
        n = np.array([int(r["n"]) for r in rr])
        factor = n * n
        value = np.array([float(r["adv"]) for r in rr]) * factor
        lo = np.array([float(r["lo"]) for r in rr]) * factor
        hi = np.array([float(r["hi"]) for r in rr]) * factor
        b.errorbar(
            n + off,
            value,
            yerr=[value - lo, hi - value],
            color=col,
            marker=marker,
            ms=4,
            capsize=3,
            ls="none",
            label=label,
        )
    b.axhline(1, lw=1.1, color="#333", ls="--", label="Planted value")
    b.axhline(0, lw=0.8, color="#777")
    b.set(
        title="Recovery of the $1/n^2$ plant",
        xlabel="Product bit length",
        ylabel="Measured agreement advantage × $n^2$",
        xticks=[64, 80, 96],
    )
    b.legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(axis="y", color="#dce2e4", lw=0.6)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(RES / "decay_and_power.png", dpi=180)


if __name__ == "__main__":
    main()
