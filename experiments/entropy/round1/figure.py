from pathlib import Path
import os

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplcache")
)
os.environ.setdefault(
    "XDG_CACHE_HOME", str(Path(__file__).resolve().parent / ".mplcache")
)
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from scoring import *


def main():
    sizes = np.asarray(SIZES)
    res = ROOT / "results"
    data = []
    for n in sizes:
        d = np.load(res / f"{n}_C_counts.npz")
        pair = np.flatnonzero((d["pairs"] == [386, 4]).all(1))[0]
        st = np.load(res / f"{n}_C_stats.npz")["stats"][pair]
        data.append(st)
    a = np.asarray(data)
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.axisbelow": True,
            "svg.fonttype": "none",
        }
    )
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 3.8), layout="constrained")
    blue = "#165c77"
    orange = "#bc592c"
    for v, offset, color, marker, label in [
        (2, -0.8, blue, "o", "All matched rows"),
        (3, 0.8, orange, "s", "Both members unflagged"),
    ]:
        av = -a[:, v, 0]
        lo = -a[:, v, 2]
        hi = -a[:, v, 1]
        ok = np.isfinite(av)
        ax.errorbar(
            sizes[ok] + offset,
            av[ok],
            yerr=[av[ok] - lo[ok], hi[ok] - av[ok]],
            color=color,
            fmt=marker,
            ms=4,
            lw=1.2,
            capsize=3,
            label=label,
        )
    ax.axhline(0, color="#777", lw=0.9)
    ax.set(ylabel="Confirmation advantage", title="Strongest 24-bit discovery hit")
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.legend(frameon=False, loc="upper right")
    rows = list(csv.DictReader(open(res / "planted.tsv"), delimiter="\t"))
    for k, color, marker, label in [(1, blue, "o", "1/n"), (2, orange, "s", "1/n²")]:
        vals = [
            d
            for d in rows
            if d["split"] == "D"
            and d["view"] == "matched_clean"
            and int(d["k"]) == k
            and int(d["m"]) > 0
        ]
        bx.plot(
            [int(d["n"]) for d in vals],
            [float(d["power_discovery_positive"]) for d in vals],
            color=color,
            marker=marker,
            ms=4,
            lw=1.4,
            label=label,
        )
    bx.set_yscale("log")
    bx.set(
        title="Power at the discovery cutoff",
        ylabel="Probability of detecting the planted direction",
        ylim=(1e-9, 2),
    )
    bx.set_yticks(
        [1e-8, 1e-6, 1e-4, 1e-2, 1], labels=["10⁻⁸", "10⁻⁶", "10⁻⁴", "0.01", "1"]
    )
    bx.legend(frameon=False, title="Planted advantage", loc="center right")
    for axis in (ax, bx):
        axis.set_xticks(sizes)
        axis.set_xlabel("Product bit length")
        axis.grid(axis="y", color="#d7dfe2", lw=0.7)
    fig.savefig(res / "confirmation_and_power.png", dpi=180)


if __name__ == "__main__":
    main()
