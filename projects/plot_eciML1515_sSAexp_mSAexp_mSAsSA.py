from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

# =========================
# paths
# =========================
folder = Path(__file__).resolve().parent

species = [
    ("E. coli", "eciML1515_growthResults.csv"),
    ("S. cerevisiae", "ecYeast_growthResults.csv"),
    ("C. glutamicum", "eciCW773_growthResults.csv"),
]

# =========================
# Nature-like colors
# =========================
blue = "#4C78A8"      # sSA >= exp
green = "#59A14F"     # mSA >= exp
red = "#D95F5F"       # below exp
line_gray = "#BFBFBF" # paired lines
zero_gray = "#4D4D4D" # zero line

rng = np.random.default_rng(12345)

# =========================
# read data first
# =========================
all_data = []
all_vals = []

for name, filename in species:
    df = pd.read_csv(folder / filename)
    df = df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["exp", "Growth_sSA", "Growth_mSA"]
    ).copy()

    ssa = (df["Growth_sSA"] - df["exp"]).to_numpy()
    msa = (df["Growth_mSA"] - df["exp"]).to_numpy()

    all_data.append((name, df, ssa, msa))
    all_vals.extend(ssa.tolist())
    all_vals.extend(msa.tolist())

# =========================
# make zero line at same vertical level
# but allow different y-scale in each panel
# =========================
global_min = min(all_vals)
global_max = max(all_vals)

# zero position as a fraction from the bottom
# e.g. 0.28 means zero line is at 28% of panel height
zero_frac = abs(global_min) / (abs(global_min) + global_max)

# if you want to tweak manually, you can set:
# zero_frac = 0.28

# =========================
# figure
# =========================
plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"

fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.8))
fig.patch.set_facecolor("white")

for i, (ax, (name, df, ssa, msa)) in enumerate(zip(axes, all_data)):
    ax.set_facecolor("white")

    n = len(df)
    jitter = rng.uniform(-0.035, 0.035, n)

    # paired lines
    for j in range(n):
        ax.plot(
            [1 + jitter[j], 2 + jitter[j]],
            [ssa[j], msa[j]],
            color=line_gray,
            lw=0.8,
            alpha=0.65,
            zorder=1
        )

    # sSA points
    ax.scatter(
        1 + jitter,
        ssa,
        s=36,
        c=np.where(ssa < 0, red, blue),
        edgecolors="white",
        linewidths=0.5,
        zorder=3
    )

    # mSA points
    ax.scatter(
        2 + jitter,
        msa,
        s=36,
        c=np.where(msa < 0, red, green),
        edgecolors="white",
        linewidths=0.5,
        zorder=3
    )

    # zero line
    ax.axhline(
        0,
        color=zero_gray,
        linestyle="--",
        linewidth=1.6,
        zorder=0
    )

    below_ssa = np.sum(ssa < 0)
    below_msa = np.sum(msa < 0)

    # titles
    ax.set_title(name, fontsize=16, pad=28)
    ax.text(
        0.5, 1.02,
        f"Below exp: {below_ssa}/{n} \u2192 {below_msa}/{n}",
        transform=ax.transAxes,
        ha="center", va="bottom",
        fontsize=11
    )

    ax.set_xlim(0.7, 2.3)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["sSA", "mSA"], fontsize=11)

    # only left panel has y label
    if i == 0:
        ax.set_ylabel(r"Growth capacity margin (h$^{-1}$)", fontsize=13)
    else:
        ax.set_ylabel("")

    # -------------------------
    # panel-specific y range
    # but same zero vertical position
    # -------------------------
    neg_need = max(0, -min(np.min(ssa), np.min(msa)))
    pos_need = max(0,  max(np.max(ssa), np.max(msa)))

    # add a little padding to data first
    neg_need *= 1.08
    pos_need *= 1.08

    # total span required to keep zero at same relative height
    total_span = max(
        neg_need / max(zero_frac, 1e-6),
        pos_need / max(1 - zero_frac, 1e-6)
    )

    ymin = -zero_frac * total_span
    ymax = (1 - zero_frac) * total_span

    ax.set_ylim(ymin, ymax)

    # make y ticks neat
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

    # styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10, width=1.1, length=4)

# legend
legend_elements = [
    Line2D([0], [0], marker="o", linestyle="none",
           markerfacecolor=blue, markeredgecolor="white",
           markeredgewidth=0.5, markersize=6, label="sSA ≥ exp"),

    Line2D([0], [0], marker="o", linestyle="none",
           markerfacecolor=green, markeredgecolor="white",
           markeredgewidth=0.5, markersize=6, label="mSA ≥ exp"),

    Line2D([0], [0], marker="o", linestyle="none",
           markerfacecolor=red, markeredgecolor="white",
           markeredgewidth=0.5, markersize=6, label="Below exp"),
]

fig.legend(
    handles=legend_elements,
    loc="lower center",
    ncol=3,
    frameon=False,
    bbox_to_anchor=(0.5, -0.02),
    fontsize=12,
    columnspacing=1.8,
    handletextpad=0.6
)

fig.subplots_adjust(
    left=0.07,
    right=0.99,
    top=0.82,
    bottom=0.22,
    wspace=0.34
)

out_svg = folder / "three_species_sSA_mSA_horizontal.svg"
out_png = folder / "three_species_sSA_mSA_horizontal.png"

fig.savefig(out_svg, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()