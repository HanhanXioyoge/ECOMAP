from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# =========================
# file paths
# =========================
folder = Path(__file__).resolve().parent
excel_path = folder / "Bayesian.xlsx"

out_png = folder / "Bayesian_RMSE_summary.png"
out_svg = folder / "Bayesian_RMSE_summary.svg"

# =========================
# read data
# =========================
df = pd.read_excel(excel_path, sheet_name="Sheet1")

# remove completely empty rows if any
df = df.dropna(how="all").copy()

# =========================
# expected species order
# =========================
species_order = ["E. coli", "S. cerevisiae", "C. glutamicum"]

# if your Excel contains slightly different names, edit here
species_display = {
    "E. coli": "E. coli",
    "S. cerevisiae": "S. cerevisiae",
    "C. glutamicum": "C. glutamicum"
}

# =========================
# appearance
# =========================
plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"

color_S = "#59A14F"   # green
color_M = "#F28E2B"   # orange
color_F = "#4E79A7"   # blue

line_gray = "#BFBFBF"
text_gray = "#3A3A3A"

# optional nicer model labels
def prettify_model_name(x):
    x = str(x)
    x = x.replace("ABC_", "ABC-")
    x = x.replace("_flux", " (Flux)")
    x = x.replace("_carbon", " (Carbon)")
    return x

# =========================
# helper function
# =========================
def plot_panel(ax, subdf, metric="flux", panel_title=""):
    """
    metric = 'flux' or 'carbon'
    """
    s_col = f"RMSE_S_{metric}"
    m_col = f"RMSE_M_{metric}"
    f_col = f"RMSE_F_{metric}"
    t_col = f"RMSE_Total_{metric}"

    # keep rows that at least have total RMSE for sorting
    plot_df = subdf.copy()

    # sort by Total RMSE (ascending = better at top)
    plot_df = plot_df.sort_values(by=t_col, ascending=True).reset_index(drop=True)

    # y positions
    n = len(plot_df)
    y = np.arange(n)

    # x values
    xS = plot_df[s_col].to_numpy(dtype=float)
    xM = plot_df[m_col].to_numpy(dtype=float)
    xF = plot_df[f_col].to_numpy(dtype=float)

    # determine x range from S/M/F only
    x_all = np.concatenate([xS, xM, xF])
    x_all = x_all[np.isfinite(x_all)]

    if len(x_all) == 0:
        xmin, xmax = 0, 1
    else:
        xmin = np.min(x_all)
        xmax = np.max(x_all)
        if xmin == xmax:
            xmin = xmin - 0.1
            xmax = xmax + 0.1
        else:
            pad = 0.08 * (xmax - xmin)
            xmin -= pad
            xmax += pad

    # plot each row
    for i in range(n):
        vals = [xS[i], xM[i], xF[i]]
        cols = [color_S, color_M, color_F]

        # light connecting line between the three dots
        finite_idx = np.isfinite(vals)
        finite_vals = np.array(vals)[finite_idx]
        if finite_vals.size >= 2:
            ax.plot(finite_vals, [y[i]] * len(finite_vals),
                    color=line_gray, lw=1.2, zorder=1)

        # points
        if np.isfinite(xS[i]):
            ax.scatter(xS[i], y[i], s=55, color=color_S, edgecolors="white", linewidths=0.6, zorder=3)
        if np.isfinite(xM[i]):
            ax.scatter(xM[i], y[i], s=55, color=color_M, edgecolors="white", linewidths=0.6, zorder=3)
        if np.isfinite(xF[i]):
            ax.scatter(xF[i], y[i], s=55, color=color_F, edgecolors="white", linewidths=0.6, zorder=3)

    # y labels
    ax.set_yticks(y)
    ax.set_yticklabels([prettify_model_name(m) for m in plot_df["model"]], fontsize=9)

    # invert y so best model is on top
    ax.invert_yaxis()

    # x axis
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("RMSE", fontsize=11)

    # title
    ax.set_title(panel_title, fontsize=13, pad=8)

    # grid and style
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # light separator between paired combinations (optional)
    # draws faint horizontal lines every 2 rows
    for k in [1.5, 3.5, 5.5]:
        if n > k:
            ax.axhline(k, color="#EFEFEF", linewidth=0.8, zorder=0)

# =========================
# create figure
# =========================
fig, axes = plt.subplots(3, 2, figsize=(14, 12))
fig.patch.set_facecolor("white")

for i, sp in enumerate(species_order):
    sub = df[df["Species"] == sp].copy()

    # left column: Flux RMSE
    plot_panel(
        axes[i, 0],
        sub,
        metric="flux",
        panel_title=f"{species_display[sp]} | Flux RMSE"
    )

    # right column: Carbon-weighted Flux RMSE
    plot_panel(
        axes[i, 1],
        sub,
        metric="carbon",
        panel_title=f"{species_display[sp]} | Carbon-weighted Flux RMSE"
    )

    # only left column shows y-axis label title
    axes[i, 0].set_ylabel("Models", fontsize=11)
    axes[i, 1].set_ylabel("")

# legend
legend_elements = [
    Line2D([0], [0], marker="o", linestyle="none",
           markerfacecolor=color_S, markeredgecolor="white",
           markeredgewidth=0.6, markersize=8, label="S"),
    Line2D([0], [0], marker="o", linestyle="none",
           markerfacecolor=color_M, markeredgecolor="white",
           markeredgewidth=0.6, markersize=8, label="M"),
    Line2D([0], [0], marker="o", linestyle="none",
           markerfacecolor=color_F, markeredgecolor="white",
           markeredgewidth=0.6, markersize=8, label="F"),
]

fig.legend(
    handles=legend_elements,
    loc="lower center",
    ncol=3,
    frameon=False,
    bbox_to_anchor=(0.5, 0.015),
    fontsize=11
)

fig.suptitle(
    "Bayesian calibration summary across species and RMSE definitions",
    fontsize=15,
    y=0.995
)

fig.subplots_adjust(
    left=0.22,
    right=0.98,
    top=0.95,
    bottom=0.07,
    wspace=0.25,
    hspace=0.32
)

fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(out_svg, bbox_inches="tight")
plt.show()