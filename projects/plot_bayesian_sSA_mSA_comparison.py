from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

# ============================================================
# Input files: place this script beside Bayesian_sSA.xlsx and
# Bayesian_mSA.xlsx. Both files should contain CalculatedRMSE.
# ============================================================
folder = Path(__file__).resolve().parent
input_files = {
    "sSA": folder / "Bayesian_sSA.xlsx",
    "mSA": folder / "Bayesian_mSA.xlsx",
}

out_effect_png = folder / "Bayesian_Calibration_Effect_by_StartModel.png"
out_effect_svg = folder / "Bayesian_Calibration_Effect_by_StartModel.svg"
out_start_png = folder / "Bayesian_StartModel_Comparison.png"
out_start_svg = folder / "Bayesian_StartModel_Comparison.svg"
out_calibration_summary = folder / "Bayesian_Calibration_Effect_Summary.csv"
out_start_summary = folder / "Bayesian_StartModel_Effect_Summary.csv"

# ============================================================
# Plot settings
# ============================================================
plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.unicode_minus"] = False

species_order = ["E. coli", "S. cerevisiae", "C. glutamicum"]
mode_order = ["flux", "carbon"]
combination_order = ["ABC_S", "ABC_SM", "ABC_SF", "ABC_SMF"]
start_order = ["sSA", "mSA"]
calibration_order = ["flux", "carbon"]

metric_cols = ["RMSE_S", "RMSE_M", "RMSE_F", "RMSE_Total"]
component_cols = ["RMSE_S", "RMSE_M", "RMSE_F"]

# Nature-style, colorblind-friendly colors and distinct marker shapes
component_style = {
    "RMSE_S": {"color": "#009E73", "marker": "o", "label": "S: Growth and substrate uptake"},
    "RMSE_M": {"color": "#E69F00", "marker": "s", "label": "M: Maximum growth"},
    "RMSE_F": {"color": "#0072B2", "marker": "^", "label": "F: 13C flux"},
}
line_gray = "#B8B8B8"
zero_gray = "#4D4D4D"


# ============================================================
# Load and standardize both workbooks
# ============================================================
frames = []

for start_model, path in input_files.items():
    df = pd.read_excel(path, sheet_name="CalculatedRMSE")
    df = df.dropna(how="all").copy()

    df["Species"] = df["Species"].astype(str).str.strip()
    df["model"] = df["model"].astype(str).str.strip()
    df["RMSEMode"] = df["RMSEMode"].astype(str).str.strip().str.lower()

    for col in metric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # model examples: sSA_ABC_S_flux, mSA_ABC_SMF_carbon
    df["StartModel"] = start_model
    df["Calibration"] = df["model"].str.rsplit("_", n=1).str[-1].str.lower()
    df["Combination"] = (
        df["model"]
        .str.replace(r"^(sSA|mSA)_", "", regex=True)
        .str.replace(r"_(flux|carbon)$", "", regex=True)
    )

    frames.append(df)

data = pd.concat(frames, ignore_index=True)

# Basic validation: only the expected calibration labels are plotted.
data = data[
    data["RMSEMode"].isin(mode_order)
    & data["Calibration"].isin(calibration_order)
].copy()


# ============================================================
# Helpers
# ============================================================
def combination_label(combo):
    return str(combo).replace("ABC_", "ABC-")


def calibration_label(calibration):
    return "Flux-calibrated" if calibration == "flux" else "Carbon-calibrated"


def panel_mode_label(mode):
    return "Flux RMSE evaluation" if mode == "flux" else "Carbon-weighted RMSE evaluation"


def make_component_legend():
    return [
        Line2D(
            [0], [0],
            marker=style["marker"],
            linestyle="none",
            markerfacecolor=style["color"],
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=7,
            label=style["label"],
        )
        for style in component_style.values()
    ]


def set_symmetric_xlim(ax, values, pad=0.16):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    max_abs = np.max(np.abs(values)) if values.size else 1.0
    if max_abs == 0:
        max_abs = 0.1
    limit = max_abs * (1 + pad)
    ax.set_xlim(-limit, limit)


def add_panel_style(ax):
    ax.axvline(0, color=zero_gray, linestyle="--", linewidth=1.25, zorder=0)
    ax.grid(axis="x", color="#E5E5E5", linewidth=0.75)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", labelsize=8.5, width=0.8, length=3)
    ax.tick_params(axis="y", labelsize=8.2, length=0)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))


def plot_component_deltas(ax, rows, delta_prefix, y_positions):
    # Small vertical offsets prevent S/M/F points from overlapping.
    offsets = {"RMSE_S": -0.13, "RMSE_M": 0.0, "RMSE_F": 0.13}
    for idx, (_, row) in enumerate(rows.iterrows()):
        for metric in component_cols:
            value = row.get(f"{delta_prefix}{metric}", np.nan)
            if pd.notna(value) and np.isfinite(value):
                style = component_style[metric]
                ax.scatter(
                    value,
                    y_positions[idx] + offsets[metric],
                    s=36,
                    marker=style["marker"],
                    color=style["color"],
                    edgecolors="white",
                    linewidths=0.45,
                    zorder=3,
                )


# ============================================================
# Figure 1 data:
# delta = carbon-calibrated RMSE - flux-calibrated RMSE
# ============================================================
pair_keys = ["Species", "StartModel", "Combination", "RMSEMode"]

flux_cal = (
    data[data["Calibration"] == "flux"]
    .set_index(pair_keys)[metric_cols]
)
carbon_cal = (
    data[data["Calibration"] == "carbon"]
    .set_index(pair_keys)[metric_cols]
)

calibration_delta = carbon_cal.add_suffix("_carbon").join(
    flux_cal.add_suffix("_flux"), how="inner"
).reset_index()

for metric in metric_cols:
    calibration_delta[f"Delta_{metric}"] = (
        calibration_delta[f"{metric}_carbon"] - calibration_delta[f"{metric}_flux"]
    )

# Summary for Figure 1: component-wise relative improvement, then
# equal-weight mean across available S/M/F dimensions within each
# Bayesian combination, then arithmetic mean across combinations.
calibration_summary_rows = []

for (species, mode, start), sub in calibration_delta.groupby(
    ["Species", "RMSEMode", "StartModel"], sort=False
):
    combo_scores = []
    total_scores = []

    for _, row in sub.iterrows():
        dimension_improvements = []
        for metric in component_cols:
            flux_value = row[f"{metric}_flux"]
            carbon_value = row[f"{metric}_carbon"]

            if (
                pd.notna(flux_value)
                and pd.notna(carbon_value)
                and np.isfinite(flux_value)
                and np.isfinite(carbon_value)
                and flux_value != 0
            ):
                improvement = (flux_value - carbon_value) / abs(flux_value) * 100
                dimension_improvements.append(improvement)

        if dimension_improvements:
            combo_scores.append(
                (row["Combination"], float(np.mean(dimension_improvements)))
            )

        total_flux = row["RMSE_Total_flux"]
        total_carbon = row["RMSE_Total_carbon"]
        if (
            pd.notna(total_flux) and pd.notna(total_carbon)
            and np.isfinite(total_flux) and np.isfinite(total_carbon)
            and total_flux != 0
        ):
            total_scores.append(
                (row["Combination"], (total_flux - total_carbon) / abs(total_flux) * 100)
            )

    scores = np.array([score for _, score in combo_scores], dtype=float)
    total_values = np.array([score for _, score in total_scores], dtype=float)
    n_combos = len(scores)
    if n_combos:
        calibration_summary_rows.append({
            "Species": species,
            "RMSEMode": mode,
            "StartModel": start,
            "MeanRelativeImprovement_SMF_pct": float(np.mean(scores)),
            "BetterCombinations_SMF": int(np.sum(scores > 0)),
            "NCombinations_SMF": int(n_combos),
            "MinCombinationImprovement_SMF_pct": float(np.min(scores)),
            "MaxCombinationImprovement_SMF_pct": float(np.max(scores)),
            "MeanRelativeImprovement_TotalRMSE_pct": float(np.mean(total_values)) if len(total_values) else np.nan,
            "BetterCombinations_TotalRMSE": int(np.sum(total_values > 0)) if len(total_values) else 0,
            "NCombinations_TotalRMSE": int(len(total_values)),
            "MinTotalRMSEImprovement_pct": float(np.min(total_values)) if len(total_values) else np.nan,
            "MaxTotalRMSEImprovement_pct": float(np.max(total_values)) if len(total_values) else np.nan,
        })

calibration_summary = pd.DataFrame(calibration_summary_rows)
calibration_summary.to_csv(out_calibration_summary, index=False, encoding="utf-8-sig")


# ============================================================
# Figure 1: effect of carbon-weighted calibration, stratified by
# starting model. 3 rows (species) x 2 columns (evaluation RMSE mode).
# ============================================================
fig1, axes1 = plt.subplots(
    3, 2,
    figsize=(15.5, 10.5),
    gridspec_kw={"height_ratios": [8, 8, 4]},
)
fig1.patch.set_facecolor("white")

for i, species in enumerate(species_order):
    for j, mode in enumerate(mode_order):
        ax = axes1[i, j]
        sub = calibration_delta[
            (calibration_delta["Species"] == species)
            & (calibration_delta["RMSEMode"] == mode)
        ].copy()

        # Keep Bayesian combinations in a fixed order; pair sSA/mSA within each.
        ordered_rows = []
        for combo in combination_order:
            for start in start_order:
                match = sub[
                    (sub["Combination"] == combo)
                    & (sub["StartModel"] == start)
                ]
                if not match.empty:
                    ordered_rows.append(match.iloc[0])

        plot_df = pd.DataFrame(ordered_rows).reset_index(drop=True)

        if plot_df.empty:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(f"{species} | {panel_mode_label(mode)}", fontsize=10)
            continue

        y = np.arange(len(plot_df))
        plot_component_deltas(ax, plot_df, "Delta_", y)

        labels = [
            f"{combination_label(row['Combination'])} ({row['StartModel']})"
            for _, row in plot_df.iterrows()
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_ylim(-0.5, len(plot_df) - 0.5)
        ax.invert_yaxis()

        # Subtle separators between Bayesian-combination pairs.
        for sep in range(1, len(plot_df)):
            previous_combo = plot_df.loc[sep - 1, "Combination"]
            current_combo = plot_df.loc[sep, "Combination"]
            if previous_combo != current_combo:
                ax.axhline(sep - 0.5, color="#E8E8E8", linewidth=0.8, zorder=0)

        delta_values = plot_df[
            [f"Delta_{m}" for m in component_cols]
        ].to_numpy().ravel()
        set_symmetric_xlim(ax, delta_values)
        add_panel_style(ax)

        # Show the summary separately for sSA and mSA, without mixing them.
        summary_parts = []
        for start in start_order:
            match = calibration_summary[
                (calibration_summary["Species"] == species)
                & (calibration_summary["RMSEMode"] == mode)
                & (calibration_summary["StartModel"] == start)
            ]
            if not match.empty:
                srow = match.iloc[0]
                summary_parts.append(
                    f"{start}: {srow['MeanRelativeImprovement_SMF_pct']:+.1f}% "
                    f"({int(srow['BetterCombinations_SMF'])}/{int(srow['NCombinations_SMF'])})"
                )

        summary_text = "  |  ".join(summary_parts)
        ax.set_title(
            f"{species} | {panel_mode_label(mode)}\n{summary_text}",
            fontsize=10.5,
            pad=7,
        )
        ax.set_xlabel(
            "ΔRMSE (carbon-calibrated − flux-calibrated)",
            fontsize=9.5,
        )
        if j == 0:
            ax.set_ylabel("Bayesian combination / starting model", fontsize=9.5)

fig1.legend(
    handles=make_component_legend(),
    loc="lower center",
    ncol=3,
    frameon=False,
    bbox_to_anchor=(0.5, 0.045),
    fontsize=9.5,
    columnspacing=1.8,
    handletextpad=0.6,
)
fig1.suptitle(
    "Effect of Carbon-Weighted Calibration by Starting Model",
    fontsize=15,
    y=0.985,
)
fig1.text(
    0.5, 0.012,
    "Negative ΔRMSE indicates lower error after carbon-weighted calibration. "
    "Summary: mean relative improvement across available S/M/F components; "
    "n/K indicates Bayesian combinations with positive improvement. "
    "C. glutamicum has no F data.",
    ha="center",
    va="bottom",
    fontsize=8,
    color="#555555",
)
fig1.subplots_adjust(
    left=0.22, right=0.985, top=0.91, bottom=0.12,
    wspace=0.38, hspace=0.48,
)
fig1.savefig(out_effect_png, dpi=300, bbox_inches="tight")
fig1.savefig(out_effect_svg, bbox_inches="tight")


# ============================================================
# Figure 2 data:
# delta = RMSE from mSA starting model - RMSE from sSA starting model
# ============================================================
start_keys = ["Species", "Combination", "Calibration", "RMSEMode"]

ssa_data = (
    data[data["StartModel"] == "sSA"]
    .set_index(start_keys)[metric_cols]
)
msa_data = (
    data[data["StartModel"] == "mSA"]
    .set_index(start_keys)[metric_cols]
)

start_delta = msa_data.add_suffix("_mSA").join(
    ssa_data.add_suffix("_sSA"), how="inner"
).reset_index()

for metric in metric_cols:
    start_delta[f"Delta_{metric}"] = (
        start_delta[f"{metric}_mSA"] - start_delta[f"{metric}_sSA"]
    )

# Summary for Figure 2: relative improvement when changing the starting
# model from sSA to mSA; positive means mSA has lower RMSE.
start_summary_rows = []

for (species, mode, calibration), sub in start_delta.groupby(
    ["Species", "RMSEMode", "Calibration"], sort=False
):
    combo_scores = []
    total_scores = []

    for _, row in sub.iterrows():
        dimension_improvements = []
        for metric in component_cols:
            ssa_value = row[f"{metric}_sSA"]
            msa_value = row[f"{metric}_mSA"]

            if (
                pd.notna(ssa_value)
                and pd.notna(msa_value)
                and np.isfinite(ssa_value)
                and np.isfinite(msa_value)
                and ssa_value != 0
            ):
                improvement = (ssa_value - msa_value) / abs(ssa_value) * 100
                dimension_improvements.append(improvement)

        if dimension_improvements:
            combo_scores.append(
                (row["Combination"], float(np.mean(dimension_improvements)))
            )

        total_ssa = row["RMSE_Total_sSA"]
        total_msa = row["RMSE_Total_mSA"]
        if (
            pd.notna(total_ssa) and pd.notna(total_msa)
            and np.isfinite(total_ssa) and np.isfinite(total_msa)
            and total_ssa != 0
        ):
            total_scores.append(
                (row["Combination"], (total_ssa - total_msa) / abs(total_ssa) * 100)
            )

    scores = np.array([score for _, score in combo_scores], dtype=float)
    total_values = np.array([score for _, score in total_scores], dtype=float)
    if len(scores):
        start_summary_rows.append({
            "Species": species,
            "RMSEMode": mode,
            "Calibration": calibration,
            "MeanRelativeImprovement_mSA_vs_sSA_SMF_pct": float(np.mean(scores)),
            "Combinations_mSA_better_SMF": int(np.sum(scores > 0)),
            "NCombinations_SMF": int(len(scores)),
            "MinCombinationImprovement_SMF_pct": float(np.min(scores)),
            "MaxCombinationImprovement_SMF_pct": float(np.max(scores)),
            "MeanRelativeImprovement_mSA_vs_sSA_TotalRMSE_pct": float(np.mean(total_values)) if len(total_values) else np.nan,
            "Combinations_mSA_better_TotalRMSE": int(np.sum(total_values > 0)) if len(total_values) else 0,
            "NCombinations_TotalRMSE": int(len(total_values)),
            "MinTotalRMSEImprovement_pct": float(np.min(total_values)) if len(total_values) else np.nan,
            "MaxTotalRMSEImprovement_pct": float(np.max(total_values)) if len(total_values) else np.nan,
        })

start_summary = pd.DataFrame(start_summary_rows)
start_summary.to_csv(out_start_summary, index=False, encoding="utf-8-sig")


# ============================================================
# Figure 2: direct comparison of sSA vs mSA starting models.
# Negative delta means mSA RMSE is lower; positive means sSA is lower.
# ============================================================
fig2, axes2 = plt.subplots(
    3, 2,
    figsize=(15.5, 10.5),
    gridspec_kw={"height_ratios": [8, 8, 4]},
)
fig2.patch.set_facecolor("white")

for i, species in enumerate(species_order):
    for j, mode in enumerate(mode_order):
        ax = axes2[i, j]
        sub = start_delta[
            (start_delta["Species"] == species)
            & (start_delta["RMSEMode"] == mode)
        ].copy()

        ordered_rows = []
        for combo in combination_order:
            for calibration in calibration_order:
                match = sub[
                    (sub["Combination"] == combo)
                    & (sub["Calibration"] == calibration)
                ]
                if not match.empty:
                    ordered_rows.append(match.iloc[0])

        plot_df = pd.DataFrame(ordered_rows).reset_index(drop=True)

        if plot_df.empty:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(f"{species} | {panel_mode_label(mode)}", fontsize=10)
            continue

        y = np.arange(len(plot_df))
        plot_component_deltas(ax, plot_df, "Delta_", y)

        labels = [
            f"{combination_label(row['Combination'])} "
            f"({calibration_label(row['Calibration'])})"
            for _, row in plot_df.iterrows()
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_ylim(-0.5, len(plot_df) - 0.5)
        ax.invert_yaxis()

        for sep in range(1, len(plot_df)):
            previous_combo = plot_df.loc[sep - 1, "Combination"]
            current_combo = plot_df.loc[sep, "Combination"]
            if previous_combo != current_combo:
                ax.axhline(sep - 0.5, color="#E8E8E8", linewidth=0.8, zorder=0)

        delta_values = plot_df[
            [f"Delta_{m}" for m in component_cols]
        ].to_numpy().ravel()
        set_symmetric_xlim(ax, delta_values)
        add_panel_style(ax)

        summary_parts = []
        for calibration in calibration_order:
            match = start_summary[
                (start_summary["Species"] == species)
                & (start_summary["RMSEMode"] == mode)
                & (start_summary["Calibration"] == calibration)
            ]
            if not match.empty:
                srow = match.iloc[0]
                summary_parts.append(
                    f"{calibration_label(calibration)}: "
                    f"{srow['MeanRelativeImprovement_mSA_vs_sSA_SMF_pct']:+.1f}% "
                    f"({int(srow['Combinations_mSA_better_SMF'])}/"
                    f"{int(srow['NCombinations_SMF'])})"
                )

        ax.set_title(
            f"{species} | {panel_mode_label(mode)}\n"
            + "  |  ".join(summary_parts),
            fontsize=9.7,
            pad=7,
        )
        ax.set_xlabel(
            "ΔRMSE (mSA starting model − sSA starting model)",
            fontsize=9.5,
        )
        if j == 0:
            ax.set_ylabel("Bayesian combination / calibration version", fontsize=9.5)

fig2.legend(
    handles=make_component_legend(),
    loc="lower center",
    ncol=3,
    frameon=False,
    bbox_to_anchor=(0.5, 0.045),
    fontsize=9.5,
    columnspacing=1.8,
    handletextpad=0.6,
)
fig2.suptitle(
    "Effect of Starting Model on Bayesian Calibration RMSE",
    fontsize=15,
    y=0.985,
)
fig2.text(
    0.5, 0.012,
    "Negative ΔRMSE indicates lower error for mSA; positive ΔRMSE indicates lower error for sSA. "
    "Summary values are mean relative improvements across available S/M/F components. "
    "C. glutamicum has no F data.",
    ha="center",
    va="bottom",
    fontsize=8,
    color="#555555",
)
fig2.subplots_adjust(
    left=0.24, right=0.985, top=0.91, bottom=0.12,
    wspace=0.34, hspace=0.48,
)
fig2.savefig(out_start_png, dpi=300, bbox_inches="tight")
fig2.savefig(out_start_svg, bbox_inches="tight")

plt.show()

print("Saved figures:")
print(f"  {out_effect_png}")
print(f"  {out_effect_svg}")
print(f"  {out_start_png}")
print(f"  {out_start_svg}")
print("Saved summary tables:")
print(f"  {out_calibration_summary}")
print(f"  {out_start_summary}")
