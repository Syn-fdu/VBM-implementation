#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced plot for Extension 02:
VBM under propensity-score misspecification versus omitted confounding.

Outputs:
  1) enhanced_vbm_ps_misspec_decomposition.png
     - the original-style R² decomposition figure;
  2) enhanced_vbm_ps_misspec_vbm_response_table.png
     - a compact table recording how VBM responds in each generated dataset.

The script is intentionally read-only with respect to the R outputs.  It only
consumes CSV files under the extension table directory and writes figures to
the enhanced-figure directory, so it does not disturb the rest of the project.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ORDER = [
    "no_disturbance",
    "model_misspecification_only",
    "confounding_only",
]
LABELS = ["No disturbance", "Model error only", "Confounding only"]

FIGSIZE = (6.8, 3.35)
TABLE_FIGSIZE = (10.8, 4)
DPI = 320


def read_csv(input_dir: Path, name: str) -> pd.DataFrame:
    path = input_dir / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def save(fig: plt.Figure, output_dir: Path, name: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_dir / name, dpi=DPI, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def prepare_summary(summary: pd.DataFrame) -> pd.DataFrame:
    summary = summary.copy()
    summary["scenario"] = pd.Categorical(summary["scenario"], ORDER, ordered=True)
    summary = summary.sort_values("scenario").reset_index(drop=True)

    hidden_present = summary["hidden_present"].astype(bool)
    ps_misspec = summary["ps_setting"].astype(str).eq("misspecified")

    summary["hidden_display"] = np.where(
        hidden_present,
        summary["hidden_only_R2"].astype(float),
        np.nan,
    )
    summary["ps_display"] = np.where(
        ps_misspec,
        summary["misspec_only_R2"].astype(float),
        np.nan,
    )
    summary["active_display"] = summary["active_error_R2"].astype(float)

    return summary


def make_decomposition_figure(summary: pd.DataFrame, output_dir: Path) -> None:
    summary = prepare_summary(summary)

    x = np.arange(len(summary), dtype=float)
    width = 0.20

    hidden_vals = summary["hidden_display"].to_numpy(dtype=float)
    ps_vals = summary["ps_display"].to_numpy(dtype=float)
    active_vals = summary["active_display"].to_numpy(dtype=float)
    vbm_vals = summary["VBM_R2_star"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    bars_hidden = ax.bar(
        x - width,
        hidden_vals,
        width=width,
        label="Hidden confounding",
        alpha=0.90,
    )
    bars_ps = ax.bar(
        x,
        ps_vals,
        width=width,
        label="PS model error",
        alpha=0.90,
    )
    bars_active = ax.bar(
        x + width,
        active_vals,
        width=width,
        label="Active-source R²",
        alpha=0.90,
    )

    pts = ax.scatter(
        x,
        vbm_vals,
        marker="D",
        s=90,
        label="VBM R²*",
        zorder=4,
    )

    for xi, yi in zip(x, vbm_vals):
        if np.isfinite(yi):
            ax.text(
                xi,
                yi + 0.012,
                f"{yi:.2f}",
                ha="center",
                va="bottom",
                fontsize=8.5,
            )

    finite = np.r_[
        hidden_vals[np.isfinite(hidden_vals)],
        ps_vals[np.isfinite(ps_vals)],
        active_vals[np.isfinite(active_vals)],
        vbm_vals[np.isfinite(vbm_vals)],
        [0.50],
    ]
    ymax = np.nanmax(finite)
    ax.set_ylim(0, max(0.50, float(ymax) + 0.06))
    ax.set_ylabel("R²")
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS)
    ax.set_title("VBM response under model misspecification and confounding", fontsize=11.5)

    ax.grid(axis="y", linewidth=0.7, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [pts, bars_hidden, bars_ps, bars_active]
    labels = ["VBM R²*", "Hidden confounding", "PS model error", "Active-source R²"]
    ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=4,
        frameon=False,
        handletextpad=0.6,
        columnspacing=1.0,
    )

    save(fig, output_dir, "enhanced_vbm_ps_misspec_decomposition.png")


def fmt_num(x: float) -> str:
    return "NA" if not np.isfinite(x) else f"{x:.3f}"


def wrap_cell(text: str, width: int = 32) -> str:
    text = str(text)
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def make_response_table(response: pd.DataFrame, output_dir: Path) -> None:
    response = response.copy()
    response["scenario"] = pd.Categorical(response["scenario"], ORDER, ordered=True)
    response = response.sort_values("scenario").reset_index(drop=True)

    rows = []
    for _, r in response.iterrows():
        rows.append([
            r.get("scenario_label", ""),
            fmt_num(float(r.get("active_error_R2", np.nan))),
            fmt_num(float(r.get("VBM_R2_star", np.nan))),
            fmt_num(float(r.get("observed_minus_oracle_ATT", np.nan))),
            wrap_cell(r.get("VBM_response", ""), 34),
            wrap_cell("VBM cannot identify whether the perturbation came from model error or confounding.", 42),
        ])

    columns = [
        "Dataset",
        "Active R²",
        "VBM R²*",
        "ATT shift",
        "VBM response",
        "Interpretation",
    ]

    fig, ax = plt.subplots(figsize=TABLE_FIGSIZE)
    ax.axis("off")
    table = ax.table(
        cellText=rows,
        colLabels=columns,
        loc="center",
        cellLoc="center",
        colLoc="center",
        colWidths=[0.16, 0.09, 0.09, 0.09, 0.28, 0.29],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 2.5)

    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(0.45)
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_height(cell.get_height() * 1.15)
        if col in (4, 5) and row > 0:
            cell.set_text_props(ha="left")

    ax.set_title("Recorded VBM response by generated dataset", fontsize=12, pad=10)
    save(fig, output_dir, "enhanced_vbm_ps_misspec_vbm_response_table.png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        required=True,
        help="Directory containing vbm_ps_misspec_factorial_summary.csv",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for the enhanced plot and response table",
    )
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)

    summary = read_csv(input_dir, "vbm_ps_misspec_factorial_summary.csv")
    make_decomposition_figure(summary, output_dir)

    response_path = input_dir / "vbm_ps_misspec_vbm_response_table.csv"
    if response_path.exists():
        response = pd.read_csv(response_path)
    else:
        # Backward-compatible fallback for old R outputs.
        response = summary.copy()
        response["VBM_response"] = "VBM response table was not produced by the R script."
    make_response_table(response, output_dir)

    print(f"Saved enhanced Extension 02 figure and VBM response table to {output_dir}")


if __name__ == "__main__":
    main()
