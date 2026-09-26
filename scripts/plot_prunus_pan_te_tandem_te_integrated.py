#!/usr/bin/env python3
"""Plot an integrated TE--tandem-repeat Figure 1."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


TARGET_ORDER = ["pziA1", "pziA2", "pziB1", "pziB2", "pspA1", "pspA2", "pspB1", "pspB2", "ppaA", "ppaB", "pcam", "ppe"]
GROUPS = ["PZI", "PSP", "PPA", "PCAM", "PPE"]
TE_CATEGORIES = ["LTR/Gypsy", "LTR/Copia", "Other LTR", "TIR/DNA", "Helitron", "LINE", "SINE", "Other TE"]
# Deep, restrained palette for a high-impact journal figure.  Classes are
# separated primarily by hue (not by near-identical blues), with sufficient
# lightness contrast for grayscale reproduction.  Avoid red/green as the only
# distinction.
TE_COLORS = {
    "LTR/Gypsy": "#214E70",   # deep navy
    "LTR/Copia": "#743C5C",   # wine
    "Other LTR": "#9A6680",   # muted plum
    "TIR/DNA": "#176B6E",     # deep teal
    "Helitron": "#6F7F3B",    # olive
    "LINE": "#A65E2E",        # burnt orange
    "SINE": "#C0834A",        # copper
    "Other TE": "#5D536F",    # charcoal violet
}
ACCOUNT_COLORS = {
    "trf_pan_te_overlap_bp": "#713B5A",  # shared/overlapping sequence
    "trf_only_bp": "#B87928",             # tandem-repeat-specific
    "pan_te_only_bp": "#176B6E",           # TE-specific
    "neither_bp": "#707070",              # unclassified/background
}
PERIOD_COLORS = {
    "1": "#B8B8B8", "2-6": "#486A83", "7-20": "#2F7F83",
    "21-100": "#6F7F3B", "101-500": "#B87928",
    "501-2000": "#8B4A45", ">2000": "#5D536F",
}

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.75, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42, "figure.facecolor": "white", "axes.facecolor": "white"})


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(fig, output_dir, stem):
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix, kwargs in (("svg", {}), ("pdf", {}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}}), ("png", {"dpi": 300})):
        fig.savefig(output_dir / f"{stem}.{suffix}", bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def label(ax, value):
    ax.text(-0.12, 1.04, value, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def te_field(category):
    return f"te_{category.replace('/', '_').replace(' ', '_')}_bp"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrated-dir", type=Path, required=True)
    parser.add_argument("--te-composition", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    target_rows = {row["target"]: row for row in read_tsv(args.integrated_dir / "integrated_tandem_te_target.tsv")}
    groups = {row["group"]: row for row in read_tsv(args.integrated_dir / "integrated_tandem_te_group.tsv")}
    contrasts = read_tsv(args.integrated_dir / "integrated_tandem_te_contrasts.tsv")
    periods = read_tsv(args.integrated_dir / "integrated_tandem_period_group.tsv")
    te_rows = read_tsv(args.te_composition)
    te_by_target = {}
    for row in te_rows:
        te_by_target.setdefault(row["target"], {})[row["category"]] = int(row["bp"]) / 1e6

    fig = plt.figure(figsize=(7.25, 6.35))
    grid = fig.add_gridspec(2, 2, width_ratios=(1.25, 1.0), height_ratios=(1.05, 0.95), wspace=0.44, hspace=0.56)

    ax = fig.add_subplot(grid[0, 0])
    y = np.arange(len(TARGET_ORDER)); left = np.zeros(len(TARGET_ORDER))
    for category in TE_CATEGORIES + ["Non-TE"]:
        values = np.array([te_by_target[target].get(category, 0.0) for target in TARGET_ORDER])
        color = "#D5D5D5" if category == "Non-TE" else TE_COLORS[category]
        ax.barh(y, values, left=left, height=0.68, color=color, label=category, linewidth=0)
        left += values
    for i, target_id in enumerate(TARGET_ORDER):
        ax.text(left[i] + 3, i, f"{left[i]:.0f}", va="center", fontsize=5.6)
    labels = [target.replace("pzi", "PZI-").replace("psp", "PSP-").replace("ppa", "PPA-").upper() if target not in ("pcam", "ppe") else target.upper() for target in TARGET_ORDER]
    ax.set_yticks(y, labels); ax.invert_yaxis(); ax.set_xlabel("Genome package size (Mb)"); ax.set_xlim(0, max(left) * 1.16); ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    for boundary in (3.5, 7.5, 9.5, 10.5): ax.axhline(boundary, color="#D9D9D9", lw=0.6)
    ax.legend(title="RepeatMasker composition", title_fontsize=5.7, fontsize=5.0, ncol=2, columnspacing=0.8, handlelength=1.0, loc="lower left", bbox_to_anchor=(0, 1.02))
    label(ax, "a")

    ax = fig.add_subplot(grid[0, 1])
    y = np.arange(len(TARGET_ORDER)); left = np.zeros(len(TARGET_ORDER))
    for field, name in (("trf_pan_te_overlap_bp", "TRF ∩ pan-TE"), ("trf_only_bp", "TRF only"), ("pan_te_only_bp", "pan-TE only"), ("neither_bp", "Neither")):
        values = np.array([int(target_rows[t][field]) / 1e6 for t in TARGET_ORDER])
        ax.barh(y, values, left=left, height=0.68, color=ACCOUNT_COLORS[field], label=name, linewidth=0)
        left += values
    ax.set_yticks(y, labels); ax.invert_yaxis(); ax.set_xlabel("Mutually exclusive sequence (Mb)"); ax.set_xlim(0, max(left) * 1.16); ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    for boundary in (3.5, 7.5, 9.5, 10.5): ax.axhline(boundary, color="#D9D9D9", lw=0.6)
    ax.legend(fontsize=5.4, ncol=1, loc="lower left", bbox_to_anchor=(0, 1.02)); label(ax, "b")
    ax.text(0.01, -0.25, "Overlap is assigned once; categories sum to genome size.", transform=ax.transAxes, fontsize=5.1, color="#555555")

    ax = fig.add_subplot(grid[1, 0])
    x = np.arange(len(GROUPS)); width = 0.25
    for offset, field, name, color in ((-width, "pan_te_total_fraction", "pan-TE total", ACCOUNT_COLORS["pan_te_only_bp"]), (0, "trf_total_fraction", "TRF total", ACCOUNT_COLORS["trf_only_bp"]), (width, "trf_pan_te_overlap_fraction", "TRF ∩ pan-TE", ACCOUNT_COLORS["trf_pan_te_overlap_bp"])):
        vals = [100 * float(groups[group][f"mean_{field}"]) for group in GROUPS]
        ax.bar(x + offset, vals, width=width, color=color, label=name)
    ax.set_xticks(x, GROUPS); ax.set_ylabel("Mean fraction of genome (%)"); ax.set_ylim(0, 100); ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True); ax.legend(fontsize=5.2, ncol=1, loc="upper right"); label(ax, "c")
    ax.text(0.01, -0.24, "TRF total and pan-TE total are overlapping views, not additive components.", transform=ax.transAxes, fontsize=5.0, color="#555555")

    ax = fig.add_subplot(grid[1, 1])
    # Signed decomposition across all comparator genomes. Positive segments
    # indicate sequence excess in PZI; negative segments indicate components
    # that are smaller in PZI than in the comparator.
    fields = [
        ("trf_pan_te_overlap_difference_bp", "TRF ∩ pan-TE", ACCOUNT_COLORS["trf_pan_te_overlap_bp"]),
        ("trf_only_difference_bp", "TRF only", ACCOUNT_COLORS["trf_only_bp"]),
        ("pan_te_only_difference_bp", "pan-TE only", ACCOUNT_COLORS["pan_te_only_bp"]),
        ("neither_difference_bp", "Neither", ACCOUNT_COLORS["neither_bp"]),
    ]
    comparator_order = ["PSP", "PPA", "PCAM", "PPE"]
    y = np.arange(len(comparator_order))
    for yi, comparator in zip(y, comparator_order):
        contrast = next(row for row in contrasts if row["comparator_group"] == comparator)
        positive = 0.0
        negative = 0.0
        for field, name, color in fields:
            value = int(contrast[field]) / 1e6
            if value >= 0:
                ax.barh(yi, value, left=positive, color=color, height=0.58, label=name if yi == 0 else None, linewidth=0)
                positive += value
            else:
                ax.barh(yi, value, left=negative, color=color, height=0.58, label=name if yi == 0 else None, linewidth=0)
                negative += value
        net = int(contrast["genome_size_difference_bp"]) / 1e6
        ax.text(max(positive, 0) + 4, yi, f"Δ {net:.0f}", va="center", fontsize=5.6)
    ax.axvline(0, color="#222222", lw=0.8)
    ax.set_yticks(y, [f"PZI − {group}" for group in comparator_order]); ax.invert_yaxis()
    ax.set_xlabel("Mean sequence difference (Mb)"); ax.set_xlim(-45, 235); ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True); ax.legend(fontsize=5.0, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.02)); label(ax, "d")
    ax.text(0.01, -0.27, "Positive values are PZI excess; negative values are comparator excess.", transform=ax.transAxes, fontsize=5.0, color="#555555")

    save(fig, args.output_dir, "figure1_integrated_te_tandem")

    # A compact supplementary period panel remains useful and avoids crowding
    # the main composition panels with a second taxonomy.
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    period_order = ["1", "2-6", "7-20", "21-100", "101-500", "501-2000", ">2000"]
    x = np.arange(len(GROUPS)); bottom = np.zeros(len(GROUPS))
    for period in period_order:
        vals = np.array([int(next((row["mean_exclusive_union_bp"] for row in periods if row["group"] == group and row["period_class"] == period), 0)) / 1e6 for group in GROUPS])
        ax.bar(x, vals, bottom=bottom, width=0.65, color=PERIOD_COLORS[period], label=period); bottom += vals
    ax.set_xticks(x, GROUPS); ax.set_ylabel("Mean non-redundant TRF span (Mb)"); ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True); ax.legend(title="TRF period (bp)", title_fontsize=5.5, fontsize=5.0, ncol=2)
    save(fig, args.output_dir, "supp_integrated_trf_period")
    print(f"Figures written to {args.output_dir}")


if __name__ == "__main__":
    main()
