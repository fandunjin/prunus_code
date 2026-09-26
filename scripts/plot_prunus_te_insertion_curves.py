#!/usr/bin/env python3
"""Plot standardized TE divergence and independent intact-LTR bars."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


TARGET_ORDER = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]
GROUP = {
    "pziA1": "PZI", "pziA2": "PZI", "pziB1": "PZI", "pziB2": "PZI",
    "pspA1": "PSP", "pspA2": "PSP", "pspB1": "PSP", "pspB2": "PSP",
    "ppaA": "PPA", "ppaB": "PPA", "pcam": "PCAM", "ppe": "PPE",
}
COLORS = {"PZI": "#C44E52", "PSP": "#2A9D8F", "PPA": "#4C78A8", "PCAM": "#777777", "PPE": "#B279A2"}
STYLE = {"PZI": "-", "PSP": "-", "PPA": "--", "PCAM": ":", "PPE": "-."}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8,
    "axes.linewidth": 0.75,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "legend.frameon": False,
})


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(fig, stem):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--divergence", type=Path, required=True)
    parser.add_argument("--intact-ltr", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-divergence", type=float, default=20.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_te = defaultdict(lambda: defaultdict(float))
    for row in read_tsv(args.divergence):
        target = row["target"]
        if target not in TARGET_ORDER:
            continue
        x = float(row["divergence_bin_mid_pct"])
        if 0 <= x <= args.max_divergence:
            all_te[target][x] += float(row["exclusive_bp"]) / 1e6

    intact = defaultdict(lambda: defaultdict(float))
    for row in read_tsv(args.intact_ltr):
        target = row["target"]
        if target not in TARGET_ORDER:
            continue
        # Bin intact-LTR divergence to the same 1% resolution as the
        # RepeatMasker landscape, avoiding visually misleading single-copy
        # spikes from exact identity values.
        x = np.floor(float(row["ltr_divergence"]) * 100.0) + 0.5
        if 0 <= x <= args.max_divergence:
            intact[target][x] += float(row["length_bp"]) / 1e6

    # One shared x-axis: lines use the left y-axis and grouped bars use the
    # right y-axis. This keeps the two measurements aligned without a second
    # long x-axis or a visually misleading log-scale line plot.
    fig, ax = plt.subplots(figsize=(7.6, 4.25))
    ax_ltr = ax.twinx()
    groups = ["PZI", "PSP", "PPA", "PCAM", "PPE"]
    x_bins = np.arange(0.5, min(args.max_divergence, 20.0) + 0.01, 1.0)

    # Intact LTR-RT bars are plotted first so the standardized TE landscapes
    # remain visually dominant. Values are averaged within each species.
    bar_width = 0.15
    offsets = (np.arange(len(groups)) - (len(groups) - 1) / 2.0) * bar_width
    for offset, group in zip(offsets, groups):
        members = [t for t in TARGET_ORDER if GROUP[t] == group and t in intact]
        if not members:
            continue
        vals = np.array([[intact[t].get(x, 0.0) for x in x_bins] for t in members])
        ax_ltr.bar(
            x_bins + offset, vals.mean(axis=0), width=bar_width * 0.92,
            color=COLORS[group], alpha=0.58, edgecolor=COLORS[group], linewidth=0.25,
            label=f"{group} LTR",
            zorder=1,
        )

    for group in groups:
        members = [t for t in TARGET_ORDER if GROUP[t] == group and t in all_te]
        if not members:
            continue
        xs = sorted({x for t in members for x in all_te[t]})
        if not xs:
            continue
        vals = np.array([[all_te[t].get(x, 0.0) for x in xs] for t in members])
        mean = vals.mean(axis=0)
        ax.plot(xs, mean, color=COLORS[group], linestyle=STYLE[group], linewidth=1.8, label=group, zorder=3)
        if len(members) > 1 and group in {"PZI", "PSP"}:
            ax.fill_between(xs, vals.min(axis=0), vals.max(axis=0), color=COLORS[group], alpha=0.10, linewidth=0, zorder=2)

    ax.set_xlim(0, min(args.max_divergence, 20.0))
    ax.set_xlabel("Sequence divergence (%)")
    ax.set_ylabel("Exclusive TE sequence (Mb per 1% bin)")
    ax_ltr.set_ylabel("Intact LTR-RT sequence (Mb)")
    # Keep the plotting field clean; the axes and ticks provide the scale.
    ax.grid(False)
    ax.set_facecolor("white")
    ax_ltr.set_facecolor("none")
    fig.patch.set_facecolor("white")
    ax.tick_params(width=0.7, length=3)
    ax_ltr.tick_params(width=0.7, length=3)
    ax.spines["top"].set_visible(False)
    ax_ltr.spines["top"].set_visible(False)
    ax_ltr.spines["right"].set_visible(True)
    ax_ltr.spines["right"].set_color("black")
    ax_ltr.spines["right"].set_linewidth(0.75)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, ncol=3, loc="upper right", fontsize=7, handlelength=2.4)
    fig.tight_layout()
    stem = args.output_dir / "prunus_te_insertion_curves"
    save(fig, stem)
    print(f"Wrote {stem}.pdf/.svg/.png/.tiff")


if __name__ == "__main__":
    main()
