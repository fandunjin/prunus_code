#!/usr/bin/env python3
"""Plot the formal two-panel Figure 2a Ks density comparison."""

from __future__ import print_function

import argparse
import csv
import math
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})


COLORS = {
    "PZI": "#B64342",
    "PSP": "#0F4D92",
    "PPA": "#42949E",
    "PCAM": "#C98223",
    "PPE": "#6F4C9B",
}
LABELS_UPPER = OrderedDict((key, key + " A-B") for key in ("PZI", "PSP", "PPA"))
LABELS_LOWER = OrderedDict((key, key + " (A+B)-GTR") for key in ("PZI", "PSP", "PPA", "PCAM", "PPE"))
DISPLAY_UPPER = {
    "PZI": "P. zippeliana A1-B1",
    "PSP": "P. spinulosa A1-B1",
    "PPA": "P. padus A-B",
}
DISPLAY_LOWER = {
    "PZI": "P. zippeliana (A+B)-G. trifoliata",
    "PSP": "P. spinulosa (A+B)-G. trifoliata",
    "PPA": "P. padus (A+B)-G. trifoliata",
    "PCAM": "P. campanulata-G. trifoliata",
    "PPE": "P. persica-G. trifoliata",
}


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_values(path, expected):
    grouped = {key: [] for key in expected}
    for row in read_tsv(path):
        comparison = row.get("comparison", "")
        if comparison not in grouped:
            raise SystemExit("[ERROR] Unexpected comparison {} in {}".format(comparison, path))
        try:
            value = float(row["Ks"])
        except (KeyError, TypeError, ValueError):
            raise SystemExit("[ERROR] Invalid Ks value in {}: {}".format(path, row))
        if not math.isfinite(value) or not 0.0 < value < 5.0:
            raise SystemExit("[ERROR] Non-filtered Ks value in {}: {}".format(path, value))
        grouped[comparison].append(value)
    for comparison, values in grouped.items():
        if not values:
            raise SystemExit("[ERROR] No Ks values for {} in {}".format(comparison, path))
    return grouped


def reflected_gaussian_kde(values, grid, bandwidth):
    values = np.asarray(values, dtype=float)
    density = np.zeros_like(grid, dtype=float)
    denominator = len(values) * bandwidth * math.sqrt(2.0 * math.pi)
    for start in range(0, len(values), 4096):
        chunk = values[start:start + 4096, None]
        direct = (grid[None, :] - chunk) / bandwidth
        reflected = (grid[None, :] + chunk) / bandwidth
        density += np.exp(-0.5 * direct * direct).sum(axis=0)
        density += np.exp(-0.5 * reflected * reflected).sum(axis=0)
    return density / denominator


def plot_panel(ax, grouped, labels, display_labels, panel_id, ylabel, x_max, bandwidth, source_rows, peak_rows):
    grid = np.linspace(0.0, x_max, 1001)
    curves = []
    for comparison, label in labels.items():
        values = grouped[comparison]
        density = reflected_gaussian_kde(values, grid, bandwidth)
        color = COLORS[comparison]
        ax.fill_between(grid, density, color=color, alpha=0.10, linewidth=0)
        ax.plot(grid, density, lw=1.25, color=color, label="{}  n={:,}".format(display_labels[comparison], len(values)))
        mode_index = int(np.argmax(density))
        curves.append((comparison, color, density, mode_index))
        peak_rows.append({
            "panel": panel_id,
            "comparison": comparison,
            "KDE_mode_Ks": "{:.8f}".format(float(grid[mode_index])),
            "mode_density": "{:.8f}".format(float(density[mode_index])),
            "n": str(len(values)),
        })
        for x, y in zip(grid, density):
            source_rows.append({
                "panel": panel_id,
                "comparison": comparison,
                "Ks_grid": "{:.8f}".format(float(x)),
                "density": "{:.8f}".format(float(y)),
                "n": str(len(values)),
                "bandwidth": "{:.4f}".format(bandwidth),
            })
    offsets = {
        "upper": {"PZI": (0, 8), "PSP": (20, 20), "PPA": (-22, 8)},
        "lower": {"PZI": (0, 8), "PSP": (22, 20), "PPA": (-22, 8),
                  "PCAM": (32, 28), "PPE": (-34, 28)},
    }
    for comparison, color, density, mode_index in curves:
        dx, dy = offsets.get(panel_id, {}).get(comparison, (0, 8))
        ax.annotate(
            r"$K_s$={:.3f}".format(float(grid[mode_index])),
            xy=(float(grid[mode_index]), float(density[mode_index])),
            xytext=(dx, dy), textcoords="offset points",
            ha="center", va="bottom", fontsize=5.3, color=color,
            arrowprops={"arrowstyle": "-", "lw": 0.45, "color": color, "alpha": 0.85},
            annotation_clip=False,
        )
    ax.set_xlim(0.0, x_max)
    ax.set_ylabel(ylabel)
    ax.legend(loc="upper right", fontsize=5.7,
              handlelength=2.0, borderpad=0.2, labelspacing=0.35)
    ax.grid(axis="y", color="#D9D9D9", lw=0.35, alpha=0.55)
    ax.set_axisbelow(True)


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upper", type=Path, required=True)
    parser.add_argument("--lower", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bandwidth", type=float, default=0.04)
    parser.add_argument("--upper-xmax", type=float, default=0.6)
    parser.add_argument("--lower-xmax", type=float, default=1.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.overwrite:
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)
    obsolete_medians = args.outdir / "panel_a_median_lines.tsv"
    if args.overwrite and obsolete_medians.is_file():
        obsolete_medians.unlink()

    upper = parse_values(args.upper, LABELS_UPPER)
    lower_tracks = parse_values(args.lower, {"PZI", "PSP", "PPA", "PCAM", "PPE"})
    lower = OrderedDict()
    for comparison in LABELS_LOWER:
        lower[comparison] = lower_tracks[comparison]
    source_rows = []
    peak_rows = []

    shared_xmax = max(args.upper_xmax, args.lower_xmax)
    fig, axes = plt.subplots(2, 1, figsize=(3.6, 6.2), sharex=True, sharey=False)
    plot_panel(axes[0], upper, LABELS_UPPER, DISPLAY_UPPER, "upper", "Paralogue density",
               shared_xmax, args.bandwidth, source_rows, peak_rows)
    plot_panel(axes[1], lower, LABELS_LOWER, DISPLAY_LOWER, "lower", "Orthologue density",
               shared_xmax, args.bandwidth, source_rows, peak_rows)
    axes[0].tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    axes[1].set_xlabel(r"Synonymous substitutions per synonymous site ($K_s$)")
    axes[0].text(-0.11, 1.03, "a", transform=axes[0].transAxes, fontsize=10, fontweight="bold", va="bottom")
    axes[1].text(-0.11, 1.03, "b", transform=axes[1].transAxes, fontsize=10, fontweight="bold", va="bottom")
    axes[0].set_ylim(0, 8.0)
    axes[1].set_ylim(0, 4.5)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.97, bottom=0.10, hspace=0.52)

    prefix = args.outdir / "figure2_panel_a_reference_ks"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    write_tsv(args.outdir / "panel_a_density_grid.tsv", ["panel", "comparison", "Ks_grid", "density", "n", "bandwidth"], source_rows)
    write_tsv(args.outdir / "panel_a_peak_labels.tsv",
              ["panel", "comparison", "KDE_mode_Ks", "mode_density", "n"], peak_rows)
    summary_rows = []
    summary_grid = np.linspace(0.0, shared_xmax, 1001)
    for panel, grouped in (("upper", upper), ("lower", lower)):
        for comparison, values in grouped.items():
            density = reflected_gaussian_kde(values, summary_grid, args.bandwidth)
            mode_index = int(np.argmax(density))
            summary_rows.append({
                "panel": panel,
                "comparison": comparison,
                "n": len(values),
                "KDE_mode_Ks": "{:.8f}".format(float(summary_grid[mode_index])),
                "mode_density": "{:.8f}".format(float(density[mode_index])),
                "bandwidth": "{:.4f}".format(args.bandwidth),
                "xmax": "{:.3f}".format(shared_xmax),
            })
    write_tsv(args.outdir / "panel_a_plot_summary.tsv",
              ["panel", "comparison", "n", "KDE_mode_Ks", "mode_density", "bandwidth", "xmax"], summary_rows)
    (args.outdir / "plot_parameters.tsv").write_text(
        "backend\tPython matplotlib\n"
        "kde\tGaussian KDE with reflection at Ks=0\n"
        "bandwidth\t{:.4f}\n"
        "filter\t0<Ks<5 before density estimation\n"
        "median_lines\tnot_drawn; medians are descriptive statistics only\n"
        "shared_xmax\t{:.3f}\n"
        "upper_xmax_requested\t{:.3f}\n"
        "lower_xmax_requested\t{:.3f}\n".format(args.bandwidth, shared_xmax, args.upper_xmax, args.lower_xmax),
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Formal Figure 2a Ks density plots completed.\n", encoding="ascii")
    print("[OK] Figure 2a density plots written: {}".format(args.outdir))


if __name__ == "__main__":
    main()
