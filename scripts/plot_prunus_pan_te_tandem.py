#!/usr/bin/env python3
"""Render publication-ready tandem-repeat composition and landscape figures.

The figures intentionally show individual haplotypes plus group ranges instead
of treating phased haplotypes as independent biological replicates.
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2", "pspA1", "pspA2", "pspB1", "pspB2", "ppaA", "ppaB", "pcam", "ppe"]
GROUPS = ["PZI", "PSP", "PPA", "PCAM", "PPE"]
PERIOD_CLASSES = ["1", "2-6", "7-20", "21-100", "101-500", "501-2000", ">2000"]
PERIOD_LABELS = {"1": "1", "2-6": "2-6", "7-20": "7-20", "21-100": "21-100", "101-500": "101-500", "501-2000": "501-2000", ">2000": ">2000"}
COLORS = {"1": "#CFCFCF", "2-6": "#99B9D6", "7-20": "#5B7FCA", "21-100": "#42949E", "101-500": "#D99042", "501-2000": "#B64342", ">2000": "#765A9D"}
GROUP_COLORS = {"PZI": "#B64342", "PSP": "#42949E", "PPA": "#5B7FCA", "PCAM": "#767676", "PPE": "#9A4D8E"}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 7, "axes.linewidth": 0.75, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.major.width": 0.75, "ytick.major.width": 0.75,
    "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
})


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(fig, output_dir, stem):
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix, kwargs in (("svg", {}), ("pdf", {}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}}), ("png", {"dpi": 300})):
        path = output_dir / f"{stem}.{suffix}"
        fig.savefig(path, bbox_inches="tight", facecolor="white", **kwargs)
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError(f"failed to export {path}")
    plt.close(fig)


def label(ax, value):
    ax.text(-0.09, 1.04, value, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def group_boundary(ax, y):
    ax.axhline(y, color="#D9D9D9", lw=0.6, zorder=0)


def plot_summary(genome_rows, period_rows, motif_rows, output_dir, top_motifs):
    genome = {r["target"]: r for r in genome_rows}
    period = defaultdict(dict)
    for row in period_rows:
        period[row["target"]][row["period_class"]] = int(row["exclusive_union_bp"]) / 1e6
    fig = plt.figure(figsize=(7.25, 6.2))
    grid = fig.add_gridspec(2, 2, width_ratios=(1.27, 1), height_ratios=(1, 1.04), wspace=0.5, hspace=0.55)

    ax = fig.add_subplot(grid[0, 0])
    y = np.arange(len(TARGETS))
    left = np.zeros(len(TARGETS))
    active_classes = [cls for cls in PERIOD_CLASSES if any(period[target].get(cls, 0) > 0 for target in TARGETS)]
    for cls in active_classes:
        vals = np.array([period[target].get(cls, 0) for target in TARGETS])
        ax.barh(y, vals, left=left, height=0.68, color=COLORS[cls], lw=0, label=PERIOD_LABELS[cls])
        left += vals
    for i, target in enumerate(TARGETS):
        ax.text(left[i] + max(left) * 0.012, i, f"{left[i]:.1f}", va="center", fontsize=5.8)
    ax.set_yticks(y, TARGETS)
    ax.invert_yaxis(); ax.set_xlabel("Non-redundant tandem repeat (Mb)")
    ax.set_xlim(0, max(left) * 1.13 if max(left) else 1)
    ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    for boundary in (3.5, 7.5, 9.5, 10.5): group_boundary(ax, boundary)
    ax.legend(title="TRF period (bp)", title_fontsize=6, ncol=2, fontsize=5.5, columnspacing=0.8, handlelength=1.0, bbox_to_anchor=(0, 1.02), loc="lower left")
    label(ax, "a")

    ax = fig.add_subplot(grid[0, 1])
    offsets = {"PZI": [-0.15, -0.05, 0.05, 0.15], "PSP": [-0.15, -0.05, 0.05, 0.15], "PPA": [-0.06, 0.06], "PCAM": [0], "PPE": [0]}
    for x, group in enumerate(GROUPS):
        rows = [r for r in genome_rows if r["group"] == group]
        values = np.array([float(r["genome_fraction"]) * 100 for r in rows])
        ax.plot([x, x], [values.min(), values.max()], color=GROUP_COLORS[group], lw=1.2, zorder=1)
        ax.scatter([x] * len(values) + np.array(offsets[group]), values, s=15, color=GROUP_COLORS[group], edgecolor="white", linewidth=0.35, zorder=2)
        ax.scatter(x, values.mean(), marker="_", s=150, color="#202020", linewidth=1.2, zorder=3)
    ax.set_xticks(range(len(GROUPS)), GROUPS)
    ax.set_ylabel("Tandem repeat fraction (%)")
    ax.set_xlim(-0.5, len(GROUPS)-0.5); ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    ax.text(0.02, 1.03, "Points: haplotypes; horizontal mark: group mean", transform=ax.transAxes, fontsize=5.5, color="#555555")
    label(ax, "b")

    ax = fig.add_subplot(grid[1, 0])
    chosen = sorted(motif_rows, key=lambda r: int(r["pzi_minus_psp_mean_bp"]), reverse=True)[:top_motifs]
    chosen = list(reversed(chosen))
    matrix = np.array([[int(row.get(f"{group.lower()}_mean_bp", 0)) / 1e6 for group in GROUPS] for row in chosen]) if chosen else np.zeros((1, len(GROUPS)))
    image = ax.imshow(matrix, aspect="auto", cmap="YlOrRd")
    labels = []
    for row in chosen:
        motif = row["canonical_motif"]
        motif = motif if len(motif) <= 18 else motif[:15] + "..."
        labels.append(f"{row['period_bp']} bp  {motif}")
    ax.set_xticks(range(len(GROUPS)), GROUPS)
    ax.set_yticks(range(len(labels)), labels, fontsize=5.5)
    cbar = fig.colorbar(image, ax=ax, fraction=0.05, pad=0.035)
    cbar.set_label("Mean non-redundant span (Mb)", fontsize=6)
    cbar.ax.tick_params(labelsize=5.5)
    ax.set_title("Motifs with largest PZI - PSP difference", loc="left", fontsize=7)
    label(ax, "c")

    ax = fig.add_subplot(grid[1, 1])
    for group in GROUPS:
        rows = [r for r in genome_rows if r["group"] == group]
        xs = np.array([int(r["genome_bp"]) / 1e6 for r in rows])
        ys = np.array([int(r["exclusive_union_bp"]) / 1e6 for r in rows])
        ax.scatter(xs, ys, s=24, color=GROUP_COLORS[group], label=group, edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Haplotype/genome size (Mb)")
    ax.set_ylabel("Tandem repeat span (Mb)")
    ax.grid(color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    ax.legend(fontsize=5.5, ncol=2, loc="upper left", handletextpad=0.35, columnspacing=0.8)
    ax.text(0.02, 0.02, "Descriptive comparison; haplotypes are not independent replicates", transform=ax.transAxes, fontsize=5.2, color="#555555", va="bottom")
    label(ax, "d")
    save(fig, output_dir, "supp_tandem_repeat_summary")


def normalized_matrix(rows, target, n_bins=100):
    selected = [r for r in rows if r["target"] == target]
    by_seq = defaultdict(list)
    for row in selected:
        by_seq[row["seqid"]].append(row)
    chroms = sorted(by_seq, key=lambda x: (-max(int(r["end"]) for r in by_seq[x]), x))[:8]
    matrix = np.zeros((len(chroms), n_bins), dtype=float)
    names = []
    for i, seqid in enumerate(chroms):
        seqrows = by_seq[seqid]
        seqend = max(int(r["end"]) for r in seqrows)
        names.append(seqid.replace("chr", "").replace("_A1", "").replace("_A2", "").replace("_B1", "").replace("_B2", ""))
        for row in seqrows:
            a, b = int(row["start"]), int(row["end"])
            value = int(row["tandem_union_bp"])
            if not value or b <= a:
                continue
            density = value / (b - a)
            lo, hi = max(0, int(a / seqend * n_bins)), min(n_bins - 1, int((b - 1) / seqend * n_bins))
            for j in range(lo, hi + 1):
                left, right = max(a, j * seqend / n_bins), min(b, (j + 1) * seqend / n_bins)
                matrix[i, j] += density * max(0, right - left) / (seqend / n_bins)
    return matrix, names


def plot_pzi_landscape(window_rows, output_dir):
    targets = ["pziA1", "pziA2", "pziB1", "pziB2"]
    matrices, names = zip(*(normalized_matrix(window_rows, target) for target in targets))
    vmax = max(float(matrix.max()) for matrix in matrices) if matrices else 1
    vmax = max(vmax, 0.01)
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 4.45), sharex=True, sharey=True)
    image = None
    for index, (ax, target, matrix, rows) in enumerate(zip(axes.flat, targets, matrices, names)):
        image = ax.imshow(matrix * 100, aspect="auto", cmap="magma", vmin=0, vmax=vmax * 100, interpolation="nearest")
        ax.set_title(target.replace("pzi", "PZI-"), loc="left", fontsize=7)
        ax.set_yticks(range(len(rows)), rows, fontsize=5.7)
        ax.set_xticks([0, 25, 50, 75, 99], ["0", "25", "50", "75", "100"])
        if index >= 2: ax.set_xlabel("Normalized chromosome position (%)")
        if index % 2 == 0: ax.set_ylabel("Chromosome")
        label(ax, chr(ord("a") + index))
    cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.022, pad=0.025)
    cbar.set_label("Tandem repeat density (%)", fontsize=6)
    cbar.ax.tick_params(labelsize=5.5)
    fig.subplots_adjust(wspace=0.28, hspace=0.4, right=0.89)
    save(fig, output_dir, "supp_tandem_repeat_pzi_chromosome_density")


def plot_repeatmasker_overlap(rows, output_dir):
    values = {r["target"]: r for r in rows}
    y = np.arange(len(TARGETS))
    shared = np.array([int(values[t]["repeatmasker_overlap_bp"]) / 1e6 for t in TARGETS])
    unique = np.array([int(values[t]["tandem_only_bp"]) / 1e6 for t in TARGETS])
    fig, ax = plt.subplots(figsize=(5.8, 4.1))
    ax.barh(y, shared, height=0.68, color="#8F8F8F", label="Covered by pan-TE RepeatMasker library")
    ax.barh(y, unique, left=shared, height=0.68, color="#D8D8D8", label="Not covered by library")
    ax.set_yticks(y, TARGETS); ax.invert_yaxis(); ax.set_xlabel("Non-redundant TRF span (Mb)")
    ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    for boundary in (3.5, 7.5, 9.5, 10.5): group_boundary(ax, boundary)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), fontsize=5.8, ncol=1)
    ax.text(0, -0.18, "Annotation-overlap diagnostic; coverage is not independent evidence of transposition.", transform=ax.transAxes, fontsize=5.5, color="#555555")
    save(fig, output_dir, "supp_tandem_repeat_repeatmasker_overlap")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--top-motifs", type=int, default=12)
    args = p.parse_args()
    needed = ["tandem_genome_summary.tsv", "tandem_period_summary.tsv", "tandem_group_motif_summary.tsv", "tandem_window_density.tsv"]
    missing = [name for name in needed if not (args.summary_dir / name).is_file()]
    if missing:
        raise SystemExit("ERROR: missing tandem summaries: " + ", ".join(missing))
    plot_summary(read_tsv(args.summary_dir / "tandem_genome_summary.tsv"), read_tsv(args.summary_dir / "tandem_period_summary.tsv"), read_tsv(args.summary_dir / "tandem_group_motif_summary.tsv"), args.output_dir, args.top_motifs)
    plot_pzi_landscape(read_tsv(args.summary_dir / "tandem_window_density.tsv"), args.output_dir)
    overlap = args.summary_dir / "tandem_repeatmasker_overlap.tsv"
    if overlap.is_file():
        plot_repeatmasker_overlap(read_tsv(overlap), args.output_dir)
    print(f"Figures written to {args.output_dir}")


if __name__ == "__main__":
    main()
