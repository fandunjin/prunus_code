#!/usr/bin/env python3
"""Export Figure-1 tandem-repeat panels and their sensitivity supplements."""

import argparse
import collections
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2", "pspA1", "pspA2", "pspB1", "pspB2", "ppaA", "ppaB", "pcam", "ppe"]
PZI_TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2"]
GROUPS = ["PZI", "PSP", "PPA", "PCAM", "PPE"]
PARTITION_COLORS = {"trf_tandem_bp": "#BF4A45", "pan_te_non_tandem_bp": "#5476B8", "other_bp": "#D5D5D5"}
PERIOD_COLORS = {"1": "#CFCFCF", "2-6": "#99B9D6", "7-20": "#5B7FCA", "21-100": "#42949E", "101-500": "#D99042", "501-2000": "#B64342", ">2000": "#765A9D"}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7, "axes.linewidth": 0.75, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.75, "ytick.major.width": 0.75, "legend.frameon": False,
    "svg.fonttype": "none", "pdf.fonttype": 42, "figure.facecolor": "white", "axes.facecolor": "white",
})


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


def plot_main(partition, periods, loci, ont, output_dir):
    by_target = {row["target"]: row for row in partition}
    fig = plt.figure(figsize=(7.25, 6.25))
    grid = fig.add_gridspec(2, 2, width_ratios=(1.22, 1.0), height_ratios=(1.05, 0.95), wspace=0.42, hspace=0.54)

    ax = fig.add_subplot(grid[0, 0])
    y = np.arange(len(TARGETS)); left = np.zeros(len(TARGETS))
    for field, name in (("trf_tandem_bp", "TRF tandem arrays"), ("pan_te_non_tandem_bp", "pan-TE outside TRF"), ("other_bp", "Other sequence")):
        values = np.array([int(by_target[target][field]) / 1e6 for target in TARGETS])
        ax.barh(y, values, left=left, height=0.68, color=PARTITION_COLORS[field], label=name, linewidth=0)
        left += values
    for index, value in enumerate(left):
        ax.text(value + 3, index, f"{value:.0f}", va="center", fontsize=5.7)
    ax.set_yticks(y, [target.replace("pzi", "PZI-").replace("psp", "PSP-").replace("ppa", "PPA-").upper() if target not in ("pcam", "ppe") else target.upper() for target in TARGETS])
    ax.invert_yaxis(); ax.set_xlabel("Genome package size (Mb)")
    ax.set_xlim(0, max(left) * 1.15); ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    for boundary in (3.5, 7.5, 9.5, 10.5): ax.axhline(boundary, color="#D9D9D9", lw=0.6)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.03), fontsize=5.5, ncol=1, handlelength=1.2)
    label(ax, "a")

    ax = fig.add_subplot(grid[0, 1])
    pairs = [("PZI-A1 - PSP-A1", "pziA1", "pspA1"), ("PZI-A2 - PSP-A2", "pziA2", "pspA2"), ("PZI-B1 - PSP-B1", "pziB1", "pspB1"), ("PZI-B2 - PSP-B2", "pziB2", "pspB2")]
    x = np.arange(len(pairs)); bottom = np.zeros(len(pairs))
    for field, name in (("trf_tandem_bp", "TRF tandem arrays"), ("pan_te_non_tandem_bp", "pan-TE outside TRF"), ("other_bp", "Other sequence")):
        values = np.array([(int(by_target[pzi][field]) - int(by_target[psp][field])) / 1e6 for _, pzi, psp in pairs])
        ax.bar(x, values, bottom=bottom, width=0.62, color=PARTITION_COLORS[field], label=name, linewidth=0)
        bottom += values
    for index, value in enumerate(bottom):
        ax.text(index, value + max(bottom) * 0.025, f"{value:.0f}", ha="center", va="bottom", fontsize=5.8)
    ax.set_xticks(x, ["A1", "A2", "B1", "B2"]); ax.set_xlabel("Matched genomic packages")
    ax.set_ylabel("PZI - PSP sequence (Mb)"); ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    label(ax, "b")

    ax = fig.add_subplot(grid[1, 0])
    by_period = collections.defaultdict(dict)
    for row in periods:
        by_period[row["target"]][row["period_class"]] = int(row["exclusive_union_bp"]) / 1e6
    group_values = {}
    active = [name for name in PERIOD_COLORS if any(by_period[target].get(name, 0) > 0 for target in TARGETS)]
    for group in GROUPS:
        members = [target for target in TARGETS if by_target[target]["group"] == group]
        group_values[group] = {period: np.mean([by_period[target].get(period, 0) for target in members]) for period in active}
    x = np.arange(len(GROUPS)); bottom = np.zeros(len(GROUPS))
    for period in active:
        values = np.array([group_values[group][period] for group in GROUPS])
        ax.bar(x, values, bottom=bottom, width=0.65, color=PERIOD_COLORS[period], label=period, linewidth=0)
        bottom += values
    ax.set_xticks(x, GROUPS); ax.set_ylabel("Mean non-redundant tandem span (Mb)")
    ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    ax.legend(title="TRF period (bp)", title_fontsize=5.5, fontsize=5.0, ncol=2, columnspacing=0.8, handlelength=1.0, loc="lower left", bbox_to_anchor=(0, 1.01))
    ax.text(0.01, -0.24, "Means summarize genomic packages; they are not biological replicates.", transform=ax.transAxes, fontsize=5.1, color="#555555")
    label(ax, "c")

    ax = fig.add_subplot(grid[1, 1])
    if loci:
        locus_by_query = collections.defaultdict(list)
        for row in loci:
            locus_by_query[row["query_target"]].append(row)
        qs = [target for target in PZI_TARGETS if target in locus_by_query]
        categories = ("PZI counterpart callable", "PSP counterpart callable", "PZI-supported / PSP-contracted")
        y = np.arange(len(qs)); height = 0.22
        colors = ("#7A7A7A", "#E39C45", "#BF4A45")
        for offset, category, color in zip((-height, 0, height), categories, colors):
            values = []
            for query in qs:
                items = locus_by_query[query]
                if category == categories[0]: value = sum(int(item["pzi_counterpart_callable"]) > 0 for item in items)
                elif category == categories[1]: value = sum(int(item["psp_counterpart_callable"]) > 0 for item in items)
                else: value = sum("PZI_homologous_array_with_PSP" in item["locus_conclusion"] for item in items)
                values.append(100 * value / len(items) if items else 0)
            ax.barh(y + offset, values, height=height, color=color, label=category)
        ax.set_yticks(y, [query.replace("pzi", "PZI-") for query in qs]); ax.set_xlabel("Candidate array loci (%)")
        ax.set_xlim(0, 100); ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
        ax.legend(fontsize=5.0, ncol=1, loc="lower left", bbox_to_anchor=(0, 1.01))
        if ont:
            passed, tested = sum(int(row["two_junctions_supported"]) for row in ont), len(ont)
            ax.text(0.01, -0.25, f"ONT: {passed}/{tested} selected loci with both array junctions supported.", transform=ax.transAxes, fontsize=5.1, color="#555555")
    else:
        ax.axis("off")
        ax.text(0.02, 0.78, "Locus-level validation pending", fontsize=8, transform=ax.transAxes)
        ax.text(0.02, 0.61, "Plot only after flank-pair and ONT QC pass.", fontsize=6.5, color="#555555", transform=ax.transAxes)
    label(ax, "d")
    save(fig, output_dir, "figure1_tandem_repeat_expansion")


def plot_sensitivity(rows, output_dir):
    order = sorted({row["trf_score"] for row in rows}, key=lambda x: float(x))
    fig, ax = plt.subplots(figsize=(5.0, 3.1))
    colors = {"PZI": "#BF4A45", "PSP": "#42949E", "PPA": "#5B7FCA", "PCAM": "#767676", "PPE": "#9A4D8E"}
    for group in GROUPS:
        subset = {row["trf_score"]: row for row in rows if row["group"] == group}
        x = [float(score) for score in order if score in subset]
        y = [int(subset[score]["mean_trf_tandem_bp"]) / 1e6 for score in order if score in subset]
        ax.plot(x, y, marker="o", markersize=3.2, lw=1.1, color=colors[group], label=group)
    ax.set_xlabel("TRF minimum score"); ax.set_ylabel("Mean tandem span (Mb)")
    ax.grid(color="#E7E7E7", lw=0.5); ax.set_axisbelow(True); ax.legend(ncol=2, fontsize=5.5)
    save(fig, output_dir, "supp_tandem_repeat_score_sensitivity")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main-dir", type=Path, required=True)
    parser.add_argument("--tandem-summary-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--locus-dir", type=Path)
    parser.add_argument("--ont-dir", type=Path)
    parser.add_argument("--sensitivity", type=Path)
    args = parser.parse_args()
    partition = read_tsv(args.main_dir / "tandem_main_partition.tsv")
    periods = read_tsv(args.tandem_summary_dir / "tandem_period_summary.tsv")
    loci, ont = [], []
    if args.locus_dir and (args.locus_dir / "pzi_tandem_locus_summary.tsv").is_file():
        loci = read_tsv(args.locus_dir / "pzi_tandem_locus_summary.tsv")
    if args.ont_dir and args.ont_dir.is_dir():
        for path in sorted(args.ont_dir.glob("*.ont_support.tsv")):
            ont.extend(read_tsv(path))
    plot_main(partition, periods, loci, ont, args.output_dir)
    if args.sensitivity and args.sensitivity.is_file():
        plot_sensitivity(read_tsv(args.sensitivity), args.output_dir)
    print(f"Figures written to {args.output_dir}")


if __name__ == "__main__":
    main()
