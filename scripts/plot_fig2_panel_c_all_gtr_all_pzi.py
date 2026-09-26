#!/usr/bin/env python3
"""Plot nine GTR retention panels with all 8+8 PZI chromosome curves."""

from __future__ import print_function

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.linewidth": 0.7, "axes.spines.right": False, "axes.spines.top": False,
    "legend.frameon": False,
})

PZI_TRACKS = ("PZI_A", "PZI_B")
OTHER_TRACKS = ("PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")
COLORS = {
    "PZI_A": "#B33C3C", "PZI_B": "#2E63A4",
    "PSP_A": "#C77C2A", "PSP_B": "#C77C2A",
    "PPA_A": "#3C8B58", "PPA_B": "#3C8B58",
    "PCAM": "#7B5AA6", "PPE": "#555555",
}


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def chrom_key(value):
    match = re.search(r"C(\d+)", value)
    if match:
        return (int(match.group(1)), value)
    return (10 ** 9, value)


def short(value):
    value = value.split("__", 1)[-1]
    match = re.search(r"C\d+_[AB](?:1|2)?", value)
    return match.group(0) if match else value


def draw_panel(ax, target, rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["track_id"], row["query_seqid"])].append(row)
    # Plot all PZI chromosomes first so they are never hidden by representative lines.
    ordered = sorted(groups, key=lambda key: (0 if key[0] == "PZI_A" else 1 if key[0] == "PZI_B" else 2, chrom_key(key[1]), key[0]))
    end_labels = []
    for track, query in ordered:
        values = sorted(groups[(track, query)], key=lambda row: float(row["window_center_gene"]))
        x = np.asarray([float(row["window_center_gene"]) for row in values])
        y = np.asarray([float(row["retention_percent"]) for row in values])
        if track == "PZI_A":
            alpha, lw, ls, z = 0.75, 0.70, "-", 4
        elif track == "PZI_B":
            alpha, lw, ls, z = 0.75, 0.70, "--", 4
        else:
            alpha, lw, ls, z = 0.65, 0.90, ":" if track == "PPE" else "-", 2
        ax.plot(x, y, color=COLORS[track], alpha=alpha, lw=lw, ls=ls, zorder=z)
        if track in PZI_TRACKS:
            end_labels.append((x[-1], y[-1], short(query), COLORS[track], ls))
    # Direct labels are restricted to PZI lines and placed in a compact right gutter.
    end_labels.sort(key=lambda item: item[1])
    adjusted = []
    for x_end, y_end, label, color, ls in end_labels:
        pos = y_end
        if adjusted:
            pos = max(pos, adjusted[-1][1] + 4.0)
        adjusted.append((label, pos, x_end, y_end, color, ls))
    if adjusted and adjusted[-1][1] > 100:
        overflow = adjusted[-1][1] - 100
        adjusted = [(label, pos-overflow, x_end, y_end, color, ls) for label, pos, x_end, y_end, color, ls in adjusted]
    for label, pos, x_end, y_end, color, ls in adjusted:
        ax.plot([x_end, x_end * 1.015], [y_end, pos], color=color, lw=0.3, clip_on=False)
        ax.text(x_end * 1.02, pos, label, color=color, fontsize=3.5, va="center", ha="left", clip_on=False)
    n_windows = max((len(v) for v in groups.values()), default=1)
    max_x = max(float(r["window_center_gene"]) for r in rows)
    ax.set_xlim(min(float(r["window_center_gene"]) for r in rows), max_x * 1.28)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.grid(axis="y", color="#D9D9D9", lw=0.35, alpha=0.7)
    ax.tick_params(length=2.0, width=0.55, pad=1.2, labelsize=5.3)
    ax.set_title(target.split("__", 1)[-1], loc="right", fontsize=6.5, color="#444444", pad=2)
    return n_windows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", required=True, type=Path)
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    rows = read_tsv(args.windows)
    selection = read_tsv(args.selection)
    if not rows:
        raise SystemExit("ERROR: empty all-PZI window table")
    if any(row.get("method") != "FraBias_local_IO_port_all_gtr_all_PZI" for row in rows):
        raise SystemExit("ERROR: windows are not from the all-PZI formal method")
    targets = sorted({r["target_seqid"] for r in rows}, key=chrom_key)
    if len(targets) != 9:
        raise SystemExit("ERROR: expected 9 GTR targets, found {}".format(len(targets)))
    for target in targets:
        present = {(r["track_id"], r["query_seqid"]) for r in rows if r["target_seqid"] == target}
        for track in PZI_TRACKS:
            if len([q for t, q in present if t == track]) != 8:
                raise SystemExit("ERROR: {} does not contain eight {} chromosome curves".format(target, track))

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 3, figsize=(12.2, 8.3), sharey=True)
    for index, (ax, target) in enumerate(zip(axes.flat, targets)):
        target_rows = [r for r in rows if r["target_seqid"] == target]
        draw_panel(ax, target, target_rows)
        if index % 3 == 0:
            ax.set_ylabel("Anchor retention (%)", fontsize=6.1)
        if index >= 6:
            ax.set_xlabel("GTR gene order", fontsize=6.1)
    handles = [
        plt.Line2D([], [], color=COLORS["PZI_A"], lw=1.0, ls="-", label="PZI-A (C01-C08)"),
        plt.Line2D([], [], color=COLORS["PZI_B"], lw=1.0, ls="--", label="PZI-B (C01-C08)"),
        plt.Line2D([], [], color=COLORS["PSP_A"], lw=1.0, label="PSP A/B representative"),
        plt.Line2D([], [], color=COLORS["PPA_A"], lw=1.0, label="PPA A/B representative"),
        plt.Line2D([], [], color=COLORS["PCAM"], lw=1.0, label="PCAM representative"),
        plt.Line2D([], [], color=COLORS["PPE"], lw=1.0, ls=":", label="PPE representative"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=5.7, frameon=False, bbox_to_anchor=(0.52, 0.005))
    fig.text(0.5, 0.035, "PZI-A and PZI-B: all eight chromosome curves retained in every GTR panel; other taxa: maximum-anchor representative", ha="center", fontsize=5.7, color="#444444")
    fig.text(0.012, 0.985, "c", fontsize=12, fontweight="bold", va="top")
    fig.subplots_adjust(left=0.06, right=0.985, top=0.96, bottom=0.10, hspace=0.32, wspace=0.25)
    prefix = outdir / "figure2_panel_c_all_gtr_all_pzi_chromosomes"
    fig.savefig(str(prefix)+".svg", bbox_inches="tight")
    fig.savefig(str(prefix)+".pdf", bbox_inches="tight")
    fig.savefig(str(prefix)+".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix)+".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Also emit one publication-ready file per GTR chromosome, preserving the nine plots separately.
    for target in targets:
        single = plt.figure(figsize=(5.6, 3.8))
        ax = single.add_subplot(111)
        draw_panel(ax, target, [r for r in rows if r["target_seqid"] == target])
        ax.set_xlabel("GTR gene order", fontsize=7)
        ax.set_ylabel("Anchor retention (%)", fontsize=7)
        single.subplots_adjust(left=0.11, right=0.83, top=0.91, bottom=0.16)
        one = outdir / ("figure2_panel_c_" + target.split("__", 1)[-1].lower() + "_all_pzi")
        single.savefig(str(one)+".pdf", bbox_inches="tight")
        single.savefig(str(one)+".png", dpi=300, bbox_inches="tight")
        plt.close(single)

    with (outdir / "all_gtr_all_pzi_source_data.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = list(rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (outdir / "all_gtr_all_pzi_selected_tracks.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = list(selection[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(selection)
    (outdir / "README.txt").write_text(
        "Nine GTR chromosome panels use formal quota-filtered anchors and overlapping 100-GTR-gene windows.\n"
        "Every panel includes PZI-A C01-C08 and PZI-B C01-C08; no-anchor chromosomes are retained as zero lines.\n"
        "PSP, PPA, PCAM and PPE each use their maximum-anchor query chromosome for that GTR target.\n"
        "The y-axis is quota-filtered anchor retention, not a validated biological gene-loss rate.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("Nine all-PZI GTR retention plots completed.\n", encoding="ascii")
    print("[OK] Wrote combined and nine individual all-PZI GTR retention plots: {}".format(outdir))


if __name__ == "__main__":
    main()
