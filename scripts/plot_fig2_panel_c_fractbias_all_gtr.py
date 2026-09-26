#!/usr/bin/env python3
"""Plot an all-nine-GTR-chromosome audit grid from formal FraBias windows."""

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
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})

TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")
STYLE = {
    "PZI_A": ("#C94C4C", "-", "PZI-A"),
    "PZI_B": ("#C94C4C", "--", "PZI-B"),
    "PSP_A": ("#3E72B8", "-", "PSP-A"),
    "PSP_B": ("#3E72B8", "--", "PSP-B"),
    "PPA_A": ("#3D8B5B", "-", "PPA-A"),
    "PPA_B": ("#3D8B5B", "--", "PPA-B"),
    "PCAM": ("#8B6AA8", "-", "PCAM"),
    "PPE": ("#555555", ":", "PPE"),
}


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def chromosome_key(seqid):
    match = re.search(r"(\d+)$", seqid)
    return (int(match.group(1)) if match else 10 ** 9, seqid)


def short_seqid(seqid):
    return seqid.split("__", 1)[-1]


def spread_labels(items, lower=4.0, upper=96.0, gap=7.5):
    ordered = sorted(items, key=lambda item: item[1])
    positions = []
    for _, value in ordered:
        positions.append(max(value, positions[-1] + gap if positions else lower))
    if positions and positions[-1] > upper:
        shift = positions[-1] - upper
        positions = [value - shift for value in positions]
        for index in range(len(positions) - 2, -1, -1):
            positions[index] = min(positions[index], positions[index + 1] - gap)
        if positions[0] < lower:
            positions = list(np.linspace(lower, upper, len(positions)))
    return {label: position for (label, _), position in zip(ordered, positions)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fractbias", required=True)
    parser.add_argument("--selection", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    rows = read_tsv(args.fractbias)
    selection = read_tsv(args.selection)
    expected_method = "FraBias_local_IO_port_all_gtr_best_track"
    if not rows or any(row.get("method") != expected_method for row in rows):
        raise SystemExit("ERROR: source rows are not all-GTR formal FraBias windows")
    targets = sorted({row["target_seqid"] for row in rows}, key=chromosome_key)
    if len(targets) != 9:
        raise SystemExit("ERROR: expected 9 chromosome-scale GTR targets, found {}".format(len(targets)))
    for target in targets:
        found = {row["track_id"] for row in rows if row["target_seqid"] == target}
        if found != set(TRACKS):
            raise SystemExit("ERROR: {} does not have all eight tracks".format(target))

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 3, figsize=(9.1, 6.65), sharey=True)
    for index, (ax, target) in enumerate(zip(axes.flat, targets)):
        chosen = [row for row in rows if row["target_seqid"] == target]
        selected_for_target = {
            row["track_id"]: row["query_seqid"]
            for row in selection if row["target_seqid"] == target
        }
        unavailable = []
        endpoint_items = []
        endpoint_series = {}
        for track in TRACKS:
            track_rows = sorted(
                (row for row in chosen if row["track_id"] == track),
                key=lambda row: float(row["window_center_gene"]),
            )
            if track_rows[0].get("status") == "no_quota_filtered_anchor":
                unavailable.append(STYLE[track][2])
                continue
            x = np.asarray([float(row["window_center_gene"]) for row in track_rows])
            y = np.asarray([float(row["retention_percent"]) for row in track_rows])
            color, linestyle, label = STYLE[track]
            line, = ax.plot(x, y, color=color, linestyle=linestyle, lw=0.75, label=label)
            direct_label = "{} {}".format(label, short_seqid(selected_for_target[track]))
            endpoint_items.append((direct_label, float(y[-1])))
            endpoint_series[direct_label] = (float(x[-1]), float(y[-1]), color)
        label_y = spread_labels(endpoint_items)
        for label, (x_end, y_end, color) in endpoint_series.items():
            ax.plot([x_end, x_end * 1.025], [y_end, label_y[label]],
                    color=color, lw=0.35, clip_on=False)
            ax.text(x_end * 1.035, label_y[label], label, color=color,
                    fontsize=4.0, ha="left", va="center", clip_on=False)
        if unavailable:
            ax.text(0.02, 0.05, "No QA anchor: " + ", ".join(unavailable),
                    transform=ax.transAxes, fontsize=5.0, color="#6B6B6B", va="bottom")
        x_values = [float(row["window_center_gene"]) for row in chosen]
        ax.set_xlim(min(x_values), max(x_values) * 1.42)
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 50, 100])
        ax.grid(axis="y", color="#D9D9D9", lw=0.35, alpha=0.7)
        ax.tick_params(length=2.0, width=0.6, pad=1.2, labelsize=5.5)
        ax.set_title(target.split("__", 1)[-1], loc="right", fontsize=6.5, color="#555555", pad=2)
        if index % 3 == 0:
            ax.set_ylabel("Anchor retention (%)", fontsize=6.1)
        if index >= 6:
            ax.set_xlabel("GTR gene order", fontsize=6.1)

    fig.subplots_adjust(left=0.07, right=0.995, top=0.97, bottom=0.08, hspace=0.34, wspace=0.38)
    prefix = outdir / "figure2_panel_c_all_gtr_chromosome_audit"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    write_tsv(outdir / "all_gtr_selected_tracks.tsv", list(selection[0].keys()), selection)
    write_tsv(outdir / "all_gtr_panel_source_data.tsv", list(rows[0].keys()), rows)
    (outdir / "README.txt").write_text(
        "Each subplot is one GTR chromosome. Each taxon/subgenome line uses the physical query chromosome with the maximum number of formal quota-filtered anchors to that GTR chromosome.\n"
        "The y-axis is collinear-anchor retention in an overlapping 100-GTR-gene window; it is not a validated gene-loss rate.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("All-nine-GTR audit grid plotted.\n", encoding="ascii")
    print("[OK] Wrote all-nine-GTR audit grid: {}".format(outdir))


if __name__ == "__main__":
    main()
