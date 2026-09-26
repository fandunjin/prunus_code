#!/usr/bin/env python3
"""Plot the formal eight-track FraBias panel in reference Figure 3c style."""

from __future__ import print_function

import argparse
import csv
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
PAIR_ORDER = ("chr01+chr13", "chr02+chr15", "chr03+chr12")
STYLE = {
    "PZI_A": ("#C94C4C", "-"),
    "PZI_B": ("#C94C4C", "--"),
    "PSP_A": ("#3E72B8", "-"),
    "PSP_B": ("#3E72B8", "--"),
    "PPA_A": ("#3D8B5B", "-"),
    "PPA_B": ("#3D8B5B", "--"),
    "PCAM": ("#8B6AA8", "-"),
    "PPE": ("#555555", ":"),
}


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def short_seqid(seqid):
    return seqid.split("__", 1)[-1]


def spread_labels(items, lower=5.0, upper=95.0, gap=5.0):
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
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    if "jcvi" in str(args.fractbias).lower() or "alltaxa_fractionation_states" in str(args.fractbias):
        raise SystemExit("ERROR: JCVI/candidate results cannot be plotted as formal FraBias")
    rows = read_tsv(args.fractbias)
    required = {
        "pair_label", "target_seqid", "window_center_gene", "track_id",
        "query_seqid", "retention_percent", "method",
    }
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("ERROR: FraBias export lacks required columns")
    if not all(row["method"] == "FraBias_local_IO_port" for row in rows):
        raise SystemExit("ERROR: mixed or unregistered method rows in formal FraBias source data")
    rows = [row for row in rows if row["pair_label"] in PAIR_ORDER and row["track_id"] in TRACKS]
    counts = defaultdict(set)
    for row in rows:
        value = float(row["retention_percent"])
        if not 0.0 <= value <= 100.0:
            raise SystemExit("ERROR: retention_percent outside [0,100]")
        counts[row["pair_label"]].add(row["track_id"])
    for pair in PAIR_ORDER:
        if counts[pair] != set(TRACKS):
            raise SystemExit("ERROR: {} does not contain all eight FraBias tracks".format(pair))

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(3.55, 3.9), squeeze=False)
    axes = axes[:, 0]
    for panel_index, (ax, pair) in enumerate(zip(axes, PAIR_ORDER)):
        selected = [row for row in rows if row["pair_label"] == pair]
        target_seqids = sorted(set(row["target_seqid"] for row in selected))
        if len(target_seqids) != 1:
            raise SystemExit("ERROR: {} maps to {} target chromosomes".format(pair, len(target_seqids)))
        endpoint_items = []
        series = {}
        for track in TRACKS:
            track_rows = sorted(
                (row for row in selected if row["track_id"] == track),
                key=lambda row: float(row["window_center_gene"]),
            )
            x = np.asarray([float(row["window_center_gene"]) for row in track_rows])
            y = np.asarray([float(row["retention_percent"]) for row in track_rows])
            color, linestyle = STYLE[track]
            ax.plot(x, y, color=color, linestyle=linestyle, lw=0.85, solid_capstyle="round")
            label = track.replace("_", "-") + " " + short_seqid(track_rows[0]["query_seqid"])
            endpoint_items.append((label, float(y[-1])))
            series[label] = (float(x[-1]), float(y[-1]), color)

        label_y = spread_labels(endpoint_items)
        for label, (x_end, y_end, color) in series.items():
            ax.plot([x_end, x_end * 1.025], [y_end, label_y[label]], color=color, lw=0.35, clip_on=False)
            ax.text(x_end * 1.032, label_y[label], label, color=color, fontsize=4.45,
                    ha="left", va="center", clip_on=False)

        x_values = [float(row["window_center_gene"]) for row in selected]
        ax.set_xlim(min(x_values), max(x_values) * 1.29)
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 50, 100])
        ax.tick_params(length=2.2, width=0.65, pad=1.5)
        ax.grid(axis="y", color="#D9D9D9", lw=0.35, alpha=0.65)
        ax.text(0.01, 0.08, pair.replace("+", " / "), transform=ax.transAxes,
                fontweight="bold", fontsize=6.4)
        ax.text(0.99, 1.02, short_seqid(target_seqids[0]), transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5.2, color="#555555")
        if panel_index < 2:
            ax.tick_params(labelbottom=False)

    axes[1].set_ylabel("Gene retention (%)")
    axes[-1].set_xlabel("GTR target chromosome gene order")
    fig.text(0.012, 0.985, "c", fontsize=9, fontweight="bold", va="top")
    fig.subplots_adjust(left=0.15, right=0.76, top=0.97, bottom=0.11, hspace=0.16)

    prefix = outdir / "figure2_panel_c_fractbias_reference_exact"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    write_tsv(outdir / "panel_c_source_data.tsv", list(rows[0].keys()), rows)
    (outdir / "COMPLETE").write_text("Formal FraBias panel plotted.\n", encoding="ascii")
    print("[OK] Wrote formal panel c: {}".format(outdir))


if __name__ == "__main__":
    main()
