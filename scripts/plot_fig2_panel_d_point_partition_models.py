#!/usr/bin/env python3
"""Plot a compact, annotation-free comparison of the five POInT models."""

from __future__ import print_function

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_rows(path):
    with Path(path).open("r", encoding="ascii", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    rows = read_rows(args.summary)
    if len(rows) != 5:
        raise SystemExit("ERROR: expected five model rows")
    # Compact labels identify only the partition; no explanatory text is added.
    labels = ["ALL", "12|3", "13|2", "23|1", "IND"]
    delta_aic = np.array([float(row["delta_AIC"]) for row in rows])
    delta_bic = np.array([float(row["delta_BIC"]) for row in rows])
    lnls = np.array([float(row["log_likelihood"]) for row in rows])
    best = int(np.argmax(lnls))
    colors = ["#4C78A8", "#72B7B2", "#F2CF5B", "#E08E79", "#B0B0B0"]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.linewidth": 0.75,
        "axes.spines.right": False,
        "axes.spines.top": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4), sharex=True, gridspec_kw={"wspace": 0.34})
    x = np.arange(len(rows))
    for ax, values, ylabel in zip(axes, (delta_aic, delta_bic), ("Delta AIC", "Delta BIC")):
        ax.bar(x, values, width=0.66, color=colors, edgecolor="none")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        ax.tick_params(length=3, width=0.7)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_linewidth(0.75)
    # Keep the best model visually distinct without a prose annotation.
    axes[0].bar(best, delta_aic[best], width=0.66, color="#1F4E79", edgecolor="none")
    axes[1].bar(best, delta_bic[best], width=0.66, color="#1F4E79", edgecolor="none")
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.22, top=0.96)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / "figure2_panel_d_point_partition_models"
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)
    (outdir / "COMPLETE").write_text("Five-model POInT comparison plotted.\n", encoding="ascii")
    print("[OK] Wrote five-model POInT panel: {}".format(outdir))


if __name__ == "__main__":
    main()
