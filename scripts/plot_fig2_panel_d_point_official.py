#!/usr/bin/env python3
"""Plot the official POInT likelihood result for Fig. 2d in Python."""

from __future__ import print_function

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.75,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def posterior_state_counts(path):
    rows = read_tsv(path)
    counts = {}
    for row in rows:
        state = row.get("best_tracking", "")
        if state:
            counts[state] = counts.get(state, 0) + 1
    return counts, len(rows)


def display_branch(branch):
    return {
        "NONE": "PZI+PSP stem",
        "Root": "ancestral stem",
    }.get(branch, branch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit", required=True)
    parser.add_argument("--transition", default=None)
    parser.add_argument("--posterior-qc", required=True)
    parser.add_argument("--posterior-summary", default=None)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    for input_path in (args.fit, args.transition, args.posterior_qc, args.posterior_summary):
        if not input_path:
            continue
        lowered = str(input_path).lower()
        if any(token in lowered for token in ("jcvi", "dollo", "candidate", "alltaxa_fractionation_states")):
            raise SystemExit("ERROR: non-POInT proxy output cannot be plotted")
    fit_rows = read_tsv(args.fit)
    transition_rows = read_tsv(args.transition) if args.transition and Path(args.transition).is_file() else []
    qc_rows = {row["field"]: row["value"] for row in read_tsv(args.posterior_qc)}
    if not fit_rows:
        raise SystemExit("ERROR: official POInT model summary is empty")
    if any(row.get("analysis_status", "").lower() not in ("converged", "complete", "converged_posterior_only") for row in fit_rows):
        raise SystemExit("ERROR: model fit is not marked converged")
    fit_names = [row["model"] for row in fit_rows]
    lnls = []
    for row in fit_rows:
        try:
            lnls.append(float(row["log_likelihood"]))
        except ValueError:
            lnls.append(float("nan"))
    if any(np.isnan(value) for value in lnls):
        raise SystemExit("ERROR: POInT log likelihood could not be parsed from stdout")
    n_sites = int(float(fit_rows[0]["pillars"]))
    # The rooted three-taxon POInT tree has five optimized branches: three
    # terminal branches, one internal branch, and the root/stem branch from
    # the fixed Dupl state to the first split.
    for row, ln_l in zip(fit_rows, lnls):
        n_params = int(row["model_parameters"]) + 5
        row["AIC"] = "{:.8g}".format(2.0 * n_params - 2.0 * ln_l)
        row["BIC"] = "{:.8g}".format(n_params * np.log(n_sites) - 2.0 * ln_l)
    min_aic = min(float(row["AIC"]) for row in fit_rows)
    for row in fit_rows:
        row["delta_AIC"] = "{:.8g}".format(float(row["AIC"]) - min_aic)
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    branches = []
    transitions = []
    matrix = {}
    for row in transition_rows:
        branch = display_branch(row["branch"])
        transition = row["transition"]
        if branch not in branches:
            branches.append(branch)
        if transition not in transitions:
            transitions.append(transition)
        matrix[(branch, transition)] = float(row["mean_probability"])

    fig = plt.figure(figsize=(7.0, 3.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.92, 1.55], wspace=0.35)
    ax_fit = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    ax_fit.set_title("model comparison", fontsize=7, pad=4)
    x = np.arange(len(fit_rows))
    deltas = [float(row["delta_AIC"]) for row in fit_rows]
    colors = ["#5B7FCA" if value == min(deltas) else "#BDBDBD" for value in deltas]
    ax_fit.bar(x, deltas, color=colors, width=0.65)
    ax_fit.set_xticks(x)
    display_names = {
        "symmetric": "symmetric",
        "reference_WGD_bias_fix": "WGD bias fix",
        "source_validated_copy_specific": "copy-specific",
    }
    ax_fit.set_xticklabels([display_names.get(name, name) for name in fit_names], rotation=25, ha="right")
    ax_fit.set_ylabel("Delta AIC")
    ax_fit.set_ylim(bottom=0)
    for idx, row in enumerate(fit_rows):
        ax_fit.text(idx, deltas[idx] + max(deltas + [1.0]) * 0.04, "lnL {}".format(row["log_likelihood"]), ha="center", va="bottom", fontsize=5.2, rotation=90)
    ax_fit.text(0.0, -0.29, "n={}; 5 POInT branches included".format(n_sites), transform=ax_fit.transAxes, fontsize=5.2)
    ax_fit.text(-0.17, 1.05, "d", transform=ax_fit.transAxes, fontsize=10, fontweight="bold", va="top")

    if transition_rows:
        values = np.array([[matrix.get((branch, transition), np.nan) for transition in transitions] for branch in branches])
        im = ax_heat.imshow(values, aspect="auto", vmin=0, vmax=1, cmap="Blues")
        ax_heat.set_title("POInT branch-state probabilities", fontsize=7, pad=4)
        ax_heat.set_yticks(np.arange(len(branches)))
        ax_heat.set_yticklabels(branches, fontsize=5.5)
        ax_heat.set_xticks(np.arange(len(transitions)))
        ax_heat.set_xticklabels(transitions, rotation=55, ha="right", fontsize=5.2)
        for iy in range(values.shape[0]):
            for ix in range(values.shape[1]):
                value = values[iy, ix]
                if not np.isnan(value):
                    ax_heat.text(ix, iy, "{:.2f}".format(value), ha="center", va="center", fontsize=4.9, color="#1B1B1B" if value < 0.65 else "white")
        cbar = fig.colorbar(im, ax=ax_heat, fraction=0.045, pad=0.03)
        cbar.set_label("mean posterior", fontsize=5.5)
        cbar.ax.tick_params(labelsize=5)
        conditional_note = "A/B labels are assembly copies; POInT reports posterior states, not verified gene loss."
    else:
        if not args.posterior_summary or not Path(args.posterior_summary).is_file():
            raise SystemExit("ERROR: no conditional matrix and no posterior summary available")
        counts, n_post = posterior_state_counts(args.posterior_summary)
        labels = list(counts)
        vals = [counts[label] / float(n_post) if n_post else 0.0 for label in labels]
        ax_heat.bar(np.arange(len(labels)), vals, color="#5B7FCA", width=0.72)
        ax_heat.set_title("POInT posterior tracking states", fontsize=7, pad=4)
        ax_heat.set_ylabel("fraction of pillars")
        ax_heat.set_ylim(0, 1)
        ax_heat.set_xticks(np.arange(len(labels)))
        ax_heat.set_xticklabels(labels, rotation=55, ha="right", fontsize=5.0)
        for ix, value in enumerate(vals):
            ax_heat.text(ix, value + 0.025, "{:.2f}".format(value), ha="center", va="bottom", fontsize=4.8)
        conditional_note = "Official POInT posterior fit; conditional branch-transition output unavailable in POInT v1.63 on this input."
    fig.text(0.01, 0.01, conditional_note, fontsize=5.4, color="#4D4D4D")
    prefix = outdir / "figure2_panel_d_point_official"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    with (outdir / "panel_d_source_data.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fit_rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(fit_rows)
    (outdir / "COMPLETE").write_text("Official POInT Fig. 2d plotted.\n", encoding="ascii")
    print("[OK] Wrote official POInT panel d: {}".format(outdir))


if __name__ == "__main__":
    main()
