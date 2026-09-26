#!/usr/bin/env python3
"""Render the integrated POInT evidence panel for Figure 2d.

The figure combines only native POInT-derived quantities:
1) symmetric versus copy-specific asymmetric-fractionation model fit;
2) the five shared/independent event partitions;
3) optional conditional branch-transition expectations;
4) posterior tracking-state frequencies for the shared three-taxon fit.

Conditional transitions are POInT model expectations.  They are never treated
as independently validated gene-loss counts.
"""

from __future__ import print_function

import argparse
import csv
import re
from collections import Counter
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


def state_label(state):
    match = re.fullmatch(r"PZI([01]{2})_PSP([01]{2})_PPA([01]{2})", state)
    if match:
        return "/".join(match.groups())
    return state


def display_branch(branch):
    return {
        "NONE": "PZI+PSP stem",
        "Root": "ancestral stem",
    }.get(branch, branch)


def require(path, label):
    if not Path(path).is_file() or Path(path).stat().st_size == 0:
        raise SystemExit("ERROR: missing {}: {}".format(label, path))


def write_source_data(path, rows):
    fields = sorted({key for row in rows for key in row})
    preferred = ["section", "model", "label", "state", "state_label", "count", "fraction",
                 "log_likelihood", "parameters", "pillars", "AIC", "BIC", "delta_AIC",
                 "delta_BIC", "max_posterior_mean"]
    fields = [field for field in preferred if field in fields] + [field for field in fields if field not in preferred]
    with Path(path).open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True, help="POInT symmetric/bias model summary")
    ap.add_argument("--partitions", required=True, help="Five shared/independent model summary")
    ap.add_argument("--posterior", required=True, help="Normalized POInT posterior summary")
    ap.add_argument("--transition", default=None, help="Optional normalized POInT conditional transition summary")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    require(args.fit, "POInT model-fit summary")
    require(args.partitions, "POInT partition summary")
    require(args.posterior, "POInT posterior summary")
    fit_rows = read_tsv(args.fit)
    partition_rows = read_tsv(args.partitions)
    posterior_rows = read_tsv(args.posterior)
    transition_rows = []
    if args.transition:
        require(args.transition, "POInT conditional transition summary")
        transition_rows = read_tsv(args.transition)
    if len(fit_rows) < 2:
        raise SystemExit("ERROR: expected at least symmetric and bias POInT fits")
    if len(partition_rows) != 5:
        raise SystemExit("ERROR: expected five partition-model rows")
    if not posterior_rows:
        raise SystemExit("ERROR: posterior summary is empty")
    if any(row.get("analysis_status", "converged").lower() not in ("converged", "complete", "converged_posterior_only") for row in fit_rows):
        raise SystemExit("ERROR: at least one POInT fit is not converged")

    # point_model_fit_all.tsv stores the model-process parameter count and
    # the native POInT likelihood.  The rooted three-taxon tree contributes
    # five optimized branch lengths to both nested models.
    fit_pillars = {int(float(row["pillars"])) for row in fit_rows}
    if len(fit_pillars) != 1:
        raise SystemExit("ERROR: model-fit rows were not evaluated on one common pillar set")
    fit_n = fit_pillars.pop()
    for row in fit_rows:
        process_parameters = int(float(row["model_parameters"]))
        total_parameters = process_parameters + 5
        log_likelihood = float(row["log_likelihood"])
        row["parameters"] = str(total_parameters)
        row["AIC"] = "{:.10f}".format(2.0 * total_parameters - 2.0 * log_likelihood)
        row["BIC"] = "{:.10f}".format(total_parameters * np.log(fit_n) - 2.0 * log_likelihood)
    fit_aic = np.array([float(row["AIC"]) for row in fit_rows], dtype=float)
    fit_min = np.min(fit_aic)
    for row, value in zip(fit_rows, fit_aic):
        row["delta_AIC"] = "{:.10f}".format(value - fit_min)

    partition_aic = np.array([float(row["AIC"]) for row in partition_rows], dtype=float)
    partition_bic = np.array([float(row["BIC"]) for row in partition_rows], dtype=float)
    for row, aic, bic in zip(partition_rows, partition_aic, partition_bic):
        row["delta_AIC"] = "{:.10f}".format(aic - np.min(partition_aic))
        row["delta_BIC"] = "{:.10f}".format(bic - np.min(partition_bic))

    counts = Counter(row.get("best_tracking", "") for row in posterior_rows)
    counts.pop("", None)
    if not counts:
        raise SystemExit("ERROR: no best_tracking states in posterior summary")
    ordered_states = sorted(counts, key=lambda state: (-counts[state], state))
    n_post = sum(counts.values())
    mean_max = float(np.mean([float(row["max_posterior"]) for row in posterior_rows if row.get("max_posterior")]))

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    if transition_rows:
        fig, axes = plt.subplots(1, 4, figsize=(10.4, 2.8), gridspec_kw={"wspace": 0.46})
        ax1, ax2, ax3, ax4 = axes
    else:
        fig, axes = plt.subplots(1, 3, figsize=(8.2, 2.8), gridspec_kw={"wspace": 0.42})
        ax1, ax2, ax4 = axes
        ax3 = None

    # The asymmetric model is POInT's nested test of unequal transitions from
    # the duplicated state to its two copy-specific singleton states.
    x1 = np.arange(len(fit_rows))
    fit_labels = ["SYM", "ASYM"]
    fit_colors = ["#B8B8B8", "#4C78A8"]
    ax1.bar(x1, fit_aic - fit_min, color=fit_colors, width=0.62)
    ax1.set_xticks(x1)
    ax1.set_xticklabels(fit_labels)
    ax1.set_ylabel("Delta AIC")
    ax1.set_ylim(bottom=0)
    ax1.text(-0.20, 1.05, "d", transform=ax1.transAxes, fontsize=10, fontweight="bold", va="top")
    ax1.text(0.98, 0.96, "n={:,}".format(n_post), transform=ax1.transAxes, ha="right", va="top", fontsize=5.5)

    # Event-partition comparison: bars are Delta AIC and open circles are
    # Delta BIC calculated on the same optimized-pillar sample.
    x2 = np.arange(len(partition_rows))
    p_labels = ["ALL", "PZI+PSP", "PZI+PPA", "PSP+PPA", "IND"]
    p_colors = ["#4C78A8", "#72B7B2", "#F2CF5B", "#E08E79", "#B0B0B0"]
    ax2.bar(x2, partition_aic - np.min(partition_aic), color=p_colors, width=0.62)
    ax2.scatter(x2, partition_bic - np.min(partition_bic), s=16, facecolors="white", edgecolors="#222222", linewidths=0.7, zorder=3)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(p_labels, rotation=32, ha="right")
    ax2.set_ylabel("Delta AIC / BIC")
    ax2.set_ylim(bottom=0)

    if ax3 is not None:
        # Summing posterior transition probabilities over sites gives an
        # expected transition count, not an observed gene-loss count.
        wanted = ("Dupl-->Copy1", "Dupl-->Copy2")
        branch_rows = {}
        branch_order = []
        for row in transition_rows:
            transition = row.get("transition", "")
            if transition not in wanted:
                continue
            branch = display_branch(row.get("branch", ""))
            if branch not in branch_rows:
                branch_rows[branch] = {}
                branch_order.append(branch)
            branch_rows[branch][transition] = float(row["mean_probability"]) * fit_n
        if not branch_order:
            raise SystemExit("ERROR: conditional summary lacks Dupl-->Copy1/Copy2 transitions")
        x3 = np.arange(len(branch_order))
        width = 0.34
        copy1 = [branch_rows[branch].get("Dupl-->Copy1", 0.0) for branch in branch_order]
        copy2 = [branch_rows[branch].get("Dupl-->Copy2", 0.0) for branch in branch_order]
        ax3.bar(x3 - width / 2, copy1, width, color="#4C78A8", label="C1")
        ax3.bar(x3 + width / 2, copy2, width, color="#E08E79", label="C2")
        ax3.set_xticks(x3)
        ax3.set_xticklabels(branch_order, rotation=32, ha="right")
        ax3.set_ylabel("Expected transitions")
        ax3.legend(loc="upper right", fontsize=5.4, handlelength=1.0)
        ax3.set_ylim(bottom=0)

    # State frequencies: labels preserve POInT's binary-copy state coding.
    fractions = np.array([counts[state] / float(n_post) for state in ordered_states])
    labels3 = [state_label(state) for state in ordered_states]
    y3 = np.arange(len(ordered_states))
    ax4.barh(y3, fractions, color="#5B7FCA", height=0.62)
    ax4.set_yticks(y3)
    ax4.set_yticklabels(labels3, fontsize=5.4)
    ax4.invert_yaxis()
    ax4.set_xlim(0, max(0.5, float(np.max(fractions)) * 1.22))
    ax4.set_xlabel("Fraction of pillars")
    ax4.text(0.98, 0.02, "mean max posterior={:.3f}".format(mean_max), transform=ax4.transAxes, ha="right", va="bottom", fontsize=5.2)

    for ax in axes:
        ax.tick_params(length=3, width=0.7)
        ax.spines["left"].set_linewidth(0.75)
        ax.spines["bottom"].set_linewidth(0.75)

    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.25, top=0.97)
    prefix = outdir / "figure2_panel_d_point_integrated"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)

    source_rows = []
    for row in fit_rows:
        source_rows.append({"section": "copy_bias_fit", **row})
    for row in partition_rows:
        source_rows.append({"section": "event_partition", **row})
    for row in transition_rows:
        source_rows.append({"section": "conditional_transition", **row})
    for state in ordered_states:
        source_rows.append({
            "section": "posterior_state",
            "state": state,
            "state_label": state_label(state),
            "count": str(counts[state]),
            "fraction": "{:.10f}".format(counts[state] / float(n_post)),
            "pillars": str(n_post),
            "max_posterior_mean": "{:.10f}".format(mean_max),
        })
    write_source_data(outdir / "panel_d_integrated_source_data.tsv", source_rows)
    boundary = [
        "Native POInT evidence only. The asymmetric model supports copy-specific asymmetric fractionation, not direct modern A/B dominance.",
        "GTR is an external anchor for pillar homology and is not a POInT likelihood terminal.",
    ]
    if transition_rows:
        boundary.append("Conditional values are model-based expected state transitions; they are not independently verified branch-specific gene-loss counts.")
    else:
        boundary.append("POInT conditional branch-transition output was unavailable; no branch-specific model transition expectations are plotted.")
    (outdir / "analysis_boundary.txt").write_text("\n".join(boundary) + "\n", encoding="ascii")
    (outdir / "COMPLETE").write_text("Integrated POInT Figure 2d completed.\n", encoding="ascii")
    print("[OK] Wrote integrated POInT panel d: {}".format(outdir))


if __name__ == "__main__":
    main()
