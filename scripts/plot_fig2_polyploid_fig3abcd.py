#!/usr/bin/env python3
"""Render Figure 2 panels analogous to Fig. 3a-d without parentage claims."""

from __future__ import print_function

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.75
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False


COLORS = {
    "PZI": "#B64342",
    "PSP": "#0F4D92",
    "PPA": "#767676",
    "both_anchored": "#5B7FCA",
    "A_only_anchored": "#D9544D",
    "B_only_anchored": "#42949E",
    "ambiguous_multi_anchor": "#B89BD9",
    "neither_anchored": "#CFCECE",
    "both_anchored_in_both": "#5B7FCA",
    "same_A_only_candidate": "#D9544D",
    "same_B_only_candidate": "#42949E",
    "opposite_single_anchor": "#B89BD9",
    "same_ambiguous_state": "#767676",
    "discordant_or_uninformative": "#CFCECE",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ks", required=True)
    parser.add_argument("--topologies", required=True)
    parser.add_argument("--support", required=True, help="Final ASTRAL/gCF/sCF support table")
    parser.add_argument("--anchor-summary", required=True)
    parser.add_argument("--joint-states", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--ks-max", type=float, default=3.0)
    parser.add_argument("--kde-bandwidth", type=float, default=0.04)
    return parser.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def panel_label(ax, label):
    ax.text(-0.14, 1.05, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def kde(values, grid, bandwidth):
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return np.zeros_like(grid)
    total = np.zeros_like(grid)
    chunk = 5000
    for start in range(0, len(values), chunk):
        block = values[start:start + chunk]
        scaled = (grid[:, None] - block[None, :]) / bandwidth
        total += np.exp(-0.5 * scaled * scaled).sum(axis=1)
    return total / (len(values) * bandwidth * np.sqrt(2.0 * np.pi))


def plot_ks(ax, rows, ks_max, bandwidth):
    comparisons = {
        "pzi_hapA1_vs_pzi_hapB1": ("PZI A-B", "PZI"),
        "psp_hapA1_vs_psp_hapB1": ("PSP A-B", "PSP"),
        "ppa_hapA_vs_ppa_hapB": ("PPA A-B", "PPA"),
    }
    grouped = defaultdict(list)
    source_rows = []
    for row in rows:
        comparison = row.get("comparison", "")
        if comparison not in comparisons or row.get("kaks_status") != "ok":
            continue
        try:
            value = float(row["Ks"])
        except (KeyError, ValueError):
            continue
        if 0.0 < value < ks_max:
            grouped[comparison].append(value)
            source_rows.append(row)
    grid = np.linspace(0.0, ks_max, 550)
    for comparison, (label, color_key) in comparisons.items():
        values = grouped[comparison]
        if not values:
            raise ValueError("No filtered Ks values for " + comparison)
        ax.plot(grid, kde(values, grid, bandwidth), color=COLORS[color_key], lw=1.45, label="%s (n=%d)" % (label, len(values)))
    ax.set_xlabel("Synonymous substitutions per synonymous site (Ks)")
    ax.set_ylabel("Density")
    ax.set_title("Syntenic A-B divergence", loc="left", fontsize=8.2, fontweight="bold", pad=7)
    ax.set_xlim(0, ks_max)
    ax.legend(fontsize=5.6, loc="upper right", handlelength=1.8)
    return source_rows


def plot_topologies(ax, rows, support_rows):
    selected = [row for row in rows if row.get("contrast") == "pzi_psp_AB"]
    if len(selected) != 3:
        raise ValueError("Expected three pzi_psp_AB topology rows, found %d" % len(selected))
    labels = []
    values = []
    lows = []
    highs = []
    colors = []
    for row in selected:
        topology = row["topology"]
        if topology == "pzi_A+psp_A|pzi_B+psp_B":
            label = "A-A and B-B\n(shared-history topology)"
            color = COLORS["PSP"]
        elif topology == "pzi_A+pzi_B|psp_A+psp_B":
            label = "PZI and PSP\npaired within species"
            color = "#CFCECE"
        else:
            label = "Crossed\nsubgenomes"
            color = "#B89BD9"
        labels.append(label)
        values.append(100.0 * float(row["frequency"]))
        lows.append(100.0 * float(row["gene_bootstrap_ci_low"]))
        highs.append(100.0 * float(row["gene_bootstrap_ci_high"]))
        colors.append(color)
    x = np.arange(3)
    errors = np.vstack([np.asarray(values) - np.asarray(lows), np.asarray(highs) - np.asarray(values)])
    bars = ax.bar(x, values, color=colors, edgecolor="#4D4D4D", linewidth=0.5, yerr=errors, capsize=2.2, error_kw={"lw": 0.7})
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 2.0, "%.1f%%" % value, ha="center", va="bottom", fontsize=6.0)
    ax.set_ylim(0, min(108, max(values) + 9))
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5.7)
    ax.set_ylabel("Decisive gene trees (%)")
    ax.set_title("Subgenome topology frequencies", loc="left", fontsize=8.2, fontweight="bold", pad=7)
    support = {}
    for row in support_rows:
        taxa = set(filter(None, row.get("displayed_descendant_taxa", "").split(",")))
        if taxa == {"psp_A", "pzi_A"}:
            support["A"] = row
        elif taxa == {"psp_B", "pzi_B"}:
            support["B"] = row
    if set(support) == {"A", "B"}:
        annotation = "ASTRAL PP / gCF / sCF: A %.2f / %.1f / %.1f; B %.2f / %.1f / %.1f" % (
            float(support["A"]["astral_pp"]), float(support["A"]["gcf"]), float(support["A"]["likelihood_scf"]),
            float(support["B"]["astral_pp"]), float(support["B"]["gcf"]), float(support["B"]["likelihood_scf"]),
        )
        ax.text(0.0, 1.16, annotation, transform=ax.transAxes, fontsize=5.0, color="#4D4D4D", va="top")
    return selected


def plot_anchor_summary(ax, rows):
    species_order = ["PZI", "PSP", "PPA"]
    selected = [row for row in rows if row.get("dataset") == "main" and row.get("species") in species_order]
    if not selected:
        raise ValueError("No main PZI/PSP/PPA anchor-summary rows")
    by_species = defaultdict(dict)
    for row in selected:
        by_species[row["species"]][row["anchor_state"]] = float(row["fraction"])
    present_species = [species for species in species_order if species in by_species]
    state_order = ["both_anchored", "A_only_anchored", "B_only_anchored", "ambiguous_multi_anchor", "neither_anchored"]
    state_labels = ["Both", "A only", "B only", "Multi", "Neither"]
    x = np.arange(len(present_species))
    bottom = np.zeros(len(present_species))
    for state, state_label in zip(state_order, state_labels):
        values = np.asarray([by_species[species].get(state, 0.0) for species in present_species])
        ax.bar(x, values * 100.0, bottom=bottom * 100.0, color=COLORS[state], width=0.62, label=state_label, edgecolor="white", linewidth=0.35)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(present_species)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Conserved external loci (%)")
    ax.set_title("Candidate A/B retention at conserved anchors", loc="left", fontsize=8.2, fontweight="bold", pad=7)
    ax.legend(fontsize=5.2, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.24), columnspacing=0.75, handlelength=0.9)
    return selected


def plot_joint_states(ax, rows):
    selected = [row for row in rows if row.get("dataset") == "main"]
    if not selected:
        raise ValueError("No main joint-anchor states")
    order = [
        "both_anchored_in_both",
        "same_A_only_candidate",
        "same_B_only_candidate",
        "opposite_single_anchor",
        "same_ambiguous_state",
        "discordant_or_uninformative",
    ]
    labels = ["Both in both", "Shared A-only\ncandidate", "Shared B-only\ncandidate", "Opposite", "Same multi", "Discordant /\nuninformative"]
    counts = Counter(row["joint_anchor_state"] for row in selected)
    total = len(selected)
    values = [100.0 * counts.get(state, 0) / total for state in order]
    x = np.arange(len(order))
    bars = ax.bar(x, values, color=[COLORS[state] for state in order], edgecolor="#4D4D4D", linewidth=0.4)
    for bar, state in zip(bars, order):
        count = counts.get(state, 0)
        if count:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.4, str(count), ha="center", va="bottom", fontsize=5.2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=28, ha="right", fontsize=5.3)
    ax.set_ylabel("Conserved GTR loci (%)")
    ax.set_ylim(0, min(108, max(values) + 9))
    ax.set_title("PZI/PSP joint candidate retention states", loc="left", fontsize=8.2, fontweight="bold", pad=7)
    ax.text(0.01, -0.40, "A-only/B-only states require flanking synteny, sequence-search, gap and read-depth validation.", transform=ax.transAxes, fontsize=5.2, color="#4D4D4D", va="top")
    return selected


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: " + str(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(7.20, 5.25))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.18, 1.0], height_ratios=[1.0, 1.0], wspace=0.40, hspace=0.68)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    ks_rows = plot_ks(ax_a, read_tsv(args.ks), args.ks_max, args.kde_bandwidth)
    topology_rows = plot_topologies(ax_b, read_tsv(args.topologies), read_tsv(args.support))
    anchor_rows = plot_anchor_summary(ax_c, read_tsv(args.anchor_summary))
    joint_rows = plot_joint_states(ax_d, read_tsv(args.joint_states))
    for ax, label in zip([ax_a, ax_b, ax_c, ax_d], ["a", "b", "c", "d"]):
        panel_label(ax, label)

    fig.subplots_adjust(left=0.085, right=0.992, top=0.95, bottom=0.16)
    prefix = outdir / "figure2_polyploid_fig3abcd"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    write_tsv(outdir / "panel_a_ks_source_data.tsv", list(ks_rows[0].keys()), ks_rows)
    write_tsv(outdir / "panel_b_topology_source_data.tsv", list(topology_rows[0].keys()), topology_rows)
    write_tsv(outdir / "panel_c_anchor_source_data.tsv", list(anchor_rows[0].keys()), anchor_rows)
    write_tsv(outdir / "panel_d_joint_source_data.tsv", list(joint_rows[0].keys()), joint_rows)
    manifest = [
        {"panel": "a", "question": "A-B duplication signal", "boundary": "Ks is relative timing only."},
        {"panel": "b", "question": "Parallel subgenome histories", "boundary": "Does not identify parents."},
        {"panel": "c", "question": "A/B reference-anchor asymmetry", "boundary": "Not anchored is not a gene loss call."},
        {"panel": "d", "question": "Joint anchor states", "boundary": "Candidate states require sequence validation."},
    ]
    write_tsv(outdir / "figure_manifest.tsv", ["panel", "question", "boundary"], manifest)
    (outdir / "COMPLETE").write_text("Figure 2 polyploid Fig. 3a-d analogue completed.\n", encoding="ascii")
    print("[OK] Wrote figure bundle: " + str(outdir))


if __name__ == "__main__":
    main()
