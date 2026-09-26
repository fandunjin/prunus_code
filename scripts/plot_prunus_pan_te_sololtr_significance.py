#!/usr/bin/env python3
"""Add species-level pairwise significance brackets to the solo-LTR ratio plot.

The observational unit is one assembled genome/haplotype.  Exact permutation
tests are used because sample sizes are small and no normality assumption is
made.  PCAM and PPE are retained in the plot and included as descriptive
singleton-group comparisons.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 7,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    }
)


PLOT_GROUPS = [
    ("PZI-A", ["pziA1", "pziA2"], "#C44E52"),
    ("PZI-B", ["pziB1", "pziB2"], "#C44E52"),
    ("PSP-A", ["pspA1", "pspA2"], "#2A9D8F"),
    ("PSP-B", ["pspB1", "pspB2"], "#2A9D8F"),
    ("PPA-A", ["ppaA"], "#4C78A8"),
    ("PPA-B", ["ppaB"], "#4C78A8"),
    ("PCAM", ["pcam"], "#767676"),
    ("PPE", ["ppe"], "#B279A2"),
]
SPECIES = ["PZI", "PSP", "PPA", "PCAM", "PPE"]
SPECIES_TARGETS = {
    "PZI": ["pziA1", "pziA2", "pziB1", "pziB2"],
    "PSP": ["pspA1", "pspA2", "pspB1", "pspB2"],
    "PPA": ["ppaA", "ppaB"],
    "PCAM": ["pcam"],
    "PPE": ["ppe"],
}
SPECIES_SPANS = {"PZI": (0, 3), "PSP": (4, 7), "PPA": (8, 9), "PCAM": (10, 10), "PPE": (11, 11)}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_ratios(path: Path) -> tuple[dict[str, float], dict[str, list[float]]]:
    ratios: dict[str, float] = {}
    by_species: dict[str, list[float]] = defaultdict(list)
    for row in read_tsv(path):
        target = row["target"]
        intact = int(row["quality_filtered_intact_records"])
        solo = int(row["final_solo_ltr_records"])
        if intact <= 0:
            raise ValueError(f"Non-positive intact count for {target}")
        ratio = solo / intact
        ratios[target] = ratio
    for species, targets in SPECIES_TARGETS.items():
        missing = [target for target in targets if target not in ratios]
        if missing:
            raise ValueError(f"Missing targets for {species}: {', '.join(missing)}")
        by_species[species] = [ratios[target] for target in targets]
    return ratios, by_species


def exact_mean_permutation_p(x: list[float], y: list[float]) -> float:
    """Two-sided exact permutation p-value for an absolute mean difference."""
    pooled = list(x) + list(y)
    n_x = len(x)
    observed = abs(float(np.mean(x)) - float(np.mean(y)))
    extreme = 0
    total = 0
    for indices in itertools.combinations(range(len(pooled)), n_x):
        selected = set(indices)
        x_perm = [value for index, value in enumerate(pooled) if index in selected]
        y_perm = [value for index, value in enumerate(pooled) if index not in selected]
        difference = abs(float(np.mean(x_perm)) - float(np.mean(y_perm)))
        if difference >= observed - 1e-12:
            extreme += 1
        total += 1
    return extreme / total if total else 1.0


def bh_adjust(p_values: dict[tuple[str, str], float]) -> dict[tuple[str, str], float]:
    ordered = sorted(p_values.items(), key=lambda item: item[1])
    m = len(ordered)
    adjusted: dict[tuple[str, str], float] = {}
    running = 1.0
    for rank, (comparison, p_value) in reversed(list(enumerate(ordered, start=1))):
        running = min(running, p_value * m / rank)
        adjusted[comparison] = min(1.0, running)
    return adjusted


def star_label(q_value: float) -> str:
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return ""


def calculate_significance(by_species: dict[str, list[float]]) -> list[dict]:
    pairs: list[tuple[str, str]] = []
    p_values: dict[tuple[str, str], float] = {}
    for first, second in itertools.combinations(SPECIES, 2):
        x, y = by_species[first], by_species[second]
        # The exact randomization calculation is defined for any two non-empty
        # groups.  Singleton species are retained as descriptive comparisons;
        # their coarse p-value resolution is recorded in the test description.
        tested = len(x) >= 1 and len(y) >= 1
        pair = (first, second)
        if tested:
            pairs.append(pair)
            p_values[pair] = exact_mean_permutation_p(x, y)
    q_values = bh_adjust(p_values)
    rows: list[dict] = []
    for first, second in itertools.combinations(SPECIES, 2):
        x, y = by_species[first], by_species[second]
        pair = (first, second)
        p_value = p_values.get(pair)
        q_value = q_values.get(pair)
        rows.append(
            {
                "species_1": first,
                "species_2": second,
                "n_1": str(len(x)),
                "n_2": str(len(y)),
                "mean_ratio_1": f"{np.mean(x):.8f}",
                "mean_ratio_2": f"{np.mean(y):.8f}",
                "difference_mean_1_minus_2": f"{np.mean(x) - np.mean(y):.8f}",
                "test": (
                    "exact two-sided mean permutation"
                    if p_value is not None and len(x) >= 2 and len(y) >= 2
                    else "exact two-sided mean permutation (singleton group; descriptive)"
                    if p_value is not None
                    else "not tested: empty group"
                ),
                "p_value": f"{p_value:.8g}" if p_value is not None else "",
                "q_value_BH": f"{q_value:.8g}" if q_value is not None else "",
                "nominal_significance": star_label(p_value) if p_value is not None else "NA",
                "significance": star_label(q_value) if q_value is not None else "NA",
            }
        )
    return rows


def draw_bracket(ax: plt.Axes, x1: float, x2: float, y: float, drop: float, label: str) -> None:
    ax.plot([x1, x1, x2, x2], [y - drop, y, y, y - drop], color="#242424", lw=0.75, clip_on=False)
    if label:
        ax.text((x1 + x2) / 2, y + 0.025, label, ha="center", va="bottom", fontsize=8, fontweight="bold")


def plot(ratios: dict[str, float], stats: list[dict], output_prefix: Path, star_field: str) -> None:
    values = np.array([np.mean([ratios[target] for target in targets]) for _, targets, _ in PLOT_GROUPS], dtype=float)
    fig, ax = plt.subplots(figsize=(4.1, 3.9), facecolor="white")
    x = np.arange(len(PLOT_GROUPS), dtype=float)
    ax.bar(x, values, width=0.70, color=[color for _, _, color in PLOT_GROUPS], edgecolor="#242424", linewidth=0.65, zorder=3)

    # Bars summarize their target haplotypes; whiskers show min-max for n=2 groups.
    for index, (_, targets, _) in enumerate(PLOT_GROUPS):
        pair_values = [ratios[target] for target in targets]
        center = float(index)
        if len(pair_values) > 1:
            ax.plot([center, center], [min(pair_values), max(pair_values)], color="#242424", lw=0.8, zorder=5)
            ax.plot([center - 0.24, center + 0.24], [min(pair_values), min(pair_values)], color="#242424", lw=0.8, zorder=5)
            ax.plot([center - 0.24, center + 0.24], [max(pair_values), max(pair_values)], color="#242424", lw=0.8, zorder=5)
        else:
            ax.plot([center - 0.24, center + 0.24], [pair_values[0], pair_values[0]], color="#242424", lw=0.7, zorder=5)

    ax.axhline(1.0, color="#8B8B8B", lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_ylim(0, 4.30)
    ax.set_yticks(np.arange(0, 4.1, 0.5))
    ax.set_ylabel("solo-LTR / intact LTR-RT", labelpad=4)
    ax.set_xticks(x, [label for label, _, _ in PLOT_GROUPS], rotation=52, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=3, labelsize=6.5)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.grid(axis="y", color="#E3E3E3", lw=0.55, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("Solo-LTR turnover", loc="left", pad=8, fontsize=8, fontweight="bold")

    # Only significant brackets carry labels.  The current data produce one
    # nominally significant comparison, so the bracket layout remains compact.
    y_levels = {("PZI", "PSP"): 3.62, ("PZI", "PPA"): 4.04, ("PSP", "PPA"): 3.83}
    for row in stats:
        if not row["p_value"] or not row[star_field]:
            continue
        first, second = row["species_1"], row["species_2"]
        first_left, first_right = SPECIES_SPANS[first]
        second_left, second_right = SPECIES_SPANS[second]
        x1 = (first_left + first_right) / 2
        x2 = (second_left + second_right) / 2
        draw_bracket(ax, x1, x2, y_levels[(first, second)], 0.07, row[star_field])

    label = "Exact permutation; BH-adjusted" if star_field == "significance" else "Exact permutation; nominal p"
    ax.text(0.01, 1.01, label, transform=ax.transAxes, fontsize=5.5, color="#555555", va="bottom")
    fig.subplots_adjust(left=0.17, right=0.98, top=0.90, bottom=0.30)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{output_prefix}.svg", bbox_inches="tight")
    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{output_prefix}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{output_prefix}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--input", type=Path, default=root / "analysis/prunus_pan_te/followup_sololtr_v2/sololtr_summary/sololtr_qc.tsv")
    parser.add_argument("--output-dir", type=Path, default=root / "analysis/prunus_pan_te/figures_v3")
    parser.add_argument("--stem", default="supp_te_sololtr_ratio_with_significance")
    parser.add_argument("--star-mode", choices=("bh", "nominal"), default="nominal", help="Threshold stars using unadjusted exact p-values (default) or BH q-values")
    args = parser.parse_args()
    input_path = args.input if args.input.is_absolute() else root / args.input
    output_dir = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    ratios, by_species = load_ratios(input_path)
    stats = calculate_significance(by_species)
    output_prefix = output_dir / args.stem
    star_field = "significance" if args.star_mode == "bh" else "nominal_significance"
    plot(ratios, stats, output_prefix, star_field)
    fields = list(stats[0]) if stats else ["species_1", "species_2"]
    write_tsv(output_dir / f"{args.stem}_statistics.tsv", stats, fields)
    (output_dir / f"{args.stem}_legend.txt").write_text(
        "Solo-LTR/intact LTR-RT ratios are shown for each phased genome or reference genome. "
        "Bars with two phased genomes show the observed min-max range; singleton genomes have no empirical range. "
        "Species-level pairwise comparisons use exact two-sided permutation tests on genome-level ratios for all non-empty species pairs; comparisons involving singleton species are descriptive because their exact p-value resolution is coarse. Benjamini-Hochberg correction is applied across all ten species pairs. "
        "The figure star mode is selected by the command-line option: BH mode uses q<0.05, q<0.01 and q<0.001; nominal mode uses unadjusted p thresholds. "
        "PCAM and PPE are included as descriptive singleton-group comparisons. "
        "The phased haplotypes are assembly-level observations from the sampled genomes and should not be interpreted as independent biological replicates. Source data are provided as a TSV file.\n",
        encoding="utf-8",
    )
    print(f"Wrote {output_prefix}.[svg|pdf|png|tiff] and pairwise statistics TSV")


if __name__ == "__main__":
    main()
