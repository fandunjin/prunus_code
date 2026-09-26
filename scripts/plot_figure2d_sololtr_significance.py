#!/usr/bin/env python3
"""Figure 2d: solo-LTR/intact LTR-RT turnover with element-level significance.

Bars summarize per-haplotype genome-level ratios (mean with observed min-max
range). Significance brackets for P. zippeliana versus each comparator use
two-sided Fisher exact tests on pooled solo and intact LTR counts, with
Benjamini-Hochberg adjustment across all ten species pairs. The haplotype-level
exact permutation p-values are retained in the statistics table for reference.
"""

from __future__ import annotations

import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact

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
BRACKET_LEVELS = {("PZI", "PSP"): 3.62}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_counts(path: Path) -> tuple[dict[str, float], dict[str, list[float]], dict[str, list[int]]]:
    ratios: dict[str, float] = {}
    by_species: dict[str, list[float]] = defaultdict(list)
    pooled: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    group = {"pziA1": "PZI", "pziA2": "PZI", "pziB1": "PZI", "pziB2": "PZI",
             "pspA1": "PSP", "pspA2": "PSP", "pspB1": "PSP", "pspB2": "PSP",
             "ppaA": "PPA", "ppaB": "PPA", "pcam": "PCAM", "ppe": "PPE"}
    for row in read_tsv(path):
        target = row["target"]
        intact = int(row["quality_filtered_intact_records"])
        solo = int(row["final_solo_ltr_records"])
        ratios[target] = solo / intact
        pooled[group[target]][0] += solo
        pooled[group[target]][1] += intact
    for species, targets in SPECIES_TARGETS.items():
        by_species[species] = [ratios[target] for target in targets]
    return ratios, by_species, pooled


def exact_mean_permutation_p(x: list[float], y: list[float]) -> float | None:
    pooled = list(x) + list(y)
    if len(x) == 0 or len(y) == 0 or len(pooled) > 12:
        return None
    n_x = len(x)
    observed = abs(float(np.mean(x)) - float(np.mean(y)))
    extreme = 0
    total = 0
    for indices in itertools.combinations(range(len(pooled)), n_x):
        selected = set(indices)
        x_perm = [value for index, value in enumerate(pooled) if index in selected]
        y_perm = [value for index, value in enumerate(pooled) if index not in selected]
        if abs(float(np.mean(x_perm)) - float(np.mean(y_perm))) >= observed - 1e-12:
            extreme += 1
        total += 1
    return extreme / total if total else None


def bh_adjust(p_values: dict[tuple[str, str], float]) -> dict[tuple[str, str], float]:
    ordered = sorted(p_values.items(), key=lambda item: item[1])
    m = len(ordered)
    adjusted: dict[tuple[str, str], float] = {}
    running = 1.0
    for rank, (comparison, p_value) in reversed(list(enumerate(ordered, start=1))):
        running = min(running, p_value * m / rank)
        adjusted[comparison] = min(1.0, running)
    return adjusted


def star_label(q_value: float | None) -> str:
    if q_value is None:
        return ""
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return ""


def format_p(p: float) -> str:
    if p <= 0:
        return r"$P < 1 \times 10^{-300}$"
    exponent = int(math.floor(math.log10(p)))
    mantissa = p / (10 ** exponent)
    if exponent >= -4:
        return rf"$P = {p:.3g}$"
    return rf"$P = {mantissa:.1f} \times 10^{{{exponent}}}$"


def draw_bracket(ax: plt.Axes, x1: float, x2: float, y: float, drop: float, label: str) -> None:
    ax.plot([x1, x1, x2, x2], [y - drop, y, y, y - drop], color="#242424", lw=0.75, clip_on=False)
    if label:
        ax.text((x1 + x2) / 2, y + 0.04, label, ha="center", va="bottom", fontsize=6.0)


def main() -> None:
    root = Path("${PROJECT_ROOT}")
    input_path = root / "analysis/prunus_pan_te/followup_sololtr_v2/sololtr_summary/sololtr_qc.tsv"
    output_dir = root / "0903/figures/figure2_abc"

    ratios, by_species, pooled = load_counts(input_path)

    fisher: dict[tuple[str, str], float] = {}
    fisher_or: dict[tuple[str, str], float] = {}
    for a, b in itertools.combinations(SPECIES, 2):
        table = [[pooled[a][0], pooled[a][1]], [pooled[b][0], pooled[b][1]]]
        odds, p = fisher_exact(table, alternative="two-sided")
        fisher[(a, b)] = float(p)
        fisher_or[(a, b)] = float(odds)
    fisher_q = bh_adjust(fisher)

    perm: dict[tuple[str, str], float] = {}
    for a, b in itertools.combinations(SPECIES, 2):
        p = exact_mean_permutation_p(by_species[a], by_species[b])
        if p is not None:
            perm[(a, b)] = p
    perm_q = bh_adjust(perm) if perm else {}

    rows = []
    for a, b in itertools.combinations(SPECIES, 2):
        rows.append({
            "species_1": a,
            "species_2": b,
            "solo_1": pooled[a][0],
            "intact_1": pooled[a][1],
            "ratio_1": f"{pooled[a][0] / pooled[a][1]:.6f}",
            "solo_2": pooled[b][0],
            "intact_2": pooled[b][1],
            "ratio_2": f"{pooled[b][0] / pooled[b][1]:.6f}",
            "odds_ratio_fisher": f"{fisher_or[(a, b)]:.6f}",
            "fisher_p": f"{fisher[(a, b)]:.6g}",
            "fisher_q_BH": f"{fisher_q[(a, b)]:.6g}",
            "fisher_significance": star_label(fisher_q[(a, b)]),
            "permutation_p": f"{perm[(a, b)]:.6g}" if (a, b) in perm else "not tested",
            "permutation_q_BH": f"{perm_q[(a, b)]:.6g}" if (a, b) in perm_q else "not tested",
        })
    write_tsv(output_dir / "figure2d_sololtr_ratio_significance_statistics.tsv", rows, list(rows[0]))

    values = np.array([np.mean([ratios[t] for t in targets]) for _, targets, _ in PLOT_GROUPS], dtype=float)
    fig, ax = plt.subplots(figsize=(4.1, 3.9), facecolor="white")
    x = np.arange(len(PLOT_GROUPS), dtype=float)
    ax.bar(x, values, width=0.70, color=[c for _, _, c in PLOT_GROUPS], edgecolor="#242424", linewidth=0.65, zorder=3)
    for index, (_, targets, _) in enumerate(PLOT_GROUPS):
        pair_values = [ratios[t] for t in targets]
        center = float(index)
        ax.plot([center, center], [min(pair_values), max(pair_values)], color="#242424", lw=0.8, zorder=5)
        ax.plot([center - 0.24, center + 0.24], [min(pair_values), min(pair_values)], color="#242424", lw=0.8, zorder=5)
        ax.plot([center - 0.24, center + 0.24], [max(pair_values), max(pair_values)], color="#242424", lw=0.8, zorder=5)
    ax.axhline(1.0, color="#8B8B8B", lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_ylim(0, 5.25)
    ax.set_yticks(np.arange(0, 5.1, 0.5))
    ax.set_ylabel("solo-LTR / intact LTR-RT", labelpad=4)
    ax.set_xticks(x, [label for label, _, _ in PLOT_GROUPS], rotation=52, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=3, labelsize=6.5)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.grid(axis="y", color="#E3E3E3", lw=0.55, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("Solo-LTR turnover", loc="left", pad=8, fontsize=8, fontweight="bold")

    for (a, b), y in BRACKET_LEVELS.items():
        p = perm.get((a, b))
        if p is None:
            continue
        left_a, right_a = SPECIES_SPANS[a]
        left_b, right_b = SPECIES_SPANS[b]
        draw_bracket(ax, (left_a + right_a) / 2, (left_b + right_b) / 2, y, 0.07, f"$P = {p:.3f}$ (nominal)")

    ax.text(0.01, 1.01, "Haplotype-level exact permutation; nominal P", transform=ax.transAxes, fontsize=5.5, color="#555555", va="bottom")
    fig.subplots_adjust(left=0.17, right=0.98, top=0.90, bottom=0.30)
    for ext in ("svg", "pdf", "png", "tiff"):
        fig.savefig(output_dir / f"figure2d_sololtr_ratio_significance.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)

    legend = (
        "Solo-LTR/intact LTR-RT ratios for each phased genome or reference genome. Bars show the mean and the "
        "observed min-max range across phased haplotypes; singleton genomes have no empirical range. The bracket "
        "shows a haplotype-level exact permutation test on genome-level ratios for P. zippeliana versus P. spinulosa "
        "(nominal P = 0.029; Benjamini-Hochberg q = 0.286 across the ten species pairs); the remaining comparisons "
        "are descriptive because only one or two genomes were available per species. Element-level Fisher exact tests "
        "on pooled LTR counts are reported in the accompanying statistics table, but they assume independent LTR "
        "elements; because elements are strongly clustered within families (overdispersion phi ~ 234-444), their "
        "unadjusted P values are anticonservative and are not used for inference.\n"
    )
    (output_dir / "figure2d_sololtr_ratio_significance_legend.txt").write_text(legend, encoding="utf-8")
    print("figure2d written to", output_dir)


if __name__ == "__main__":
    main()
