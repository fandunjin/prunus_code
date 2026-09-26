#!/usr/bin/env python3
"""Create publication-ready Prunus pan-TE main and supplementary panels."""

import argparse
import collections
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


TARGET_ORDER = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]

COMPOSITION_ORDER = [
    "LTR/Gypsy", "LTR/Copia", "Other LTR", "TIR/DNA",
    "Helitron", "LINE", "SINE", "Other TE", "Non-TE",
]

COMPOSITION_COLORS = {
    "LTR/Gypsy": "#4C78A8",
    "LTR/Copia": "#72A6C8",
    "Other LTR": "#A6C8DF",
    "TIR/DNA": "#E39C45",
    "Helitron": "#D56A5D",
    "LINE": "#59A14F",
    "SINE": "#8BC17C",
    "Other TE": "#B8A1CF",
    "Non-TE": "#D8D8D8",
}

GROUP_COLORS = {
    "PZI": "#C44E52",
    "PZI-A": "#C44E52",
    "PZI-B": "#8E3B5B",
    "PSP": "#2A9D8F",
    "PSP-A": "#2A9D8F",
    "PSP-B": "#187B73",
    "PPA": "#4C78A8",
    "PCAM": "#7A7A7A",
    "PPE": "#B279A2",
}

AGE_GROUP_TARGETS = {
    "PZI": ["pziA1", "pziA2", "pziB1", "pziB2"],
    "PSP": ["pspA1", "pspA2", "pspB1", "pspB2"],
    "PPA": ["ppaA", "ppaB"],
    "PCAM": ["pcam"],
    "PPE": ["ppe"],
}


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_figure(fig, output_dir, stem):
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix, kwargs in (
        ("svg", {}),
        ("pdf", {}),
        ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}}),
        ("png", {"dpi": 300}),
    ):
        path = output_dir / f"{stem}.{suffix}"
        fig.savefig(path, bbox_inches="tight", facecolor="white", **kwargs)
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"figure export failed: {path}")
        paths.append(path)
    plt.close(fig)
    return paths


def selected_drivers(rows, limit=None):
    selected = [row for row in rows if row["driver_selected"] == "1"]
    selected.sort(key=lambda row: int(row["positive_rank"]))
    return selected[:limit] if limit else selected


def family_label(row):
    class_name = row["class_family"].replace("LTR/", "").replace("DNA/", "")
    return f"{row['pan_family']}  {class_name}"


def panel_label(ax, label):
    ax.text(-0.10, 1.04, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def plot_composition(rows, output_dir):
    values = collections.defaultdict(dict)
    genome = {}
    for row in rows:
        values[row["target"]][row["category"]] = int(row["bp"]) / 1e6
        genome[row["target"]] = int(row["genome_bp"]) / 1e6

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = np.arange(len(TARGET_ORDER))
    left = np.zeros(len(TARGET_ORDER))
    for category in COMPOSITION_ORDER:
        data = np.array([values[target].get(category, 0.0) for target in TARGET_ORDER])
        ax.barh(y, data, left=left, height=0.72, color=COMPOSITION_COLORS[category], label=category, linewidth=0)
        left += data
    for index, target in enumerate(TARGET_ORDER):
        ax.text(genome[target] + max(genome.values()) * 0.012, index, f"{genome[target]:.0f}", va="center", fontsize=6)
    ax.set_yticks(y, TARGET_ORDER)
    ax.invert_yaxis()
    ax.set_xlabel("Sequence length (Mb)")
    ax.set_xlim(0, max(genome.values()) * 1.11)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(ncol=5, loc="lower center", bbox_to_anchor=(0.5, 1.01), columnspacing=0.9, handlelength=1.2)
    panel_label(ax, "b")
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.83)
    return save_figure(fig, output_dir, "figure1_te_composition")


def plot_drivers(rows, output_dir, max_drivers):
    drivers = selected_drivers(rows, max_drivers)
    if not drivers:
        raise ValueError("no driver families were selected")
    drivers = list(reversed(drivers))
    y = np.arange(len(drivers))
    delta_a = np.array([int(row["delta_a_bp"]) / 1e6 for row in drivers])
    delta_b = np.array([int(row["delta_b_bp"]) / 1e6 for row in drivers])

    fig_height = max(3.2, 0.28 * len(drivers) + 1.5)
    fig, ax = plt.subplots(figsize=(7.2, fig_height))
    for index, (a_value, b_value) in enumerate(zip(delta_a, delta_b)):
        ax.plot([a_value, b_value], [index, index], color="#B7B7B7", linewidth=1.1, zorder=1)
    ax.scatter(delta_a, y, s=23, color=GROUP_COLORS["PZI-A"], label="A: PZI - PSP", zorder=3)
    ax.scatter(delta_b, y, s=23, color=GROUP_COLORS["PZI-B"], label="B: PZI - PSP", zorder=3)
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_yticks(y, [family_label(row) for row in drivers])
    ax.set_xlabel("Excess pan-TE component sequence in PZI (Mb per haploid genome)")
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right")
    panel_label(ax, "c")
    fig.subplots_adjust(left=0.30, right=0.98, bottom=0.14, top=0.96)
    return save_figure(fig, output_dir, "figure1_te_driver_excess")


def aggregate_divergence(rows, family, max_divergence):
    by_target = collections.defaultdict(lambda: collections.Counter())
    bins = set()
    for row in rows:
        if row["pan_family"] != family:
            continue
        x_value = float(row["divergence_bin_mid_pct"])
        if x_value > max_divergence:
            continue
        by_target[row["target"]][x_value] += int(row["exclusive_bp"])
        bins.add(x_value)
    ordered_bins = sorted(bins)
    grouped = {}
    for group, targets in AGE_GROUP_TARGETS.items():
        grouped[group] = np.array([
            np.mean([by_target[target][x_value] for target in targets]) / 1e6
            for x_value in ordered_bins
        ])
    return np.array(ordered_bins), grouped


def add_age_axis(ax, substitution_rate):
    if not substitution_rate:
        return
    def divergence_to_age(value):
        return value / 100 / (2 * substitution_rate) / 1e6
    def age_to_divergence(value):
        return value * 1e6 * 2 * substitution_rate * 100
    secondary = ax.secondary_xaxis("top", functions=(divergence_to_age, age_to_divergence))
    secondary.set_xlabel("Estimated insertion time (Ma)")


def plot_divergence(divergence_rows, driver_rows, output_dir, max_families, max_divergence, substitution_rate):
    drivers = [
        row for row in selected_drivers(driver_rows)
        if not row["class_family"].startswith("LTR/")
        or int(row.get("library_sequence_length_bp") or 0) >= 500
    ]
    families = [row["pan_family"] for row in drivers[:max_families]]
    if not families:
        raise ValueError("no driver families available for divergence plots")
    ncols = 2
    nrows = math.ceil(len(families) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.2, 2.35 * nrows), sharex=True, squeeze=False)
    for axis, family in zip(axes.flat, families):
        x_values, grouped = aggregate_divergence(divergence_rows, family, max_divergence)
        for group in AGE_GROUP_TARGETS:
            axis.plot(x_values, grouped[group], color=GROUP_COLORS[group], linewidth=1.25, label=group)
        driver = next(row for row in drivers if row["pan_family"] == family)
        axis.set_title(family_label(driver), fontsize=7, loc="left")
        axis.set_ylabel("Exclusively assigned sequence (Mb/bin)")
        axis.grid(color="#ECECEC", linewidth=0.45)
        add_age_axis(axis, substitution_rate)
    for axis in axes[-1, :]:
        axis.set_xlabel("RepeatMasker divergence (%)")
    for axis in axes.flat[len(families):]:
        axis.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=5, loc="upper center", bbox_to_anchor=(0.53, 1.01))
    axes.flat[0].text(
        0.01, 0.97, "d", transform=axes.flat[0].transAxes,
        fontsize=9, fontweight="bold", va="top",
    )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.10, top=0.88, hspace=0.45, wspace=0.28)
    return save_figure(fig, output_dir, "figure1_te_divergence")


def plot_driver_heatmap(driver_rows, output_dir, max_drivers):
    drivers = selected_drivers(driver_rows, max_drivers)
    if not drivers:
        return []
    matrix = np.array([
        [int(row[f"{target}_bp"]) / 1e6 for target in TARGET_ORDER]
        for row in drivers
    ])
    transformed = np.log10(matrix + 0.001)
    fig_height = max(3.5, 0.27 * len(drivers) + 1.6)
    fig, ax = plt.subplots(figsize=(7.2, fig_height))
    image = ax.imshow(transformed, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(TARGET_ORDER)), TARGET_ORDER, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(drivers)), [family_label(row) for row in drivers])
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("log10(Mb + 0.001)")
    ax.set_xlabel("Haploid genome")
    fig.subplots_adjust(left=0.28, right=0.93, bottom=0.20, top=0.97)
    return save_figure(fig, output_dir, "supp_te_driver_heatmap")


def plot_mapping_qc(rows, output_dir):
    mapped = [float(row["mapped_fraction"]) * 100 for row in rows]
    identity = [float(row["median_mapped_identity"]) * 100 if row["median_mapped_identity"] else np.nan for row in rows]
    x = np.arange(len(rows))
    colors = [GROUP_COLORS.get(row["group"], "#777777") for row in rows]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.8), sharex=True)
    axes[0].bar(x, mapped, color=colors, width=0.72)
    axes[0].set_ylabel("Families mapped (%)")
    axes[0].set_ylim(0, 105)
    axes[1].scatter(x, identity, color=colors, s=25)
    axes[1].set_ylabel("Median identity (%)")
    axes[1].set_xticks(x, [row["target"] for row in rows], rotation=45, ha="right")
    for axis in axes:
        axis.grid(axis="y", color="#E8E8E8", linewidth=0.5)
        axis.set_axisbelow(True)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.22, top=0.97, hspace=0.16)
    return save_figure(fig, output_dir, "supp_te_library_mapping_qc")


def plot_pan_vs_independent(rows, output_dir):
    x_values = np.array([int(row["independent_edta_union_bp"]) / 1e6 for row in rows])
    y_values = np.array([int(row["pan_union_bp"]) / 1e6 for row in rows])
    valid = (x_values > 0) | (y_values > 0)
    x_values = x_values[valid]
    y_values = y_values[valid]
    valid_rows = [row for row, keep in zip(rows, valid) if keep]
    if not valid_rows:
        raise ValueError("pan-versus-independent table contains no non-zero comparisons")
    upper = max(float(x_values.max()), float(y_values.max())) * 1.05
    fig, ax = plt.subplots(figsize=(4.8, 4.5))
    for group in ("PZI", "PSP", "PPA", "PCAM", "PPE"):
        indices = [index for index, row in enumerate(valid_rows) if row["group"] == group]
        if not indices:
            continue
        ax.scatter(x_values[indices], y_values[indices], s=22, alpha=0.78, color=GROUP_COLORS[group], label=group)
    ax.plot([0, upper], [0, upper], color="#555555", linestyle="--", linewidth=0.8)
    ax.set_xlim(0, upper)
    ax.set_ylim(0, upper)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Independent EDTA annotation (Mb)")
    ax.set_ylabel("Standardized pan-library annotation (Mb)")
    ax.grid(color="#ECECEC", linewidth=0.45)
    ax.legend(loc="upper left")
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.14, top=0.97)
    return save_figure(fig, output_dir, "supp_te_pan_vs_independent")


def plot_intact_ltr(
    intact_rows, driver_rows, output_dir, max_families, max_divergence,
    min_overlap_fraction, min_intact_elements, min_ltr_library_length,
):
    ltr_drivers = [
        row for row in selected_drivers(driver_rows)
        if row["class_family"].startswith("LTR/")
        and int(row.get("library_sequence_length_bp") or 0) >= min_ltr_library_length
    ]
    supported = collections.Counter()
    for row in intact_rows:
        if (
            row["assignment_method"] == "maximum_coordinate_overlap"
            and float(row["overlap_fraction"]) >= min_overlap_fraction
        ):
            supported[row["pan_family"]] += 1
    families = [
        row["pan_family"] for row in ltr_drivers
        if supported[row["pan_family"]] >= min_intact_elements
    ][:max_families]
    if not families:
        print(
            "WARNING: no selected LTR driver meets the library-length and intact-element thresholds; "
            "supp_te_intact_ltr_divergence was not exported",
            file=sys.stderr,
        )
        return []
    ncols = 2
    nrows = math.ceil(len(families) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.2, 2.35 * nrows), sharex=True, squeeze=False)
    bins = np.linspace(0, max_divergence, 31)
    for axis, family in zip(axes.flat, families):
        for group, targets in AGE_GROUP_TARGETS.items():
            values = [
                float(row["ltr_divergence"]) * 100
                for row in intact_rows
                if row["pan_family"] == family and row["target"] in targets
                and row["assignment_method"] == "maximum_coordinate_overlap"
                and float(row["overlap_fraction"]) >= min_overlap_fraction
                and float(row["ltr_divergence"]) * 100 <= max_divergence
            ]
            if values:
                counts, edges = np.histogram(values, bins=bins)
                axis.step((edges[:-1] + edges[1:]) / 2, counts / len(targets), where="mid", color=GROUP_COLORS[group], linewidth=1.2, label=group)
        driver = next(row for row in ltr_drivers if row["pan_family"] == family)
        axis.set_title(f"{family_label(driver)}  n={supported[family]}", fontsize=7, loc="left")
        axis.set_ylabel("Mean intact LTR count/genome")
        axis.grid(color="#ECECEC", linewidth=0.45)
    for axis in axes[-1, :]:
        axis.set_xlabel("LTR-pair divergence (%)")
    for axis in axes.flat[len(families):]:
        axis.set_visible(False)
    legend_by_label = {}
    for axis in axes.flat[:len(families)]:
        handles, labels = axis.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            legend_by_label.setdefault(label, handle)
    if legend_by_label:
        fig.legend(
            legend_by_label.values(), legend_by_label.keys(), ncol=5,
            loc="upper center", bbox_to_anchor=(0.53, 1.01),
        )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.10, top=0.88, hspace=0.40, wspace=0.28)
    return save_figure(fig, output_dir, "supp_te_intact_ltr_divergence")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--postprocess-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-drivers", type=int, default=18)
    parser.add_argument("--age-families", type=int, default=4)
    parser.add_argument("--max-divergence", type=float, default=50.0)
    parser.add_argument("--min-ltr-overlap-fraction", type=float, default=0.50)
    parser.add_argument("--min-intact-elements", type=int, default=20)
    parser.add_argument("--min-ltr-library-length", type=int, default=500)
    parser.add_argument("--substitution-rate", type=float)
    args = parser.parse_args()
    if not 0 <= args.min_ltr_overlap_fraction <= 1:
        parser.error("--min-ltr-overlap-fraction must be between 0 and 1")
    if args.min_intact_elements < 1:
        parser.error("--min-intact-elements must be positive")
    if args.min_ltr_library_length < 1:
        parser.error("--min-ltr-library-length must be positive")

    composition = read_tsv(args.postprocess_dir / "te_composition.tsv")
    drivers = read_tsv(args.postprocess_dir / "te_driver_families.tsv")
    divergence = read_tsv(args.postprocess_dir / "te_repeatmasker_divergence.tsv")
    intact_ltr = read_tsv(args.postprocess_dir / "te_intact_ltr.tsv")
    mapping_qc = read_tsv(args.postprocess_dir / "te_library_mapping_qc.tsv")
    pan_vs_independent = read_tsv(args.postprocess_dir / "te_pan_vs_independent.tsv")

    manifest_path = args.postprocess_dir / "postprocess_manifest.json"
    with open(manifest_path) as handle:
        manifest = json.load(handle)
    substitution_rate = args.substitution_rate or manifest.get("substitution_rate_substitutions_per_site_per_year")

    exported = []
    exported += plot_composition(composition, args.output_dir)
    exported += plot_drivers(drivers, args.output_dir, args.max_drivers)
    exported += plot_divergence(divergence, drivers, args.output_dir, args.age_families, args.max_divergence, substitution_rate)
    exported += plot_driver_heatmap(drivers, args.output_dir, args.max_drivers)
    exported += plot_mapping_qc(mapping_qc, args.output_dir)
    exported += plot_pan_vs_independent(pan_vs_independent, args.output_dir)
    exported += plot_intact_ltr(
        intact_ltr, drivers, args.output_dir, args.age_families,
        min(args.max_divergence, 20.0), args.min_ltr_overlap_fraction,
        args.min_intact_elements, args.min_ltr_library_length,
    )

    print("=== Prunus pan-TE figure export completed ===")
    print(f"Exported files: {len(exported)}")
    for path in exported:
        print(path)


if __name__ == "__main__":
    main()
