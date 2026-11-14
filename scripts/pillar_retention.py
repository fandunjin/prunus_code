#!/usr/bin/env python3
"""Plot Fig. 4a: GTR-centred ancestral syntenic pillar retention.

The main panel uses contrast-stretched variable-height density tracks.  Each
column is a local window of approximately ten GTR-centred pillars; its height
expands the 75–100% retained-fraction interval for physical A or B anchors.
The paired A/B state summary retains the exact pillar-level counts at right.
"""

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd


# Nature-figure requirements: editable text in the SVG and TrueType PDF text.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 7.0,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "legend.frameon": False,
    }
)


import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pillars", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
INPUT = args.pillars
FIGURE_DIR = args.out
OUT_BASE = FIGURE_DIR / "figure4_panel_a_gtr_pillars"
SOURCE = FIGURE_DIR / "figure4_panel_a_source_data.tsv"
DENSITY_SOURCE = FIGURE_DIR / "figure4_panel_a_density_source.tsv"

SPECIES = ["PZI", "PSP", "PPA"]
STATE_LABELS = {
    "both": "A and B retained",
    "a_only": "A retained only",
    "b_only": "B retained only",
    "none": "Both lost",
}
COLORS = {
    "a_row": "#D97B70",
    "b_row": "#4E78A7",
    "lost": "#E7E7E7",
    "both": "#395C78",
    "a_only": "#D97B70",
    "b_only": "#4E78A7",
    "none": "#C6C6C6",
    "chrom_a": "#D9E3EA",
    "chrom_b": "#EEF1F3",
    "outline": "#383838",
}
DISPLAY_FLOOR = 0.75
SPECIES_COLORS = {
    "PZI": "#C65F5A",
    "PSP": "#4E78A7",
    "PPA": "#4F8A6E",
}


def state_name(a: int, b: int) -> str:
    if a and b:
        return "both"
    if a and not b:
        return "a_only"
    if b and not a:
        return "b_only"
    return "none"


def main() -> None:
    data = pd.read_csv(INPUT, sep="\t")
    required = [
        "gtr_global_index",
        "gtr_seqid",
        *[f"{sp}_{copy}_state" for sp in SPECIES for copy in ["A", "B"]],
    ]
    missing = [c for c in required if c not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    data = data.sort_values("gtr_global_index", kind="stable").reset_index(drop=True)
    n = len(data)
    if n == 0:
        raise ValueError("No syntenic pillars found")
    expected = os.environ.get("EXPECTED_PILLARS")
    if expected is not None and n != int(expected):
        raise ValueError(f"Expected {int(expected):,} GTR-centred pillars, found {n:,}")

    # Normalize the paired states and retain the source coordinates for auditability.
    for sp in SPECIES:
        data[f"{sp}_A"] = data[f"{sp}_A_state"].fillna(0).astype(int).clip(0, 1)
        data[f"{sp}_B"] = data[f"{sp}_B_state"].fillna(0).astype(int).clip(0, 1)
        data[f"{sp}_paired_state"] = [
            state_name(a, b) for a, b in zip(data[f"{sp}_A"], data[f"{sp}_B"])
        ]

    source_cols = [
        "pillar_id",
        "gtr_global_index",
        "gtr_seqid",
        "gtr_gene",
        *[f"{sp}_{copy}" for sp in SPECIES for copy in ["A", "B"]],
        *[f"{sp}_paired_state" for sp in SPECIES],
    ]
    SOURCE.parent.mkdir(parents=True, exist_ok=True)
    data[source_cols].to_csv(SOURCE, sep="\t", index=False)

    # Main panel: variable-height A/B density tracks, ordered by GTR coordinate.
    # A bin contains approximately 10 pillars. The display height expands the
    # 75–100% interval to make local differences visible without changing the
    # exact fractions written to the source-data table and summary panel.
    n_bins = int(np.ceil(n / 10))
    edges = np.linspace(0, n, n_bins + 1, dtype=int)
    density_rows = []
    for sp in SPECIES:
        for copy in ["A", "B"]:
            values = data[f"{sp}_{copy}"].to_numpy(dtype=float)
            for bin_id, (start, end) in enumerate(zip(edges[:-1], edges[1:]), start=1):
                n_in_bin = end - start
                retained_count = int(values[start:end].sum())
                density_rows.append(
                    {
                        "species": sp,
                        "track": copy,
                        "bin": bin_id,
                        "pillar_start": int(start),
                        "pillar_end": int(end),
                        "n_pillars": n_in_bin,
                        "retained_count": retained_count,
                        "retention_fraction": retained_count / n_in_bin,
                    }
                )
    pd.DataFrame(density_rows).to_csv(DENSITY_SOURCE, sep="\t", index=False)

    fig = plt.figure(figsize=(8.25, 4.65), facecolor="white")
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[5.75, 1.75],
        height_ratios=[1.0, 1.18],
        left=0.095,
        right=0.985,
        bottom=0.12,
        top=0.80,
        hspace=0.45,
        wspace=0.38,
    )
    ax = fig.add_subplot(gs[:, 0])
    ax_diff = fig.add_subplot(gs[0, 1])
    ax_bar = fig.add_subplot(gs[1, 1])

    track_height = 0.76
    baselines = {(sp, "A"): 4.95 - i for i, sp in enumerate(SPECIES)}
    baselines.update({(sp, "B"): 1.55 - i for i, sp in enumerate(SPECIES)})
    bin_mid = (edges[:-1] + edges[1:]) / 2.0
    smooth_bins = 11

    def edge_smooth(values: np.ndarray, width: int) -> np.ndarray:
        if width <= 1:
            return values
        radius = width // 2
        padded = np.pad(values, (radius, radius), mode="edge")
        return np.convolve(padded, np.ones(width) / width, mode="valid")

    wave_values = {}
    # A light base and a waveform make both the local profile and the direction
    # of change readable without treating the smoothed display as new counts.
    for sp in SPECIES:
        for copy in ["A", "B"]:
            baseline = baselines[(sp, copy)]
            ax.add_patch(
                Rectangle(
                    (0, baseline),
                    n,
                    track_height,
                    facecolor="#F1F3F4",
                    edgecolor="none",
                    zorder=1,
                )
            )
            ax.plot(
                [0, n],
                [baseline + track_height * 0.5] * 2,
                color="#A6AEB6",
                linewidth=0.35,
                linestyle=(0, (1.5, 1.5)),
                alpha=0.75,
                zorder=2,
            )
            values = data[f"{sp}_{copy}"].to_numpy(dtype=float)
            local_values = np.array(
                [values[start:end].mean() for start, end in zip(edges[:-1], edges[1:])]
            )
            smoothed = edge_smooth(local_values, smooth_bins)
            wave_values[(sp, copy)] = smoothed
            display_fraction = np.clip(
                (smoothed - DISPLAY_FLOOR) / (1.0 - DISPLAY_FLOOR), 0.0, 1.0
            )
            y = baseline + track_height * display_fraction
            color = SPECIES_COLORS[sp]
            ax.fill_between(
                bin_mid, baseline, y, color=color, alpha=0.27, linewidth=0, zorder=3
            )
            ax.plot(
                bin_mid,
                y,
                color=color,
                linewidth=0.85,
                linestyle="-" if copy == "A" else "--",
                zorder=4,
            )

    ax.set_xlim(0, n)
    ax.set_ylim(-0.55, 6.25)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlabel("GTR-ordered syntenic pillars", labelpad=9, fontsize=7.2)

    # Species labels and group labels make the A-first/B-second layout explicit.
    for i, sp in enumerate(SPECIES):
        for copy in ["A", "B"]:
            center = baselines[(sp, copy)] + track_height / 2
            ax.text(
                -0.012,
                center,
                sp,
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=6.5,
                fontweight="bold",
                color=SPECIES_COLORS[sp],
            )

    ax.text(
        -0.068,
        baselines[("PZI", "A")] + 0.45,
        "A",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=8.2,
        fontweight="bold",
        color=COLORS["a_row"],
    )
    ax.text(
        -0.068,
        baselines[("PZI", "B")] + 0.45,
        "B",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=8.2,
        fontweight="bold",
        color=COLORS["b_row"],
    )
    for y in [2.62, 0.62]:
        ax.axhline(y, color=COLORS["outline"], linewidth=0.3, alpha=0.65, zorder=6)

    # GTR chromosome strip. Pillars are already ordered by gtr_global_index.
    chrom_groups = data.groupby("gtr_seqid", sort=False).size()
    start = 0
    for idx, (chrom, count) in enumerate(chrom_groups.items()):
        ax.add_patch(
            Rectangle(
                (start, 5.97),
                count,
                0.25,
                facecolor=COLORS["chrom_a" if idx % 2 == 0 else "chrom_b"],
                edgecolor=COLORS["outline"],
                linewidth=0.25,
                clip_on=False,
                zorder=8,
            )
        )
        label = chrom.replace("GTR__", "").replace("Chr", "Chr ")
        ax.text(
            start + count / 2,
            6.095,
            label,
            ha="center",
            va="center",
            fontsize=5.8,
            color=COLORS["outline"],
            zorder=9,
        )
        if start > 0:
            ax.axvline(
                start, color=COLORS["outline"], linewidth=0.45, alpha=0.65, zorder=7
            )
        start += count

    ax.text(
        0,
        6.17,
        "Smoothed local retention waves",
        ha="left",
        va="bottom",
        fontsize=7.4,
        fontweight="bold",
        color=COLORS["outline"],
    )
    ax.text(
        n,
        6.17,
        f"n = {n:,}",
        ha="right",
        va="bottom",
        fontsize=6.8,
        color=COLORS["outline"],
    )
    fig.text(
        0.095,
        0.815,
        "Wave height expands 75–100% retention; ~10 pillars per bin; 11-bin mean; solid A, dashed B",
        ha="left",
        va="bottom",
        fontsize=6.0,
        color="#5B6570",
    )

    # Difference heatmap: positive values indicate greater B retention in the
    # same local window, making the A/B contrast visible at a glance.
    diff_matrix = np.array(
        [wave_values[(sp, "B")] - wave_values[(sp, "A")] for sp in SPECIES]
    )
    diff_limit = 0.25
    diff_norm = mpl.colors.TwoSlopeNorm(vmin=-diff_limit, vcenter=0.0, vmax=diff_limit)
    diff_image = ax_diff.imshow(
        diff_matrix,
        origin="upper",
        interpolation="nearest",
        aspect="auto",
        extent=[0, n, 0, len(SPECIES)],
        cmap="RdBu",
        norm=diff_norm,
    )
    ax_diff.set_xlim(0, n)
    ax_diff.set_yticks([2.5, 1.5, 0.5])
    ax_diff.set_yticklabels(SPECIES, fontsize=6.2)
    ax_diff.tick_params(axis="y", length=0, pad=2)
    ax_diff.set_xticks([])
    ax_diff.set_title(
        "B − A local retention (blue = B higher)",
        fontsize=7.0,
        loc="left",
        pad=6,
    )
    chrom_starts = np.cumsum(chrom_groups.to_numpy())[:-1]
    for start in chrom_starts:
        ax_diff.axvline(start, color="white", linewidth=0.4, alpha=0.8)
    cbar = fig.colorbar(
        diff_image,
        ax=ax_diff,
        fraction=0.08,
        pad=0.06,
        ticks=[-diff_limit, 0, diff_limit],
    )
    cbar.ax.set_yticklabels(["−25", "0", "+25"], fontsize=5.5)
    cbar.set_label("percentage points\n(clipped at ±25)", fontsize=5.8, labelpad=3)
    cbar.outline.set_linewidth(0.3)

    # Summary panel: paired A/B states and exact global A/B retention.
    categories = ["both", "a_only", "b_only", "none"]
    y_positions = np.arange(len(SPECIES))[::-1]
    left = np.zeros(len(SPECIES))
    for cat in categories:
        vals = np.array(
            [(data[f"{sp}_paired_state"] == cat).sum() / n for sp in SPECIES]
        )
        ax_bar.barh(
            y_positions,
            vals,
            left=left,
            height=0.52,
            color=COLORS[cat],
            edgecolor="white",
            linewidth=0.45,
        )
        for y, x0, value in zip(y_positions, left, vals):
            # The small one-sided-loss segments are encoded by color and are
            # reported exactly in the source-data summary; labeling only the
            # dominant retained state avoids collisions at print size.
            if value >= 0.045 and cat == "both":
                ax_bar.text(
                    x0 + value / 2,
                    y,
                    f"{value * 100:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=5.4,
                    color=(
                        "white"
                        if cat in {"both", "a_only", "b_only"}
                        else COLORS["outline"]
                    ),
                )
        left += vals

    # Leave room for labels in the small terminal segments.
    ax_bar.set_xlim(0, 1.43)
    ax_bar.set_ylim(-0.6, 2.6)
    ax_bar.set_yticks(y_positions)
    ax_bar.set_yticklabels(SPECIES, fontsize=6.8)
    ax_bar.set_xticks([0, 0.5, 1.0])
    ax_bar.set_xticklabels(["0", "50", "100"], fontsize=6.2)
    ax_bar.set_xlabel("Pillars (%)", fontsize=6.8, labelpad=5)
    ax_bar.set_title("Paired A/B state", fontsize=7.3, loc="left", pad=8)
    ax_bar.grid(axis="x", color="#D8D8D8", linewidth=0.45, zorder=0)
    ax_bar.set_axisbelow(True)
    ax_bar.spines["left"].set_visible(False)
    ax_bar.tick_params(axis="y", length=0, pad=3)
    for y, sp in zip(y_positions, SPECIES):
        a_retained = int(data[f"{sp}_A"].sum())
        b_retained = int(data[f"{sp}_B"].sum())
        ax_bar.text(
            1.04,
            y,
            f"A {a_retained / n * 100:.1f}%\nB {b_retained / n * 100:.1f}%",
            ha="left",
            va="center",
            fontsize=5.5,
            linespacing=1.2,
            color=COLORS["outline"],
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=COLORS["outline"],
            linewidth=1.0,
            linestyle="-",
            label="A track",
        ),
        Line2D(
            [0],
            [0],
            color=COLORS["outline"],
            linewidth=1.0,
            linestyle="--",
            label="B track",
        ),
        Patch(facecolor=COLORS["both"], label="A and B retained"),
        Patch(facecolor=COLORS["a_only"], label="A retained only"),
        Patch(facecolor=COLORS["b_only"], label="B retained only"),
        Patch(facecolor=COLORS["none"], label="Both lost"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.54, 0.985),
        ncol=6,
        fontsize=6.0,
        handlelength=1.0,
        handleheight=0.7,
        columnspacing=1.1,
        borderaxespad=0,
    )

    fig.text(
        0.010,
        0.985,
        "a",
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
        color="black",
    )

    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

    # Lightweight audit record for the manuscript/source-data bundle.
    summary = []
    for sp in SPECIES:
        counts = (
            data[f"{sp}_paired_state"].value_counts().reindex(categories, fill_value=0)
        )
        a_retained = int(data[f"{sp}_A"].sum())
        b_retained = int(data[f"{sp}_B"].sum())
        summary.append(
            {
                "species": sp,
                "n_pillars": n,
                "a_retained_count": a_retained,
                "b_retained_count": b_retained,
                "a_retained_fraction": float(a_retained / n),
                "b_retained_fraction": float(b_retained / n),
                **{f"{cat}_count": int(counts[cat]) for cat in categories},
                **{f"{cat}_fraction": float(counts[cat] / n) for cat in categories},
            }
        )
    pd.DataFrame(summary).to_csv(
        OUT_BASE.with_name(OUT_BASE.name + "_summary.tsv"), sep="\t", index=False
    )


if __name__ == "__main__":
    main()
