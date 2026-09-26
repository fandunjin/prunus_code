#!/usr/bin/env python3
"""Draw the requested Figure 2 Ks panels on smooth nonlinear Ks axes.

The density is always calculated on the original Ks scale. Only the display
coordinate is smoothly compressed using ``asinh(Ks / scale)``; no artificial
breaks or discontinuities are drawn. Panel c optionally applies one uniform
lower Ks cutoff to all self-syntenic datasets before density estimation.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import OrderedDict, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.75,
    "legend.frameon": False,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
})


DISPLAY = {
    "PZI": r"$P.$ zippeliana",
    "PSP": r"$P.$ spinulosa",
    "PPA": r"$P.$ padus",
    "PCAM": r"$P.$ campanulata",
    "PPE": r"$P.$ persica",
    "GTR": r"$G.$ trifoliata",
}
ABBR = {"PZI": "PZI", "PSP": "PSP", "PPA": "PPA", "PCAM": "PCAM", "PPE": "PPE", "GTR": "GTR"}
SPECIES_COLORS = {"PZI": "#B64643", "PSP": "#1E5A9E", "PPA": "#39836E"}
SELF_COLORS = {
    "PZI_A": "#B64643", "PZI_B": "#B64643",
    "PSP_A": "#1E5A9E", "PSP_B": "#1E5A9E",
    "PPA_A": "#39836E", "PPA_B": "#39836E",
    "PCAM": "#C17A27", "PPE": "#76539B", "GTR": "#404040",
}
LINE_STYLES = {"A": "-", "B": (0, (4.0, 2.1))}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Missing header in {}".format(path))
        return list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class SmoothKsAxis:
    """Continuously compress high Ks values without a discontinuous axis break."""

    def __init__(self, scale: float, xmax: float = 5.0):
        if scale <= 0.0 or xmax <= 0.0:
            raise ValueError("Smooth Ks axis scale and xmax must be positive")
        self.scale = float(scale)
        self.xmax = float(xmax)

    def forward(self, values: np.ndarray | float) -> np.ndarray | float:
        array = np.asarray(values, dtype=float)
        clipped = np.clip(array, 0.0, self.xmax)
        result = np.arcsinh(clipped / self.scale)
        if np.isscalar(values):
            return float(result)
        return result


class FocusedKsAxis(SmoothKsAxis):
    """Smoothly devote more display width to a scientifically relevant Ks window."""

    def __init__(self, scale: float, focus_center: float, focus_width: float,
                 focus_gain: float, xmax: float = 5.0):
        super().__init__(scale=scale, xmax=xmax)
        if focus_center <= 0.0 or focus_width <= 0.0 or focus_gain < 0.0:
            raise ValueError("Focused Ks axis parameters must be positive")
        self.focus_center = float(focus_center)
        self.focus_width = float(focus_width)
        self.focus_gain = float(focus_gain)

    def forward(self, values: np.ndarray | float) -> np.ndarray | float:
        array = np.asarray(values, dtype=float)
        clipped = np.clip(array, 0.0, self.xmax)
        base = np.arcsinh(clipped / self.scale)
        zero_offset = math.erf(-self.focus_center / self.focus_width)
        focus = self.focus_gain * (
            np.vectorize(math.erf)((clipped - self.focus_center) / self.focus_width) - zero_offset
        )
        result = base + focus
        if np.isscalar(values):
            return float(result)
        return result


def gaussian_kde(values: np.ndarray, grid: np.ndarray, bandwidth: float,
                 reflect_at_zero: bool = True) -> np.ndarray:
    density = np.zeros_like(grid, dtype=float)
    denominator = len(values) * bandwidth * math.sqrt(2.0 * math.pi)
    for start in range(0, len(values), 4096):
        block = values[start:start + 4096, None]
        direct = (grid[None, :] - block) / bandwidth
        density += np.exp(-0.5 * direct * direct).sum(axis=0)
        if reflect_at_zero:
            reflected = (grid[None, :] + block) / bandwidth
            density += np.exp(-0.5 * reflected * reflected).sum(axis=0)
    return density / denominator


def read_grouped(path: Path, expected: tuple[str, ...], min_ks: float = 0.0) -> tuple[dict[str, list[float]], list[dict[str, object]]]:
    if min_ks < 0.0 or min_ks >= 5.0:
        raise ValueError("min_ks must satisfy 0<=min_ks<5")
    result = {name: [] for name in expected}
    before = {name: 0 for name in expected}
    for row in read_tsv(path):
        name = row.get("comparison", "")
        if name not in result:
            continue
        try:
            value = float(row["Ks"])
        except (KeyError, TypeError, ValueError):
            raise ValueError("Invalid Ks row in {}".format(path))
        if not math.isfinite(value) or not 0.0 < value < 5.0:
            raise ValueError("Unfiltered Ks in {}: {}".format(path, value))
        before[name] += 1
        if value < min_ks:
            continue
        result[name].append(value)
    missing = [name for name, values in result.items() if not values]
    if missing:
        raise ValueError("Missing or empty comparisons in {}: {}".format(path, ", ".join(missing)))
    audit_rows = [
        {
            "comparison": name,
            "pairs_before_minimum": before[name],
            "pairs_removed_below_minimum": before[name] - len(values),
            "pairs_retained": len(values),
            "minimum_Ks": "{:.4f}".format(min_ks),
            "retained_percent": "{:.4f}".format(100.0 * len(values) / before[name]),
        }
        for name, values in result.items()
    ]
    return result, audit_rows


def true_tick_label(value: float) -> str:
    if value == 0:
        return "0"
    if value < 0.1:
        return "{:.3f}".format(value).rstrip("0").rstrip(".")
    if value < 1:
        return "{:.2f}".format(value).rstrip("0").rstrip(".")
    return "{:.1f}".format(value).rstrip("0").rstrip(".")


def decorate_smooth_axis(ax: plt.Axes, axis: SmoothKsAxis, ticks: tuple[float, ...]) -> None:
    transformed_ticks = [axis.forward(value) for value in ticks]
    ax.set_xlim(0.0, axis.forward(axis.xmax))
    ax.set_xticks(transformed_ticks)
    ax.set_xticklabels([true_tick_label(value) for value in ticks])
    ax.tick_params(axis="x", direction="out", length=2.5, pad=1.5, labelsize=4.8)
    ax.set_xlabel(r"$K_s$", labelpad=1.5, fontsize=6)


def panel_curves(grouped: dict[str, list[float]], order: tuple[str, ...], axis: SmoothKsAxis,
                 bandwidth: float, grid_min: float = 0.0,
                 reflect_at_zero: bool = True) -> tuple[np.ndarray, list[tuple[str, np.ndarray, int]]]:
    # 0.00167 Ks resolution across the full 0-5 calculation range is enough
    # for a figure annotated to three decimals while keeping the 18-curve
    # orthologue panel practical on a local workstation.
    if not 0.0 <= grid_min < axis.xmax:
        raise ValueError("grid_min must satisfy 0<=grid_min<xmax")
    grid = np.linspace(grid_min, axis.xmax, 3001)
    curves: list[tuple[str, np.ndarray, int]] = []
    for name in order:
        density = gaussian_kde(
            np.asarray(grouped[name], dtype=float), grid, bandwidth,
            reflect_at_zero=reflect_at_zero,
        )
        curves.append((name, density, int(np.argmax(density))))
    return grid, curves


def annotate_peak_rack(ax: plt.Axes, axis: SmoothKsAxis, grid: np.ndarray,
                       curves: list[tuple[str, np.ndarray, int]],
                       slots: dict[str, tuple[float, float]], peak_rows: list[dict[str, object]],
                       panel: str, facet: str = "") -> None:
    """Use fixed, separated label slots so closely spaced modes remain legible."""
    for name, density, peak_index in curves:
        true_mode = float(grid[peak_index])
        x_mode = axis.forward(true_mode)
        try:
            x_fraction, y_fraction = slots[name]
        except KeyError as error:
            raise ValueError("Missing annotation slot for {}".format(name)) from error
        ax.annotate(
            r"${:.3f}$".format(true_mode),
            xy=(x_mode, float(density[peak_index])),
            xycoords="data", xytext=(x_fraction, y_fraction), textcoords="axes fraction",
            ha="center", va="center", fontsize=4.4,
            color=peak_color(name),
            arrowprops={"arrowstyle": "->", "lw": 0.38, "color": "#555555", "shrinkA": 1.2, "shrinkB": 1.2},
            annotation_clip=False,
        )
        peak_rows.append({
            "panel": panel, "facet": facet, "comparison": name,
            "KDE_mode_Ks": "{:.8f}".format(true_mode),
            "mode_density": "{:.8f}".format(float(density[peak_index])),
            "display_x": "{:.8f}".format(float(x_mode)),
            "annotation_display_x": "{:.8f}".format(float(x_fraction)),
            "annotation_display_y": "{:.8f}".format(float(y_fraction)),
            "annotation_layout": "fixed_nonoverlapping_rack",
        })


def peak_color(name: str) -> str:
    """Return the line color associated with one panel comparison."""
    species = name.split("_", 1)[0]
    if species in SPECIES_COLORS:
        return SPECIES_COLORS[species]
    return SELF_COLORS.get(name, "#202020")


def plot_a(ax: plt.Axes, grouped: dict[str, list[float]], peak_rows: list[dict[str, object]],
           density_rows: list[dict[str, object]]) -> None:
    order = ("PZI", "PSP", "PPA")
    axis = SmoothKsAxis(scale=0.06)
    grid, curves = panel_curves(grouped, order, axis, 0.035)
    transformed_grid = axis.forward(grid)
    ymax = max(float(density.max()) for _, density, _ in curves) * 1.30
    for name, density, _ in curves:
        color = SPECIES_COLORS[name]
        ax.fill_between(transformed_grid, density, color=color, alpha=0.065, linewidth=0)
        ax.plot(transformed_grid, density, color=color, lw=1.25, label=ABBR[name])
        for true_x, display_x, y in zip(grid, transformed_grid, density):
            density_rows.append({"panel": "a", "facet": "", "comparison": name,
                                 "Ks_grid": "{:.8f}".format(float(true_x)), "display_x": "{:.8f}".format(float(display_x)),
                                 "density": "{:.8f}".format(float(y)), "n": len(grouped[name]), "bandwidth": "0.0350"})
    ax.set_ylabel("Paralogue density")
    ax.set_ylim(0.0, ymax)
    annotate_peak_rack(
        ax, axis, grid, curves,
        {"PPA": (0.25, 0.90), "PZI": (0.43, 0.82), "PSP": (0.61, 0.96)},
        peak_rows, "a",
    )
    decorate_smooth_axis(ax, axis, (0.0, 0.05, 0.10, 0.20, 1.0, 5.0))
    ax.legend(loc="upper right", fontsize=4.1, handlelength=1.3, labelspacing=0.18)
    ax.text(-0.105, 1.015, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def plot_b(axes: list[plt.Axes], grouped: dict[str, list[float]], peak_rows: list[dict[str, object]],
           density_rows: list[dict[str, object]]) -> None:
    order = ("PZI", "PSP", "PPA")
    per_reference = {
        "PCAM": (SmoothKsAxis(scale=0.04), (0.0, 0.04, 0.08, 0.12, 0.16, 0.30, 0.50, 1.0, 2.0, 5.0)),
        "PPE": (SmoothKsAxis(scale=0.04), (0.0, 0.04, 0.08, 0.12, 0.16, 0.30, 0.50, 1.0, 2.0, 5.0)),
        "GTR": (
            FocusedKsAxis(scale=0.12, focus_center=0.35, focus_width=0.12, focus_gain=1.00),
            (0.0, 0.24, 0.30, 0.36, 0.44, 0.52, 1.0, 5.0),
        ),
    }
    for index, (reference, ax) in enumerate(zip(("PCAM", "PPE", "GTR"), axes)):
        axis, ticks = per_reference[reference]
        curve_order = tuple("{}_{}_{}".format(species, subgenome, reference) for species in order for subgenome in ("A", "B"))
        grid, curves = panel_curves(grouped, curve_order, axis, 0.040)
        transformed_grid = axis.forward(grid)
        ymax = max(float(density.max()) for _, density, _ in curves) * 1.36
        for name, density, _ in curves:
            species, subgenome, _ = name.split("_", 2)
            color = SPECIES_COLORS[species]
            ax.plot(transformed_grid, density, color=color, lw=1.02, ls=LINE_STYLES[subgenome],
                    label="{} {}".format(ABBR[species], subgenome))
            for true_x, display_x, y in zip(grid, transformed_grid, density):
                density_rows.append({"panel": "b", "facet": reference, "comparison": name,
                                     "Ks_grid": "{:.8f}".format(float(true_x)), "display_x": "{:.8f}".format(float(display_x)),
                                     "density": "{:.8f}".format(float(y)), "n": len(grouped[name]), "bandwidth": "0.0400"})
        ax.set_ylabel("Orthologue density")
        ax.set_ylim(0.0, ymax)
        annotate_peak_rack(
            ax, axis, grid, curves,
            {
                "PZI_A_" + reference: (0.30, 0.96), "PZI_B_" + reference: (0.45, 0.96),
                "PSP_A_" + reference: (0.60, 0.96), "PSP_B_" + reference: (0.30, 0.84),
                "PPA_A_" + reference: (0.45, 0.84), "PPA_B_" + reference: (0.60, 0.84),
            }, peak_rows, "b", reference,
        )
        decorate_smooth_axis(ax, axis, ticks)
        ax.text(0.012, 0.94, ABBR[reference], transform=ax.transAxes, fontsize=5.2, va="top")
        if index == 0:
            ax.legend(loc="upper right", fontsize=3.8, ncol=2, handlelength=1.25, columnspacing=0.35, labelspacing=0.15)
            ax.text(-0.105, 1.015, "b", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def plot_c(ax: plt.Axes, grouped: dict[str, list[float]], peak_rows: list[dict[str, object]],
           density_rows: list[dict[str, object]]) -> None:
    order = ("PZI_A", "PSP_A", "PPA_A", "PZI_B", "PSP_B", "PPA_B", "PCAM", "PPE", "GTR")
    axis = SmoothKsAxis(scale=0.10)
    grid, curves = panel_curves(
        grouped, order, axis, 0.035,
        grid_min=min(min(grouped[name]) for name in order),
        reflect_at_zero=False,
    )
    transformed_grid = axis.forward(grid)
    ymax = max(float(density.max()) for _, density, _ in curves) * 1.42
    for name, density, _ in curves:
        if name.endswith("_A"):
            line_style = LINE_STYLES["A"]
            label = "{} A".format(ABBR[name.split("_", 1)[0]])
        elif name.endswith("_B"):
            line_style = LINE_STYLES["B"]
            label = "{} B".format(ABBR[name.split("_", 1)[0]])
        else:
            line_style = "-"
            label = ABBR[name]
        color = SELF_COLORS[name]
        ax.plot(transformed_grid, density, color=color, lw=1.0, ls=line_style, label=label)
        for true_x, display_x, y in zip(grid, transformed_grid, density):
            density_rows.append({"panel": "c", "facet": "", "comparison": name,
                                 "Ks_grid": "{:.8f}".format(float(true_x)), "display_x": "{:.8f}".format(float(display_x)),
                                 "density": "{:.8f}".format(float(y)), "n": len(grouped[name]), "bandwidth": "0.0350"})
    ax.set_ylabel("Paralogue density")
    ax.set_ylim(0.0, ymax)
    annotate_peak_rack(
        ax, axis, grid, curves,
        {
            "PZI_A": (0.17, 0.96), "PSP_A": (0.32, 0.96), "PPA_A": (0.47, 0.96), "PZI_B": (0.62, 0.96),
            "PSP_B": (0.17, 0.84), "PPA_B": (0.32, 0.84), "PCAM": (0.47, 0.84), "PPE": (0.62, 0.84),
            "GTR": (0.84, 0.62),
        }, peak_rows, "c",
    )
    decorate_smooth_axis(ax, axis, (0.0, 0.05, 0.10, 0.20, 0.75, 1.5, 5.0))
    ax.legend(loc="upper right", fontsize=3.55, ncol=2, handlelength=1.15, columnspacing=0.28, labelspacing=0.13)
    ax.text(-0.105, 1.015, "c", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-a", type=Path, required=True)
    parser.add_argument("--panel-b", type=Path, required=True)
    parser.add_argument("--panel-c", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--panel-c-min-ks", type=float, default=0.02,
                        help="Uniform lower Ks cutoff for every panel c self-synteny curve")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.overwrite:
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.panel_c_min_ks < 0.0 or args.panel_c_min_ks >= 5.0:
        raise SystemExit("[ERROR] --panel-c-min-ks must satisfy 0<=value<5")
    a, _ = read_grouped(args.panel_a, ("PZI", "PSP", "PPA"))
    b_order = tuple("{}_{}_{}".format(species, subgenome, reference)
                    for species in ("PZI", "PSP", "PPA") for subgenome in ("A", "B")
                    for reference in ("PCAM", "PPE", "GTR"))
    b, _ = read_grouped(args.panel_b, b_order)
    c, c_filter_rows = read_grouped(
        args.panel_c,
        ("PZI_A", "PSP_A", "PPA_A", "PZI_B", "PSP_B", "PPA_B", "PCAM", "PPE", "GTR"),
        min_ks=args.panel_c_min_ks,
    )

    fig = plt.figure(figsize=(1.80, 11.25))
    grid_spec = fig.add_gridspec(5, 1, height_ratios=(1.65, 1.35, 1.35, 1.35, 2.05), hspace=0.72)
    ax_a = fig.add_subplot(grid_spec[0])
    axes_b = [fig.add_subplot(grid_spec[index]) for index in (1, 2, 3)]
    ax_c = fig.add_subplot(grid_spec[4])
    density_rows: list[dict[str, object]] = []
    peak_rows: list[dict[str, object]] = []
    plot_a(ax_a, a, peak_rows, density_rows)
    plot_b(axes_b, b, peak_rows, density_rows)
    plot_c(ax_c, c, peak_rows, density_rows)
    for axis in [ax_a] + axes_b + [ax_c]:
        axis.tick_params(axis="y", direction="out", length=3, pad=2)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.985, bottom=0.055)

    prefix = args.outdir / "figure2_ks_smooth_axis"
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)

    write_tsv(args.outdir / "density_grid.tsv",
              ["panel", "facet", "comparison", "Ks_grid", "display_x", "density", "n", "bandwidth"], density_rows)
    write_tsv(args.outdir / "peak_labels.tsv",
              ["panel", "facet", "comparison", "KDE_mode_Ks", "mode_density", "display_x", "annotation_display_x", "annotation_display_y", "annotation_layout"], peak_rows)
    write_tsv(args.outdir / "panel_c_low_ks_filter_audit.tsv",
              ["comparison", "pairs_before_minimum", "pairs_removed_below_minimum", "pairs_retained", "minimum_Ks", "retained_percent"],
              c_filter_rows)
    (args.outdir / "plot_parameters.tsv").write_text(
        "backend\tPython matplotlib\n"
        "density\tGaussian KDE computed on original Ks scale; reflection at Ks=0 for panels a/b, no reflection after panel-c lower truncation\n"
        "display_axis\tcontinuous asinh(Ks/scale) coordinate with true Ks tick labels; no discontinuous breaks\n"
        "panel_a\tCodeML F3X4 A-B collinear homoeologues; bandwidth=0.035; smooth_scale=0.06\n"
        "panel_b\tYN A/B-to-reference collinear orthologues; bandwidth=0.040; smooth_scale=0.04 for PCAM/PPE and GTR asinh scale=0.12 with a continuous focus window centered at Ks=0.35\n"
        "panel_c\tCodeML cross-chromosome self-syntenic paralogues; bandwidth=0.035; smooth_scale=0.10; uniform lower cutoff Ks>={:.4f}; density drawn only over the retained range\n"
        "filter\tAll source rows satisfy 0<Ks<5 before the panel-c lower cutoff\n"
        "peak_annotation\tKDE modes; labels placed in predefined non-overlapping panel racks and connected to the peak by arrows\n"
        "axis_notice\tcontinuous nonlinear display coordinate; reduced true Ks tick labels retained; horizontal axis title is Ks\n".format(args.panel_c_min_ks),
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Smooth-axis three-panel Ks figure plotted.\n", encoding="ascii")
    print("[OK] Wrote smooth-axis three-panel Ks figure: {}".format(args.outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
