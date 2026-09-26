#!/usr/bin/env python3
"""Aggregate and plot five corrected, phase-aware Figure 2 Ks distributions.

The script only consumes the newly generated Ks tables below the corrected
``figure2_repartition_ks`` root.  It does not read any earlier Figure 2 A/B
tables, anchors, or figures.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PANEL_ORDER = ["orthologue_pcam", "orthologue_ppe", "orthologue_gtr", "ab_paralogue", "self_paralogue"]
PANEL_META = {
    "orthologue_pcam": ("a", "PCAM", "Orthologue density"),
    "orthologue_ppe": ("b", "PPE", "Orthologue density"),
    "orthologue_gtr": ("c", "GTR", "Orthologue density"),
    "ab_paralogue": ("d", "A-B collinear pairs", "Paralogue density"),
    "self_paralogue": ("e", "Self-syntenic pairs", "Paralogue density"),
}
COLORS = {
    "PZI": "#B45657", "PSP": "#4D78A8", "PPA": "#3D8873",
    "PCAM": "#8A659D", "PPE": "#B77935", "GTR": "#545454",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--codeml-source", type=Path,
                        help="Combined CodeML source table from combine_fig2_repartition_codeml_results.py")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--xmax", type=float, default=5.0,
                        help="Full Ks calculation/display range; the display is continuously compressed")
    parser.add_argument("--self-plot-xmax", type=float, default=3.0,
                        help="Upper Ks limit used for the self-synteny panel; values above this are saturation-prone")
    parser.add_argument("--self-axis-xmax", type=float, default=5.0,
                        help="KDE display limit for the self-synteny panel; extends beyond the data cutoff for a smooth tail")
    parser.add_argument("--self-min-ks", type=float, default=0.02,
                        help="Uniform lower cutoff for all self-synteny curves, matching the old panel-c display")
    parser.add_argument("--points", type=int, default=2400,
                        help="Number of Ks grid points used for each density curve")
    parser.add_argument("--min-curve-pairs", type=int, default=50)
    parser.add_argument("--bandwidth", type=float, default=0.035)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError("No header: {}".format(path))
        return list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def valid_ks(row: dict[str, str]) -> float | None:
    try:
        value = float(row.get("Ks", ""))
    except ValueError:
        return None
    if not math.isfinite(value) or not (0.0 < value < 5.0):
        return None
    return value


def plan_lookup(path: Path, include_panel: set[str]) -> dict[str, dict[str, str]]:
    values: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        if row.get("panel") not in include_panel:
            continue
        key = row.get("comparison", "")
        if not key or key in values:
            raise RuntimeError("Invalid or duplicate planned comparison in {}: {}".format(path, key))
        values[key] = row
    return values


def self_plan_lookup(path: Path) -> dict[str, dict[str, str]]:
    values: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        label = row.get("label", "")
        if not label or label in values:
            raise RuntimeError("Invalid or duplicate self plan label: {}".format(label))
        values[label + "_self"] = row
    return values


def estimate_bandwidth(values: np.ndarray, minimum: float) -> float:
    if len(values) < 2:
        return minimum
    sd = float(np.std(values, ddof=1))
    iqr = float(np.quantile(values, 0.75) - np.quantile(values, 0.25))
    scale = min(sd, iqr / 1.34) if iqr > 0.0 else sd
    if scale <= 0.0 or not math.isfinite(scale):
        return minimum
    return max(minimum, 0.9 * scale * len(values) ** (-0.2))


def reflected_kde(values: np.ndarray, grid: np.ndarray, minimum_bandwidth: float,
                  reflect_at_zero: bool = True) -> tuple[np.ndarray, float]:
    # The old panel-c figure used a fixed 0.035 Ks bandwidth.  Use exactly the
    # same fixed bandwidth throughout these five panels so curve height and
    # apparent peak width are not driven by sample-size-specific smoothing.
    bandwidth = minimum_bandwidth
    density = np.zeros_like(grid)
    norm = 1.0 / (math.sqrt(2.0 * math.pi) * bandwidth * len(values))
    # Chunking keeps memory bounded for curves with many syntenic pairs.
    for start in range(0, len(values), 2000):
        chunk = values[start:start + 2000]
        direct = (grid[:, None] - chunk[None, :]) / bandwidth
        density += np.exp(-0.5 * direct * direct).sum(axis=1)
        if reflect_at_zero:
            reflected = (grid[:, None] + chunk[None, :]) / bandwidth
            density += np.exp(-0.5 * reflected * reflected).sum(axis=1)
    return density * norm, bandwidth


def curve_style(panel: str, curve: str) -> tuple[str, str, float]:
    if curve.startswith("PZI"):
        color = COLORS["PZI"]
    elif curve.startswith("PSP"):
        color = COLORS["PSP"]
    elif curve.startswith("PPA"):
        color = COLORS["PPA"]
    else:
        color = COLORS[curve]
    linestyle = "--" if curve.endswith("_B") else "-"
    return color, linestyle, 1.15 if panel != "self_paralogue" else 1.0


def sort_curves(panel: str, curves: set[str]) -> list[str]:
    order = {
        "orthologue_pcam": ["PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B"],
        "orthologue_ppe": ["PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B"],
        "orthologue_gtr": ["PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B"],
        "ab_paralogue": ["PZI", "PSP", "PPA"],
        "self_paralogue": ["PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE", "GTR"],
    }[panel]
    return [item for item in order if item in curves]


class SmoothKsAxis:
    """Continuous nonlinear display coordinate with true Ks tick labels."""

    def __init__(self, scale: float, xmax: float = 5.0):
        if scale <= 0.0 or xmax <= 0.0:
            raise ValueError("axis scale and xmax must be positive")
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
    """Smoothly allocates extra display width around a target Ks region."""

    def __init__(self, scale: float, focus_center: float, focus_width: float,
                 focus_gain: float, xmax: float = 5.0):
        super().__init__(scale, xmax)
        if focus_center <= 0.0 or focus_width <= 0.0 or focus_gain < 0.0:
            raise ValueError("invalid focus parameters")
        self.focus_center = float(focus_center)
        self.focus_width = float(focus_width)
        self.focus_gain = float(focus_gain)

    def forward(self, values: np.ndarray | float) -> np.ndarray | float:
        array = np.asarray(values, dtype=float)
        clipped = np.clip(array, 0.0, self.xmax)
        base = np.arcsinh(clipped / self.scale)
        zero_offset = math.erf(-self.focus_center / self.focus_width)
        erf_values = np.vectorize(math.erf)((clipped - self.focus_center) / self.focus_width)
        result = base + self.focus_gain * (erf_values - zero_offset)
        if np.isscalar(values):
            return float(result)
        return result


def true_tick_label(value: float) -> str:
    if value == 0.0:
        return "0"
    if value < 0.1:
        return "{:.3f}".format(value).rstrip("0").rstrip(".")
    if value < 1.0:
        return "{:.2f}".format(value).rstrip("0").rstrip(".")
    return "{:.1f}".format(value).rstrip("0").rstrip(".")


def decorate_axis(ax: plt.Axes, axis: SmoothKsAxis, ticks: tuple[float, ...], xlabel: bool = True) -> None:
    ax.set_xlim(0.0, float(axis.forward(axis.xmax)))
    ax.set_xticks([float(axis.forward(value)) for value in ticks])
    ax.set_xticklabels([true_tick_label(value) for value in ticks])
    ax.tick_params(axis="x", direction="out", length=2.5, width=0.7, pad=1.5, labelsize=5.2)
    if xlabel:
        ax.set_xlabel("Ks", labelpad=1.5, fontsize=6.5)


def ks_grid(points: int, xmax: float) -> np.ndarray:
    """Use dense grids at both low-Ks and self-synteny peak regions."""
    if points < 400:
        raise ValueError("points must be >=400")
    weights = (0.34, 0.10, 0.42, 0.14)
    counts = [max(20, int(points * weight)) for weight in weights]
    # Clip every breakpoint to the requested range.  Without this guard, a
    # self-synteny axis ending below 4.4 would create grid points beyond xmax;
    # the coordinate transform would then collapse them onto the right border
    # and produce an artificial vertical tail.
    breakpoints = [0.0, min(0.6, xmax), min(2.5, xmax), min(4.4, xmax), xmax]
    edges = [(start, stop) for start, stop in zip(breakpoints, breakpoints[1:]) if stop > start]
    parts = []
    for index, (start, stop) in enumerate(edges):
        count = counts[min(index, len(counts) - 1)]
        parts.append(np.linspace(start, stop, count, endpoint=False))
    parts.append(np.asarray([xmax], dtype=float))
    return np.unique(np.concatenate(parts))


def panel_axis(panel: str, xmax: float) -> tuple[SmoothKsAxis, tuple[float, ...], float, bool]:
    if panel in {"orthologue_pcam", "orthologue_ppe", "ab_paralogue"}:
        # Keep the low-Ks region readable while retaining the long tail.
        return SmoothKsAxis(scale=0.060, xmax=xmax), (0.0, 0.04, 0.08, 0.12, 0.20, 0.30, 0.50, 1.0, 2.0, 3.0, 5.0), 0.035, True
    if panel == "orthologue_gtr":
        return FocusedKsAxis(scale=0.12, focus_center=0.30, focus_width=0.12, focus_gain=1.0, xmax=xmax), (0.0, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.50, 1.0, 2.0, 3.0, 5.0), 0.035, True
    return FocusedKsAxis(scale=0.10, focus_center=1.65, focus_width=0.50, focus_gain=1.10, xmax=xmax), (0.0, 0.05, 0.10, 0.20, 0.50, 1.0, 1.5, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0), 0.035, False


def label_for(panel: str, curve: str) -> str:
    if panel.startswith("orthologue_"):
        return "{} {}".format(curve.split("_", 1)[0], curve.split("_")[1])
    if panel == "ab_paralogue":
        return curve
    if "_" in curve:
        return "{} {}".format(*curve.split("_", 1))
    return curve


def annotation_slots(panel: str, curves: list[str]) -> dict[str, tuple[float, float]]:
    if panel in {"orthologue_pcam", "orthologue_ppe"}:
        # Two compact rows leave the upper-right legend clear.
        xs = (0.16, 0.33, 0.50, 0.16, 0.33, 0.50)
        ys = (0.95, 0.95, 0.95, 0.82, 0.82, 0.82)
        return {name: (xs[index], ys[index]) for index, name in enumerate(curves)}
    if panel == "orthologue_gtr":
        xs = (0.16, 0.33, 0.50, 0.16, 0.33, 0.50)
        ys = (0.95, 0.95, 0.95, 0.82, 0.82, 0.82)
        return {name: (xs[index], ys[index]) for index, name in enumerate(curves)}
    if panel == "ab_paralogue":
        xs = (0.24, 0.43, 0.62)
        return {name: (xs[index], 0.93) for index, name in enumerate(curves)}
    # Self-synteny modes are evaluated below the saturation-prone Ks>3 tail.
    # Keep labels inside the axes in two compact racks, while arrows point back
    # to the corresponding internal modes.
    slots = {}
    xs_top = (0.10, 0.23, 0.36, 0.49, 0.62)
    xs_bottom = (0.16, 0.32, 0.48, 0.64)
    for index, name in enumerate(curves):
        if index < 5:
            slots[name] = (xs_top[index], 0.96)
        else:
            slots[name] = (xs_bottom[index - 5], 0.82)
    return slots


def annotate_modes(ax: plt.Axes, axis: SmoothKsAxis, grid: np.ndarray,
                   curves: list[tuple[str, np.ndarray, int]], panel: str,
                   peak_rows: list[dict[str, object]]) -> None:
    # Assign text slots by increasing peak position so arrows remain ordered
    # even when several curves have nearly identical modes (notably GTR).
    ordered = sorted(curves, key=lambda item: float(grid[item[2]]))
    slots = annotation_slots(panel, [name for name, _, _ in ordered])
    for name, density, peak_index in curves:
        mode_ks = float(grid[peak_index])
        x_mode = float(axis.forward(mode_ks))
        x_fraction, y_fraction = slots[name]
        color, _, _ = curve_style(panel, name)
        ax.annotate(
            "{:.3f}".format(mode_ks), xy=(x_mode, float(density[peak_index])),
            xycoords="data", xytext=(x_fraction, y_fraction), textcoords="axes fraction",
            ha="center", va="center", fontsize=4.2, color=color,
            arrowprops={"arrowstyle": "->", "lw": 0.35, "color": "#555555", "shrinkA": 1.0, "shrinkB": 1.0},
            annotation_clip=False,
        )
        peak_rows.append({
            "panel": panel, "curve": name, "peak_Ks": "{:.8g}".format(mode_ks),
            "peak_density": "{:.8g}".format(float(density[peak_index])),
        })


def collect(base: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    jcvi_plan = plan_lookup(base / "plans" / "jcvi_pair_plan.tsv", set(PANEL_ORDER[:4]))
    self_plan = self_plan_lookup(base / "plans" / "self_pair_plan.tsv")
    tables = [
        ("jcvi", base / "03_yn_kaks" / "jcvi" / "summary" / "all_merged_4dtv_kaks.tsv", jcvi_plan),
        ("self", base / "03_yn_kaks" / "self" / "summary" / "all_merged_4dtv_kaks.tsv", self_plan),
    ]
    data: list[dict[str, object]] = []
    audits: list[dict[str, object]] = []
    for source_class, table, lookup in tables:
        if not table.is_file() or table.stat().st_size == 0:
            raise RuntimeError("Missing corrected Ks table: {}".format(table))
        raw = read_tsv(table)
        seen_comparisons: set[str] = set()
        per_comparison: defaultdict[str, dict[str, int]] = defaultdict(lambda: {"rows": 0, "valid": 0, "outside_xmax": 0})
        for row in raw:
            comparison = row.get("comparison", "")
            if comparison not in lookup:
                raise RuntimeError("Unexpected comparison in corrected {} Ks table: {}".format(source_class, comparison))
            seen_comparisons.add(comparison)
            per_comparison[comparison]["rows"] += 1
            ks = valid_ks(row)
            if ks is None:
                continue
            per_comparison[comparison]["valid"] += 1
            plan = lookup[comparison]
            panel = plan["panel"]
            data.append({
                "source_class": source_class,
                "panel": panel,
                "curve": plan["curve"],
                "comparison": comparison,
                "pair_id": row.get("pair_id", ""),
                "query_label": row.get("query_label", ""),
                "ref_label": row.get("ref_label", ""),
                "query_gene": row.get("query_gene", ""),
                "ref_gene": row.get("ref_gene", ""),
                "query_seqid": row.get("query_seqid", ""),
                "ref_seqid": row.get("ref_seqid", ""),
                "Ks": "{:.10g}".format(ks),
            })
        missing = sorted(set(lookup) - seen_comparisons)
        if missing:
            raise RuntimeError("No KaKs rows for planned corrected comparisons: {}".format(",".join(missing[:8])))
        for comparison in sorted(lookup):
            plan = lookup[comparison]
            counts = per_comparison[comparison]
            audits.append({
                "source_class": source_class, "panel": plan["panel"], "curve": plan["curve"],
                "comparison": comparison, "kaks_rows": counts["rows"], "valid_0_lt_Ks_lt_5": counts["valid"],
            })
    return data, audits


def collect_codeml(source: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Read only independently estimated CodeML Ks values for the five panels."""
    table = source / "five_panel_ks.tsv"
    audit = source / "comparison_audit.tsv"
    if not table.is_file() or table.stat().st_size == 0:
        raise RuntimeError("Missing combined CodeML Ks table: {}".format(table))
    if not audit.is_file() or audit.stat().st_size == 0:
        raise RuntimeError("Missing CodeML comparison audit: {}".format(audit))
    allowed = set(PANEL_ORDER)
    data = []
    for row in read_tsv(table):
        panel = row.get("panel", "")
        curve = row.get("curve", "")
        if panel not in allowed or not curve:
            raise RuntimeError("Invalid CodeML plotting metadata in {}".format(table))
        value = valid_ks(row)
        if value is None:
            raise RuntimeError("CodeML source contains a value outside 0<Ks<5: {}".format(table))
        data.append({
            "source_class": "codeml_pairwise",
            "panel": panel,
            "curve": curve,
            "comparison": row.get("source_comparison", ""),
            "pair_id": row.get("family_id", ""),
            "query_label": row.get("query_label", ""),
            "ref_label": row.get("ref_label", ""),
            "query_gene": row.get("gene_1", ""),
            "ref_gene": row.get("gene_2", ""),
            "query_seqid": "",
            "ref_seqid": "",
            "Ks": "{:.10g}".format(value),
        })
    if not data:
        raise RuntimeError("No CodeML Ks rows in {}".format(table))
    audits = []
    for row in read_tsv(audit):
        panel = row.get("panel", "")
        curve = row.get("curve", "")
        if panel not in allowed or not curve:
            raise RuntimeError("Invalid CodeML comparison audit in {}".format(audit))
        audits.append({
            "source_class": "codeml_pairwise",
            "panel": panel,
            "curve": curve,
            "comparison": row.get("source_comparison", ""),
            "codeml_rows_0_lt_Ks_lt_5": row.get("retained_0_lt_Ks_lt_5", ""),
        })
    return data, audits


def main() -> int:
    args = parse_args()
    if (args.xmax < 5.0 or args.self_plot_xmax <= 0.0 or args.self_plot_xmax > args.xmax
            or args.self_axis_xmax < args.self_plot_xmax or args.self_axis_xmax > args.xmax
            or args.self_min_ks < 0.0 or args.self_min_ks >= args.self_plot_xmax
            or args.points < 400 or args.bandwidth <= 0):
        raise RuntimeError("Invalid Ks range, self lower cutoff, grid size, or bandwidth")
    base = args.base.resolve()
    outdir = args.outdir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        if not args.force:
            raise RuntimeError("Output directory is not empty: {}".format(outdir))
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    source_dir = outdir / "source_data"
    source_dir.mkdir()

    if args.codeml_source is None:
        data, audits = collect(base)
        estimator = "MAFFT protein alignment, PAL2NAL codon backtranslation, KaKs_Calculator 3.0 YN"
    else:
        data, audits = collect_codeml(args.codeml_source.resolve())
        estimator = "MUSCLE 3.8.1551 protein alignment, PAL2NAL codon backtranslation, PAML CodeML pairwise F3X4"
    write_tsv(source_dir / "ks_plot_source.tsv", data, list(data[0]))
    write_tsv(source_dir / "comparison_audit.tsv", audits, list(audits[0]))
    grouped: defaultdict[tuple[str, str], list[float]] = defaultdict(list)
    for row in data:
        grouped[(str(row["panel"]), str(row["curve"]))].append(float(str(row["Ks"])))
    for panel in PANEL_ORDER:
        present = {curve for current_panel, curve in grouped if current_panel == panel}
        if not present:
            raise RuntimeError("No valid corrected Ks values for panel {}".format(panel))
        for curve in present:
            data_max = args.self_plot_xmax if panel == "self_paralogue" else args.xmax
            lower = args.self_min_ks if panel == "self_paralogue" else 0.0
            count = sum(lower <= value <= data_max for value in grouped[(panel, curve)])
            if count < args.min_curve_pairs:
                raise RuntimeError("{} {} has only {} valid Ks pairs (< {})".format(
                    panel, curve, count, args.min_curve_pairs))

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 6.5,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "legend.frameon": False,
        "savefig.facecolor": "white",
    })
    curve_rows: list[dict[str, object]] = []
    peak_rows: list[dict[str, object]] = []
    density_rows: list[dict[str, object]] = []
    panel_densities: dict[str, dict[str, object]] = {}

    # The reference-style asymmetric layout has three external comparisons on
    # the left, the A-B comparison at upper right, and self-synteny below it.
    fig = plt.figure(figsize=(5.95, 5.65), facecolor="white")
    outer = fig.add_gridspec(3, 2, width_ratios=(1.0, 1.03), height_ratios=(1, 1, 1),
                             wspace=0.46, hspace=0.76)
    axes_by_panel = {
        "orthologue_pcam": fig.add_subplot(outer[0, 0]),
        "orthologue_ppe": fig.add_subplot(outer[1, 0]),
        "orthologue_gtr": fig.add_subplot(outer[2, 0]),
        "ab_paralogue": fig.add_subplot(outer[0, 1]),
        "self_paralogue": fig.add_subplot(outer[1:, 1]),
    }

    for panel in PANEL_ORDER:
        ax = axes_by_panel[panel]
        letter, title, ylabel = PANEL_META[panel]
        curves = sort_curves(panel, {curve for current, curve in grouped if current == panel})
        data_max = args.self_plot_xmax if panel == "self_paralogue" else args.xmax
        axis_max = args.self_axis_xmax if panel == "self_paralogue" else args.xmax
        axis, ticks, bandwidth, reflect_at_zero = panel_axis(panel, axis_max)
        grid = ks_grid(args.points, axis_max)
        density_pairs: list[tuple[str, np.ndarray, int]] = []
        max_y = 0.0
        for curve in curves:
            raw_values = np.asarray(grouped[(panel, curve)], dtype=float)
            lower = args.self_min_ks if panel == "self_paralogue" else 0.0
            values = raw_values[(raw_values >= lower) & (raw_values <= data_max)]
            density, used_bandwidth = reflected_kde(
                values, grid, bandwidth, reflect_at_zero=reflect_at_zero,
            )
            maximum = int(np.argmax(density))
            density_pairs.append((curve, density, maximum))
            max_y = max(max_y, float(density[maximum]))
            color, linestyle, linewidth = curve_style(panel, curve)
            ax.plot(axis.forward(grid), density, color=color, linestyle=linestyle,
                    linewidth=linewidth, label=label_for(panel, curve))
            for true_x, display_x, y in zip(grid, axis.forward(grid), density):
                density_rows.append({
                    "panel": panel, "curve": curve, "Ks_grid": "{:.8f}".format(float(true_x)),
                    "display_x": "{:.8f}".format(float(display_x)), "density": "{:.8f}".format(float(y)),
                    "n": len(values), "bandwidth": "{:.8f}".format(used_bandwidth),
                })
            curve_rows.append({
                "panel": panel, "curve": curve, "n_pairs_0_lt_Ks_lt_5": len(values),
                "n_pairs_plotted_Ks_le_xmax": len(values),
                "n_pairs_removed_below_self_min": int(sum(raw_values < lower)),
                "n_pairs_excluded_saturation_tail": int(len(raw_values) - len(values)),
                "ks_min": "{:.8g}".format(float(values.min())),
                "ks_median": "{:.8g}".format(float(np.median(values))),
                "ks_max": "{:.8g}".format(float(values.max())),
                "bandwidth": "{:.8g}".format(used_bandwidth),
            })
        panel_densities[panel] = {
            "axis": axis, "ticks": ticks, "grid": grid,
            "curves": density_pairs, "max_y": max_y,
            "bandwidth": bandwidth, "reflect_at_zero": reflect_at_zero,
        }
        ax.set_ylim(0.0, max_y * (1.32 if panel == "self_paralogue" else 1.28))
        ax.set_ylabel(ylabel, fontsize=6.2, labelpad=2)
        ax.text(-0.13, 1.02, letter, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom")
        ax.text(0.0, 1.02, title, transform=ax.transAxes, fontsize=7.5, va="bottom")
        decorate_axis(ax, axis, ticks, xlabel=True)
        ax.tick_params(axis="y", direction="out", length=2.5, width=0.7, labelsize=5.5, pad=2)
        columns = 2 if len(curves) >= 6 else len(curves)
        legend_loc = "upper right" if panel == "self_paralogue" else "lower right"
        ax.legend(ncol=columns, loc=legend_loc, fontsize=4.6, handlelength=1.45,
                  columnspacing=0.55, labelspacing=0.18, borderaxespad=0.15)
        annotate_modes(ax, axis, grid, density_pairs, panel, peak_rows)

    # Keep the three left x axes independent, as in the reference panel.
    for panel in ("orthologue_pcam", "orthologue_ppe"):
        axes_by_panel[panel].tick_params(axis="x", labelbottom=False)

    figures = outdir / "figures"
    figures.mkdir()
    stem = figures / "figure2_repartition_five_ks"
    fig.savefig(str(stem) + ".pdf", bbox_inches="tight")
    fig.savefig(str(stem) + ".svg", bbox_inches="tight")
    fig.savefig(str(stem) + ".png", dpi=args.dpi, bbox_inches="tight")
    fig.savefig(str(stem) + ".tiff", dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    # Individual exports make each of the five Ks panels independently usable.
    for panel in PANEL_ORDER:
        letter, title, ylabel = PANEL_META[panel]
        info = panel_densities[panel]
        axis = info["axis"]
        ticks = info["ticks"]
        grid = info["grid"]
        curves = [name for name, _, _ in info["curves"]]
        max_y = float(info["max_y"])
        single, ax = plt.subplots(figsize=(3.55, 2.25), constrained_layout=True)
        for curve, density, _ in info["curves"]:
            color, linestyle, linewidth = curve_style(panel, curve)
            ax.plot(axis.forward(grid), density, color=color, linestyle=linestyle,
                    linewidth=linewidth, label=label_for(panel, curve))
        ax.set_ylim(0.0, max_y * (1.32 if panel == "self_paralogue" else 1.28))
        ax.set_ylabel(ylabel, fontsize=6.2, labelpad=2)
        ax.text(-0.13, 1.02, letter, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom")
        ax.text(0.0, 1.02, title, transform=ax.transAxes, fontsize=7.5, va="bottom")
        decorate_axis(ax, axis, ticks, xlabel=True)
        ax.tick_params(axis="y", direction="out", length=2.5, width=0.7, labelsize=5.5, pad=2)
        columns = 3 if len(curves) >= 6 else len(curves)
        legend_loc = "upper right" if panel == "self_paralogue" else "lower right"
        ax.legend(ncol=columns, loc=legend_loc, fontsize=4.6, handlelength=1.45,
                  columnspacing=0.55, labelspacing=0.18, borderaxespad=0.15)
        annotate_modes(ax, axis, grid, info["curves"], panel, [])
        single_stem = figures / "figure2_repartition_ks_{}".format(letter)
        single.savefig(str(single_stem) + ".pdf", bbox_inches="tight")
        single.savefig(str(single_stem) + ".svg", bbox_inches="tight")
        single.savefig(str(single_stem) + ".png", dpi=args.dpi, bbox_inches="tight")
        single.savefig(str(single_stem) + ".tiff", dpi=args.dpi, bbox_inches="tight")
        plt.close(single)

    write_tsv(source_dir / "curve_summary.tsv", curve_rows, list(curve_rows[0]))
    write_tsv(source_dir / "kde_peak_summary.tsv", peak_rows, list(peak_rows[0]))
    write_tsv(source_dir / "density_grid.tsv", density_rows, list(density_rows[0]))
    (outdir / "method_contract.txt").write_text(
        "All source pairs derive from corrected haploid_repartition packages or canonical PCAM/PPE/GTR inputs.\n"
        "External and A-B pairs: JCVI phase-specific anchors; self pairs: i-ADHoRe real cross-chromosome anchorpoints.\n"
        "Ks: {}.\n"
        "Retained values: 0<Ks<5. PZI/PSP phases are pooled only after pair-level Ks calculation.\n"
        "KDE is calculated on the original Ks scale; the display coordinate is a smooth monotonic nonlinear map.\n"
        "External and A-B panels use boundary reflection at Ks=0; self-synteny is shown without artificial low-Ks reflection.\n"
        "For the self-synteny panel, only {}<=Ks<=3 is used as data for the KDE; the KDE is evaluated to Ks={} so the retained-data tail decays smoothly, while values outside the displayed range are excluded.\n"
        "No median lines are drawn; numeric KDE modes are recorded in source_data/kde_peak_summary.tsv and labelled in the plots.\n".format(
            estimator, args.self_min_ks, args.self_axis_xmax),
        encoding="ascii",
    )
    (outdir / "plot_parameters.tsv").write_text(
        "parameter\tvalue\n"
        "xmax\t{}\nself_plot_xmax\t{}\nself_axis_xmax\t{}\nself_min_ks\t{}\npoints\t{}\nfixed_bandwidth\t{}\nminimum_curve_pairs\t{}\naxis\tcontinuous asinh plus smooth focus term; external/A-B range 0-5; self-synteny data range lower-cutoff-to-3 and KDE display range 0-{}\n".format(
            args.xmax, args.self_plot_xmax, args.self_axis_xmax, args.self_min_ks, args.points, args.bandwidth, args.min_curve_pairs, args.self_axis_xmax),
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("Corrected five-panel Ks figure completed.\n", encoding="ascii")
    print("[OK] Corrected five-panel Ks figure: {}".format(figures))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print("[ERROR] {}".format(error))
        raise SystemExit(2)
