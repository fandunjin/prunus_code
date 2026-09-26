#!/usr/bin/env python3
"""Create GTR-centred and PZI chromosome-centred syntenic retention heatmaps."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
import numpy as np


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 7,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    }
)


TRACK_ORDER = ["PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE"]
TRACK_LABELS = ["PZI-A1", "PZI-B1", "PSP-A1", "PSP-B1", "PPA-A", "PPA-B", "PCAM", "PPE"]
PZI_CHROMS = [f"C{i}A" for i in range(1, 9)] + [f"C{i}B" for i in range(1, 9)]

RETENTION_CMAP = LinearSegmentedColormap.from_list(
    "prunus_retention", ["#F3F1EC", "#D1DCE4", "#718BA2", "#7C2947"], N=256
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fieldnames:
        raise ValueError(f"Cannot infer columns for empty table: {path}")
    fields = fieldnames or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig: plt.Figure, prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{prefix}.svg", bbox_inches="tight")
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{prefix}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{prefix}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def require_columns(rows: list[dict[str, str]], required: set[str], path: Path) -> None:
    found = set(rows[0]) if rows else set()
    missing = required - found
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")


def build_gtr_matrix(
    gene_order_path: Path,
    presence_path: Path,
    window_size: int,
    step: int,
) -> tuple[np.ndarray, list[str], list[int], list[dict], list[dict]]:
    """Build a complete GTR chromosome-ordered matrix from binary track presence."""
    gene_rows = read_tsv(gene_order_path)
    require_columns(gene_rows, {"gtr_global_index", "gtr_seqid", "gtr_gene"}, gene_order_path)
    gene_rows = sorted(gene_rows, key=lambda r: int(r["gtr_global_index"]))
    by_chrom: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in gene_rows:
        by_chrom[row["gtr_seqid"]].append(row)
    chrom_order = sorted(by_chrom, key=lambda x: min(int(r["gtr_global_index"]) for r in by_chrom[x]))

    presence_rows = read_tsv(presence_path)
    require_columns(presence_rows, {"track_id", "gtr_global_index", "anchored"}, presence_path)
    presence: dict[str, set[int]] = defaultdict(set)
    for row in presence_rows:
        if row["track_id"] in TRACK_ORDER and row["anchored"] == "1":
            presence[row["track_id"]].add(int(row["gtr_global_index"]))

    matrix_rows: list[dict] = []
    matrix_values: list[list[float]] = []
    chromosome_centers: list[int] = []
    chromosome_boundaries: list[int] = []
    row_index = 0
    for chrom in chrom_order:
        genes = by_chrom[chrom]
        if len(genes) < window_size:
            raise ValueError(f"{chrom} has only {len(genes)} genes; cannot use {window_size}-gene windows")
        first_row = row_index
        starts = range(0, len(genes) - window_size + 1, step)
        for window_index, start in enumerate(starts, start=1):
            window = genes[start : start + window_size]
            indices = [int(r["gtr_global_index"]) for r in window]
            values = [sum(idx in presence[track] for idx in indices) / len(indices) for track in TRACK_ORDER]
            matrix_values.append(values)
            center = window[len(window) // 2]
            matrix_rows.append(
                {
                    "gtr_seqid": chrom,
                    "window_index": window_index,
                    "window_start_gene_order": window[0]["gtr_global_index"],
                    "window_end_gene_order": window[-1]["gtr_global_index"],
                    "window_center_gene": center["gtr_gene"],
                    "window_size": window_size,
                    "window_step": step,
                    **{f"{track}_retention": f"{value:.8f}" for track, value in zip(TRACK_ORDER, values)},
                }
            )
            row_index += 1
        chromosome_centers.append((first_row + row_index - 1) // 2)
        chromosome_boundaries.append(row_index - 0.5)

    matrix = np.asarray(matrix_values, dtype=float)
    summary_rows: list[dict] = []
    for chrom in chrom_order:
        chrom_rows = [r for r in matrix_rows if r["gtr_seqid"] == chrom]
        row_values = matrix[[r["gtr_seqid"] == chrom for r in matrix_rows]]
        summary_rows.append(
            {
                "gtr_seqid": chrom,
                "n_genes": len(by_chrom[chrom]),
                "n_windows": len(chrom_rows),
                "window_size": window_size,
                "window_step": step,
                **{f"{track}_mean_retention": f"{np.mean(row_values[:, j]):.8f}" for j, track in enumerate(TRACK_ORDER)},
            }
        )
    return matrix, chrom_order, chromosome_centers, chromosome_boundaries, matrix_rows, summary_rows


def build_pzi_matrix(
    anchors_path: Path,
    window_size: int,
) -> tuple[np.ndarray, list[dict], list[dict]]:
    """Build chromosome-centred rolling retention from quota-filtered anchors."""
    rows = read_tsv(anchors_path)
    require_columns(
        rows,
        {"track_id", "query_gene", "query_seqid", "query_start", "gtr_global_index", "quota_score"},
        anchors_path,
    )

    best: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        if row["track_id"] not in TRACK_ORDER:
            continue
        key = (row["track_id"], row["query_gene"])
        old = best.get(key)
        if old is None or int(row["quota_score"]) > int(old["quota_score"]):
            best[key] = row

    target_presence: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["track_id"] in TRACK_ORDER:
            target_presence[row["track_id"]].add(row["gtr_global_index"])

    source_by_chrom: dict[str, list[dict[str, str]]] = defaultdict(list)
    for (track, _query_gene), row in best.items():
        if track not in {"PZI_A", "PZI_B"}:
            continue
        prefix = f"{track}__"
        if row["query_seqid"].startswith(prefix):
            chrom = row["query_seqid"][len(prefix) :]
            if chrom in PZI_CHROMS:
                source_by_chrom[chrom].append(row)

    values: list[list[float]] = []
    summary_rows: list[dict] = []
    window_rows: list[dict] = []
    for chrom in PZI_CHROMS:
        source_track = "PZI_A" if chrom.endswith("A") else "PZI_B"
        ordered = sorted(source_by_chrom[chrom], key=lambda r: (int(r["query_start"]), r["query_gene"]))
        if len(ordered) < window_size:
            raise ValueError(f"{chrom} has only {len(ordered)} anchored genes; cannot use {window_size}-gene windows")
        gtr_indices = [r["gtr_global_index"] for r in ordered]
        row_values: list[float] = []
        for track in TRACK_ORDER:
            indicator = np.asarray([idx in target_presence[track] for idx in gtr_indices], dtype=float)
            rolling = np.convolve(indicator, np.ones(window_size), mode="valid") / window_size
            row_values.append(float(np.mean(rolling)))
            for index, density in enumerate(rolling, start=1):
                center = ordered[index - 1 + window_size // 2]
                window_rows.append(
                    {
                        "source_chromosome": chrom,
                        "source_track": source_track,
                        "track_id": track,
                        "window_index": index,
                        "window_size": window_size,
                        "window_center_query_gene": center["query_gene"],
                        "window_center_query_start": center["query_start"],
                        "retention_proportion": f"{density:.8f}",
                    }
                )
        values.append(row_values)
        summary_rows.append(
            {
                "source_chromosome": chrom,
                "source_track": source_track,
                "window_size": window_size,
                "n_source_anchored_genes": len(ordered),
                "n_complete_windows": len(ordered) - window_size + 1,
                **{f"{track}_mean_retention": f"{value:.8f}" for track, value in zip(TRACK_ORDER, row_values)},
            }
        )
    return np.asarray(values, dtype=float), summary_rows, window_rows


def style_heatmap_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)


def plot_gtr_heatmap(
    matrix: np.ndarray,
    chrom_order: list[str],
    centers: list[int],
    boundaries: list[int],
    output_prefix: Path,
) -> None:
    height = max(5.5, min(11.5, 4.2 + matrix.shape[0] * 0.018))
    fig, ax = plt.subplots(figsize=(5.35, height), facecolor="white")
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=RETENTION_CMAP, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(TRACK_LABELS)), TRACK_LABELS, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(centers, chrom_order)
    ax.tick_params(axis="x", labelsize=6, pad=3)
    ax.tick_params(axis="y", labelsize=6, pad=3)
    ax.set_xlabel("Subgenome track", labelpad=6)
    ax.set_ylabel("GTR reference chromosome", labelpad=6)
    for boundary in boundaries[:-1]:
        ax.axhline(boundary, color="white", linewidth=1.0)
    for x in [1.5, 3.5, 5.5]:
        ax.axvline(x, color="white", linewidth=1.0)
    ax.set_title("GTR-centred syntenic retention", loc="left", pad=8, fontsize=8, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.045, ticks=[0, 0.25, 0.5, 0.75, 1.0])
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x * 100:.0f}%"))
    cbar.ax.tick_params(length=2, width=0.5, labelsize=5.5)
    cbar.outline.set_linewidth(0.5)
    cbar.set_label("Syntenic gene retention", fontsize=6, labelpad=4)
    ax.text(-0.14, 1.02, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    style_heatmap_axis(ax)
    fig.subplots_adjust(left=0.18, right=0.89, top=0.93, bottom=0.18)
    save_figure(fig, output_prefix)


def plot_gtr_summary_heatmap(
    matrix: np.ndarray,
    chrom_order: list[str],
    output_prefix: Path,
) -> None:
    """Plot one mean-retention row per GTR chromosome for a compact main-panel option."""
    fig, ax = plt.subplots(figsize=(5.35, 3.05), facecolor="white")
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=RETENTION_CMAP, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(TRACK_LABELS)), TRACK_LABELS, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(len(chrom_order)), [label.replace("GTR__", "") for label in chrom_order])
    ax.tick_params(axis="x", labelsize=6, pad=3)
    ax.tick_params(axis="y", labelsize=6, pad=3)
    ax.set_xlabel("Subgenome track", labelpad=6)
    ax.set_ylabel("GTR reference chromosome", labelpad=6)
    for x in [1.5, 3.5, 5.5]:
        ax.axvline(x, color="white", linewidth=1.0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            ax.text(j, i, f"{value * 100:.0f}", ha="center", va="center", fontsize=5.2,
                    color="white" if value >= 0.42 else "#26313A")
    ax.set_title("GTR-centred chromosome summary", loc="left", pad=8, fontsize=8, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.045, ticks=[0, 0.25, 0.5, 0.75, 1.0])
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x * 100:.0f}%"))
    cbar.ax.tick_params(length=2, width=0.5, labelsize=5.5)
    cbar.outline.set_linewidth(0.5)
    cbar.set_label("Syntenic gene retention", fontsize=6, labelpad=4)
    ax.text(-0.14, 1.02, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    style_heatmap_axis(ax)
    fig.subplots_adjust(left=0.18, right=0.89, top=0.88, bottom=0.23)
    save_figure(fig, output_prefix)


def plot_pzi_heatmap(matrix: np.ndarray, output_prefix: Path, window_size: int) -> None:
    fig, ax = plt.subplots(figsize=(5.35, 4.7), facecolor="white")
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=RETENTION_CMAP, vmin=0.5, vmax=1.0)
    ax.set_xticks(np.arange(len(TRACK_LABELS)), TRACK_LABELS, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(len(PZI_CHROMS)), PZI_CHROMS)
    ax.tick_params(axis="x", labelsize=6, pad=3)
    ax.tick_params(axis="y", labelsize=6, pad=3)
    ax.set_xlabel("Subgenome track", labelpad=6)
    ax.set_ylabel("PZI haplotype chromosome", labelpad=6)
    for x in [1.5, 3.5, 5.5]:
        ax.axvline(x, color="white", linewidth=1.0)
    ax.axhline(7.5, color="white", linewidth=1.2)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            ax.text(j, i, f"{value * 100:.0f}", ha="center", va="center", fontsize=5.5,
                    color="white" if value >= 0.78 else "#26313A")
    ax.set_title(f"PZI chromosome-centred retention ({window_size}-gene windows)", loc="left", pad=8, fontsize=8, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.045, ticks=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x * 100:.0f}%"))
    cbar.ax.tick_params(length=2, width=0.5, labelsize=5.5)
    cbar.outline.set_linewidth(0.5)
    cbar.set_label("Syntenic gene retention", fontsize=6, labelpad=4)
    ax.text(-0.14, 1.02, "b", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    style_heatmap_axis(ax)
    fig.subplots_adjust(left=0.18, right=0.89, top=0.90, bottom=0.20)
    save_figure(fig, output_prefix)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--gtr-gene-order", type=Path, default=Path("0903/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact/01_catalog/gtr_gene_order.tsv"))
    parser.add_argument("--gtr-presence", type=Path, default=Path("0903/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact/04_quota_align/normalized/gtr_track_presence.tsv"))
    parser.add_argument("--anchors", type=Path, default=Path("0903/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact/04_quota_align/normalized/quota_filtered_anchors.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("0903/figures/figure2_gtr_pzi_heatmaps"))
    parser.add_argument("--gtr-window-size", type=int, default=100)
    parser.add_argument("--gtr-step", type=int, default=100)
    parser.add_argument("--pzi-window-size", type=int, default=601)
    return parser.parse_args()


def resolve(root: Path, path: Path) -> Path:
    return path if path.is_absolute() else root / path


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    output_dir = resolve(root, args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gtr_order = resolve(root, args.gtr_gene_order)
    gtr_presence = resolve(root, args.gtr_presence)
    anchors = resolve(root, args.anchors)
    for path in (gtr_order, gtr_presence, anchors):
        if not path.exists():
            raise FileNotFoundError(path)

    gtr_matrix, gtr_chroms, centers, boundaries, gtr_windows, gtr_summary = build_gtr_matrix(
        gtr_order, gtr_presence, args.gtr_window_size, args.gtr_step
    )
    # `centers` and `boundaries` are returned separately to keep the plotting code explicit.
    gtr_center_rows = centers
    gtr_boundary_rows = boundaries
    pzi_matrix, pzi_summary, pzi_windows = build_pzi_matrix(anchors, args.pzi_window_size)
    gtr_summary_matrix = np.asarray(
        [[float(row[f"{track}_mean_retention"]) for track in TRACK_ORDER] for row in gtr_summary],
        dtype=float,
    )

    # Source data are wide matrices plus long window tables for reproducibility.
    gtr_matrix_rows = []
    for index, row in enumerate(gtr_windows):
        gtr_matrix_rows.append({"matrix_row": index, **row})
    write_tsv(output_dir / "figure2a_gtr_windows.tsv", gtr_matrix_rows)
    write_tsv(output_dir / "figure2a_gtr_chromosome_summary.tsv", gtr_summary)
    write_tsv(output_dir / "figure2b_pzi_chromosome_summary.tsv", pzi_summary)
    write_tsv(output_dir / "figure2b_pzi_windows.tsv", pzi_windows)
    write_tsv(
        output_dir / "figure2a_gtr_matrix.tsv",
        [{"gtr_seqid": chrom, "matrix_row": i, **{track: f"{value:.8f}" for track, value in zip(TRACK_ORDER, row)}} for i, (chrom, row) in enumerate(zip([r["gtr_seqid"] for r in gtr_windows], gtr_matrix))],
    )
    write_tsv(
        output_dir / "figure2b_pzi_matrix.tsv",
        [{"source_chromosome": chrom, **{track: f"{value:.8f}" for track, value in zip(TRACK_ORDER, row)}} for chrom, row in zip(PZI_CHROMS, pzi_matrix)],
    )

    plot_gtr_heatmap(gtr_matrix, gtr_chroms, gtr_center_rows, gtr_boundary_rows, output_dir / "figure2a_gtr_centred")
    plot_gtr_summary_heatmap(gtr_summary_matrix, gtr_chroms, output_dir / "figure2a_gtr_centred_summary")
    plot_pzi_heatmap(pzi_matrix, output_dir / "figure2b_pzi_centred", args.pzi_window_size)

    manifest = [
        {"item": "figure2a", "value": "GTR-centred syntenic retention heatmap"},
        {"item": "figure2a_gene_order", "value": str(gtr_order)},
        {"item": "figure2a_presence", "value": str(gtr_presence)},
        {"item": "figure2a_window", "value": f"{args.gtr_window_size}-gene non-overlapping bins; step={args.gtr_step}"},
        {"item": "figure2a_rows", "value": f"{len(gtr_matrix)} GTR windows across {len(gtr_chroms)} reference chromosomes"},
        {"item": "figure2b", "value": "PZI chromosome-centred syntenic retention heatmap"},
        {"item": "figure2b_anchors", "value": str(anchors)},
        {"item": "figure2b_window", "value": f"{args.pzi_window_size}-gene rolling windows; complete windows only"},
        {"item": "figure2b_rows", "value": "C1A-C8A and C1B-C8B from PZI-A1 and PZI-B1"},
        {"item": "tracks", "value": ",".join(TRACK_ORDER)},
        {"item": "exports", "value": "SVG, PDF, 600-dpi PNG and TIFF"},
    ]
    write_tsv(output_dir / "figure2_gtr_pzi_heatmaps_manifest.tsv", manifest, ["item", "value"])
    (output_dir / "figure2_gtr_pzi_heatmaps_methods.txt").write_text(
        "Version A (GTR-centred): GTR genes were ordered by gtr_global_index and grouped by GTR chromosome. "
        f"For each complete non-overlapping {args.gtr_window_size}-gene bin (step {args.gtr_step}), the value for a track is the proportion of GTR genes with anchored=1 in gtr_track_presence.tsv. "
        "Version B (PZI chromosome-centred): quota_filtered_anchors.tsv was deduplicated to the highest quota_score per query gene, then PZI-A1 and PZI-B1 genes were ordered by query coordinate. "
        f"For each PZI chromosome, complete rolling windows of {args.pzi_window_size} source genes were calculated; each target value is the proportion of those source genes whose GTR locus is present in that target track. "
        "Both panels measure quota-filtered syntenic retention, not validated gene loss. GTR is a low-copy reference gene-order catalog and is not treated as a validated ancestral karyotype.\n",
        encoding="utf-8",
    )
    (output_dir / "figure2_gtr_pzi_heatmaps_legend.txt").write_text(
        "Figure 2 | GTR-centred and PZI chromosome-centred syntenic gene retention. "
        "a, GTR reference genes are grouped into chromosome-ordered 100-gene bins; colors show the proportion of loci with a quota-filtered anchor in each Prunus track. "
        "b, PZI-A1 and PZI-B1 chromosomes are represented by complete 601-gene rolling windows; colors and cell labels show the corresponding target-track retention density. "
        "The eight tracks are PZI-A1, PZI-B1, PSP-A1, PSP-B1, PPA-A, PPA-B, PPE and PCAM. "
        "Values represent observable syntenic anchor retention and should not be interpreted as confirmed gene-loss calls or as a complete Prunus ancestral chromosome map. Source data are provided as Source Data files.\n",
        encoding="utf-8",
    )
    print(f"Wrote GTR-centred and PZI-centred heatmaps to {output_dir}")


if __name__ == "__main__":
    main()
