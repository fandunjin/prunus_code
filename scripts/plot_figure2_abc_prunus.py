#!/usr/bin/env python3
"""Create Nature-style Figure 2a-c for Prunus subgenome retention and TE composition."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter, MultipleLocator
import numpy as np


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    }
)


TRACK_ORDER = ["PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE"]
TRACK_LABELS = ["PZI-A1", "PZI-B1", "PSP-A1", "PSP-B1", "PPA-A", "PPA-B", "PCAM", "PPE"]
CHROM_ORDER = [f"C{i}A" for i in range(1, 9)] + [f"C{i}B" for i in range(1, 9)]

GENE_COUNT_ORDER = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "ppe", "pcam",
]
GENE_COUNT_LABELS = [
    "PZI-A1", "PZI-A2", "PZI-B1", "PZI-B2",
    "PSP-A1", "PSP-A2", "PSP-B1", "PSP-B2",
    "PPA-A", "PPA-B", "PPE", "PCAM",
]

TE_ORDER = [
    "Copia", "Gypsy", "CACTA", "hAT", "Mutator", "PIF", "Tc1", "Helitron", "Other TEs",
]
TE_COLORS = {
    "Copia": "#334E7D",
    "Gypsy": "#2E8A78",
    "CACTA": "#C94C43",
    "hAT": "#D28B2D",
    "Mutator": "#8065A1",
    "PIF": "#B65E8A",
    "Tc1": "#B8A14A",
    "Helitron": "#5F6268",
    "Other TEs": "#A3A3A3",
}

DISPLAY_ORDER = ["PZI-A", "PZI-B", "PSP-A", "PSP-B", "PPA-A", "PPA-B", "PPE", "PCAM"]
DISPLAY_TARGETS = {
    "PZI-A": ["pziA1", "pziA2"],
    "PZI-B": ["pziB1", "pziB2"],
    "PSP-A": ["pspA1", "pspA2"],
    "PSP-B": ["pspB1", "pspB2"],
    "PPA-A": ["ppaA"],
    "PPA-B": ["ppaB"],
    "PPE": ["ppe"],
    "PCAM": ["pcam"],
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def count_gff_genes(path: Path) -> int:
    """Count unique gene features, using the feature ID where available."""
    ids: set[str] = set()
    n_without_id = 0
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2].lower() != "gene":
                continue
            attrs = fields[8].split(";")
            gene_id = ""
            for attr in attrs:
                key, sep, value = attr.partition("=")
                if sep and key in {"ID", "gene_id", "locus_tag"}:
                    gene_id = value
                    break
                key, sep, value = attr.strip().partition(" ")
                if sep and key in {"gene_id", "locus_tag"}:
                    gene_id = value.strip('"')
                    break
            if gene_id:
                ids.add(gene_id)
            else:
                n_without_id += 1
    return len(ids) + n_without_id


def compute_chromosome_fractionation(
    anchors_path: Path, window_size: int
) -> tuple[list[dict], list[dict], np.ndarray]:
    """Calculate 601-gene retention windows along all eight PZI A/B chromosomes."""
    rows = read_tsv(anchors_path)
    required = {"track_id", "query_gene", "query_seqid", "query_start", "gtr_global_index", "quota_score"}
    missing = required.difference(rows[0]) if rows else required
    if missing:
        raise ValueError(f"Missing columns in {anchors_path}: {sorted(missing)}")

    # Keep the highest-scoring ancestral anchor for each query gene.
    best: dict[tuple[str, str], dict] = {}
    for row in rows:
        if row["track_id"] not in TRACK_ORDER:
            continue
        key = (row["track_id"], row["query_gene"])
        score = int(row["quota_score"])
        old = best.get(key)
        if old is None or score > int(old["quota_score"]):
            best[key] = row

    target_presence: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["track_id"] in TRACK_ORDER:
            target_presence[row["track_id"]].add(row["gtr_global_index"])

    source_by_chrom: dict[str, list[dict]] = defaultdict(list)
    for (track, _query_gene), row in best.items():
        if track not in {"PZI_A", "PZI_B"}:
            continue
        prefix = "PZI_A__" if track == "PZI_A" else "PZI_B__"
        if not row["query_seqid"].startswith(prefix):
            continue
        chrom_label = row["query_seqid"].removeprefix(prefix)
        if chrom_label in CHROM_ORDER:
            source_by_chrom[chrom_label].append(row)

    summary_rows: list[dict] = []
    window_rows: list[dict] = []
    matrix = np.full((len(CHROM_ORDER), len(TRACK_ORDER)), np.nan, dtype=float)
    for i, chrom_label in enumerate(CHROM_ORDER):
        source_track = "PZI_A" if chrom_label.endswith("A") else "PZI_B"
        ordered = sorted(source_by_chrom[chrom_label], key=lambda x: (int(x["query_start"]), x["query_gene"]))
        if len(ordered) < window_size:
            raise ValueError(f"{chrom_label} has only {len(ordered)} anchored genes; cannot use {window_size}-gene windows")
        gtr_indices = [row["gtr_global_index"] for row in ordered]
        n_windows = len(ordered) - window_size + 1
        for j, track in enumerate(TRACK_ORDER):
            anchored = np.array([idx in target_presence[track] for idx in gtr_indices], dtype=float)
            rolling = np.convolve(anchored, np.ones(window_size, dtype=float), mode="valid") / window_size
            matrix[i, j] = float(np.mean(rolling))
            for k, density in enumerate(rolling, start=1):
                center = ordered[k - 1 + window_size // 2]
                window_rows.append(
                    {
                        "source_chromosome": chrom_label,
                        "source_track": source_track,
                        "track_id": track,
                        "source_query_seqid": center["query_seqid"],
                        "window_index": k,
                        "window_size": window_size,
                        "window_center_query_gene": center["query_gene"],
                        "window_center_query_start": center["query_start"],
                        "retention_proportion": f"{density:.8f}",
                    }
                )
            summary_rows.append(
                {
                    "source_chromosome": chrom_label,
                    "source_track": source_track,
                    "track_id": track,
                    "source_query_seqid": ordered[0]["query_seqid"],
                    "window_size": window_size,
                    "n_source_anchored_genes": len(ordered),
                    "n_complete_windows": n_windows,
                    "mean_retention_proportion": f"{np.mean(rolling):.8f}",
                    "median_retention_proportion": f"{np.median(rolling):.8f}",
                    "min_retention_proportion": f"{np.min(rolling):.8f}",
                    "max_retention_proportion": f"{np.max(rolling):.8f}",
                }
            )
    return summary_rows, window_rows, matrix


def load_gene_counts(root: Path) -> tuple[list[dict], dict[str, int]]:
    paths = {
        "pziA1": root / "0903/haploid_repartition/inputs/pzi/gff/pzi_hapA1.gff3",
        "pziA2": root / "0903/haploid_repartition/inputs/pzi/gff/pzi_hapA2.gff3",
        "pziB1": root / "0903/haploid_repartition/inputs/pzi/gff/pzi_hapB1.gff3",
        "pziB2": root / "0903/haploid_repartition/inputs/pzi/gff/pzi_hapB2.gff3",
        "pspA1": root / "0903/haploid_repartition/inputs/psp/gff/psp_hapA1.gff3",
        "pspA2": root / "0903/haploid_repartition/inputs/psp/gff/psp_hapA2.gff3",
        "pspB1": root / "0903/haploid_repartition/inputs/psp/gff/psp_hapB1.gff3",
        "pspB2": root / "0903/haploid_repartition/inputs/psp/gff/psp_hapB2.gff3",
        "ppaA": root / "0903/haploid_repartition/inputs/ppa/gff/ppa_hapA.gff3",
        "ppaB": root / "0903/haploid_repartition/inputs/ppa/gff/ppa_hapB.gff3",
        "ppe": root / "data/ppe/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.gff.gz",
        "pcam": root / "data/pcam/GCA_050084535.1_PCAv2.0_genomic.gff.gz",
    }
    missing = [f"{name}: {path}" for name, path in paths.items() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing GFF inputs:\n" + "\n".join(missing))
    counts = {target: count_gff_genes(path) for target, path in paths.items()}
    rows = [
        {"target": target, "display_label": label, "gene_count": counts[target], "source": str(paths[target])}
        for target, label in zip(GENE_COUNT_ORDER, GENE_COUNT_LABELS)
    ]
    return rows, counts


def load_te_summary(path: Path) -> tuple[list[dict], list[dict]]:
    rows = read_tsv(path)
    required = {"target", "category", "bp", "genome_bp", "genome_fraction"}
    missing = required.difference(rows[0]) if rows else required
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    raw_targets = {row["target"] for row in rows}
    expected_targets = {target for targets in DISPLAY_TARGETS.values() for target in targets}
    missing_targets = expected_targets.difference(raw_targets)
    if missing_targets:
        raise ValueError(f"Missing TE targets in {path}: {sorted(missing_targets)}")

    broad_map = {
        "LTR/Copia": "Copia",
        "LTR/Gypsy": "Gypsy",
        "Helitron": "Helitron",
        "TIR/DNA": "Other TEs",
        "LINE": "Other TEs",
        "SINE": "Other TEs",
        "Other LTR": "Other TEs",
        "Other TE": "Other TEs",
    }
    by_target_category: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        mapped = broad_map.get(row["category"])
        if mapped in TE_ORDER:
            by_target_category[(row["target"], mapped)].append(row)

    summary: list[dict] = []
    for display in DISPLAY_ORDER:
        targets = DISPLAY_TARGETS[display]
        target_rows = [row for row in rows if row["target"] in targets]
        genome_bp = float(np.mean([float(row["genome_bp"]) for row in target_rows]))
        for category in TE_ORDER:
            category_rows = [row for target in targets for row in by_target_category[(target, category)]]
            bp = float(np.mean([float(row["bp"]) for row in category_rows])) if category_rows else 0.0
            summary.append(
                {
                    "display_label": display,
                    "category": category,
                    "bp": f"{bp:.2f}",
                    "mb": f"{bp / 1e6:.8f}",
                    "genome_bp": f"{genome_bp:.2f}",
                    "genome_fraction": f"{bp / genome_bp:.8f}" if genome_bp else "0",
                }
            )
    return summary, rows


def load_te_family_summary(
    composition_path: Path, family_path: Path
) -> tuple[list[dict], list[dict]]:
    """Aggregate exclusive pan-family bases into the reference paper's TE classes."""
    composition_rows = read_tsv(composition_path)
    family_rows = read_tsv(family_path)
    required = {"target", "class_family", "exclusive_bp"}
    missing = required.difference(family_rows[0]) if family_rows else required
    if missing:
        raise ValueError(f"Missing columns in {family_path}: {sorted(missing)}")

    def map_class(class_family: str) -> str:
        if class_family == "LTR/Copia":
            return "Copia"
        if class_family == "LTR/Gypsy":
            return "Gypsy"
        if class_family in {"DNA/CMC-EnSpm", "DNA/DTC", "MITE/DTC"}:
            return "CACTA"
        if class_family in {"DNA/DTA", "MITE/DTA", "DNA/hAT-Ac", "DNA/hAT-Tag1", "DNA/hAT-Tip100"}:
            return "hAT"
        if class_family in {"DNA/DTM", "MITE/DTM", "DNA/MULE-MuDR"}:
            return "Mutator"
        if class_family in {"DNA/DTH", "MITE/DTH", "DNA/PIF-Harbinger"}:
            return "PIF"
        if class_family in {"DNA/DTT", "MITE/DTT"}:
            return "Tc1"
        if class_family in {"DNA/Helitron", "RC/Helitron"}:
            return "Helitron"
        return "Other TEs"

    expected_targets = {target for targets in DISPLAY_TARGETS.values() for target in targets}
    raw_targets = {row["target"] for row in family_rows}
    missing_targets = expected_targets.difference(raw_targets)
    if missing_targets:
        raise ValueError(f"Missing TE family targets in {family_path}: {sorted(missing_targets)}")
    genome_bp_by_target = {}
    for row in composition_rows:
        genome_bp_by_target.setdefault(row["target"], float(row["genome_bp"]))

    totals: dict[tuple[str, str], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in family_rows:
        display_targets = [display for display, targets in DISPLAY_TARGETS.items() if row["target"] in targets]
        if not display_targets:
            continue
        display = display_targets[0]
        totals[(display, map_class(row["class_family"]))][row["target"]] += float(row["exclusive_bp"])

    summary: list[dict] = []
    for display in DISPLAY_ORDER:
        targets = DISPLAY_TARGETS[display]
        genome_bp = float(np.mean([genome_bp_by_target[target] for target in targets]))
        for category in TE_ORDER:
            values = [totals[(display, category)].get(target, 0.0) for target in targets]
            bp = float(np.mean(values))
            summary.append(
                {
                    "display_label": display,
                    "category": category,
                    "bp": f"{bp:.2f}",
                    "mb": f"{bp / 1e6:.8f}",
                    "genome_bp": f"{genome_bp:.2f}",
                    "genome_fraction": f"{bp / genome_bp:.8f}" if genome_bp else "0",
                }
            )
    return summary, family_rows


def thousands(value: float, _position: int) -> str:
    return f"{value / 1000:.0f}"


def mb(value: float, _position: int) -> str:
    return f"{value / 1e6:.0f}"


def pct(value: float, _position: int) -> str:
    return f"{value * 100:.0f}%"


def plot_figure(
    matrix: np.ndarray,
    fractionation_summary: list[dict],
    gene_rows: list[dict],
    te_summary: list[dict],
    output_prefix: Path,
) -> None:
    fig = plt.figure(figsize=(7.2, 4.55), facecolor="white")
    outer = fig.add_gridspec(
        1, 3, width_ratios=[2.25, 2.05, 3.85], wspace=0.54,
        left=0.07, right=0.985, top=0.90, bottom=0.22,
    )
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    cgrid = outer[0, 2].subgridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.16)
    ax_c_len = fig.add_subplot(cgrid[0, 0])
    ax_c_pct = fig.add_subplot(cgrid[0, 1], sharey=ax_c_len)

    # Panel a: 601-gene rolling retention density.
    retention_cmap = LinearSegmentedColormap.from_list(
        "retention", ["#F4F1EC", "#CBD8E6", "#6283A2", "#8B334B"], N=256
    )
    im = ax_a.imshow(matrix, aspect="auto", interpolation="nearest", cmap=retention_cmap, vmin=0, vmax=0.70)
    ax_a.set_xticks(range(len(TRACK_LABELS)), TRACK_LABELS, rotation=48, ha="right", rotation_mode="anchor")
    ax_a.set_yticks(range(len(CHROM_ORDER)), CHROM_ORDER)
    ax_a.tick_params(axis="x", pad=3, length=0, labelsize=5.5)
    ax_a.tick_params(axis="y", length=0, labelsize=6)
    ax_a.set_xlabel("Subgenome track", labelpad=7)
    ax_a.set_ylabel("PZI haplotype chromosome", labelpad=5)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax_a.text(j, i, f"{matrix[i, j] * 100:.0f}", ha="center", va="center", fontsize=5.5,
                      color="white" if matrix[i, j] > 0.42 else "#26313A")
    for x in [1.5, 3.5, 5.5]:
        ax_a.axvline(x, color="white", linewidth=1.2)
    ax_a.axhline(7.5, color="white", linewidth=1.4)
    ax_a.set_title("Syntenic gene density", loc="left", pad=8, fontsize=8, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax_a, fraction=0.045, pad=0.04, ticks=[0, 0.2, 0.4, 0.6])
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(pct))
    cbar.ax.tick_params(length=2, width=0.5, labelsize=5.5)
    cbar.outline.set_linewidth(0.5)
    cbar.set_label("Syntenic genes", fontsize=6, labelpad=4)

    # Panel b: gene counts from annotation GFF3 files.
    gene_values = np.array([next(row["gene_count"] for row in gene_rows if row["target"] == target) for target in GENE_COUNT_ORDER], dtype=float)
    y = np.arange(len(gene_values))
    colors = ["#7A3E48", "#B77A82", "#1C5965", "#6D9AA3"] * 2 + ["#7A3E48", "#1C5965", "#4E5962", "#4E5962"]
    ax_b.barh(y, gene_values, height=0.68, color=colors, edgecolor="white", linewidth=0.35)
    ax_b.set_yticks(y, GENE_COUNT_LABELS)
    ax_b.invert_yaxis()
    ax_b.tick_params(axis="y", length=0, labelsize=6, pad=3)
    ax_b.tick_params(axis="x", labelsize=6)
    ax_b.xaxis.set_major_formatter(FuncFormatter(thousands))
    ax_b.xaxis.set_major_locator(MultipleLocator(5000))
    ax_b.set_xlim(0, max(gene_values) * 1.16)
    ax_b.set_xlabel("Genes (thousands)", labelpad=5)
    ax_b.set_title("Annotated gene number", loc="left", pad=8, fontsize=8, fontweight="bold")
    for boundary in [3.5, 7.5, 9.5, 10.5]:
        ax_b.axhline(boundary, color="#D7D7D7", linewidth=0.55, zorder=0)
    for yi, val in zip(y, gene_values):
        ax_b.text(val + max(gene_values) * 0.018, yi, f"{int(val):,}", va="center", fontsize=5.2, color="#36414A")

    # Panel c: TE length and fraction with a shared y-axis.
    te_by_display: dict[str, dict[str, float]] = defaultdict(dict)
    for row in te_summary:
        te_by_display[row["display_label"]][row["category"]] = float(row["bp"])
    cy = np.arange(len(DISPLAY_ORDER))
    ax_c_len.set_yticks(cy, DISPLAY_ORDER)
    ax_c_len.invert_yaxis()
    ax_c_len.tick_params(axis="y", length=0, labelsize=6, pad=3)
    ax_c_pct.tick_params(axis="y", left=False, labelleft=False)
    for ax in (ax_c_len, ax_c_pct):
        ax.tick_params(axis="x", labelsize=6)
        ax.grid(axis="x", color="#E4E5E6", linewidth=0.45, zorder=0)
        ax.set_axisbelow(True)
    left_base = np.zeros(len(DISPLAY_ORDER))
    right_base = np.zeros(len(DISPLAY_ORDER))
    for category in TE_ORDER:
        vals_bp = np.array([te_by_display[label][category] for label in DISPLAY_ORDER], dtype=float)
        vals_fraction = vals_bp / np.array([
            float(next(row["genome_bp"] for row in te_summary if row["display_label"] == label))
            for label in DISPLAY_ORDER
        ])
        ax_c_len.barh(cy, vals_bp, left=left_base, height=0.68, color=TE_COLORS[category], edgecolor="white", linewidth=0.25, label=category)
        ax_c_pct.barh(cy, vals_fraction, left=right_base, height=0.68, color=TE_COLORS[category], edgecolor="white", linewidth=0.25)
        left_base += vals_bp
        right_base += vals_fraction
    ax_c_len.xaxis.set_major_formatter(FuncFormatter(mb))
    ax_c_len.set_xlim(0, float(left_base.max()) * 1.10)
    ax_c_len.set_xlabel("TE length (Mb)", labelpad=5)
    ax_c_pct.xaxis.set_major_formatter(FuncFormatter(pct))
    ax_c_pct.set_xlim(0, 0.82)
    ax_c_pct.set_xlabel("TE / genome", labelpad=5)
    ax_c_len.set_title("TE composition", loc="left", pad=8, fontsize=8, fontweight="bold")
    ax_c_pct.spines["left"].set_visible(False)
    for boundary in [1.5, 3.5, 5.5]:
        ax_c_len.axhline(boundary, color="#D7D7D7", linewidth=0.55, zorder=0)
        ax_c_pct.axhline(boundary, color="#D7D7D7", linewidth=0.55, zorder=0)

    # Panel labels and one compact shared legend.
    for ax, label in [(ax_a, "a"), (ax_b, "b"), (ax_c_len, "c")]:
        ax.text(-0.16, 1.07, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    handles, labels = ax_c_len.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.64, 0.035), ncol=4,
               frameon=False, fontsize=5.8, handlelength=1.0, columnspacing=1.0, handletextpad=0.35)

    # Keep c axes aligned and avoid duplicate labels from the shared y-axis.
    ax_c_pct.tick_params(axis="y", labelsize=0)
    for ax in (ax_a, ax_b, ax_c_len, ax_c_pct):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.savefig(f"{output_prefix}.svg", bbox_inches="tight")
    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{output_prefix}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(f"{output_prefix}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--anchors",
        type=Path,
        default=Path("0903/downstream_analysis/figure2_polyploid/strict_wgd_mapping/04_fractbias_8track_reference_exact/04_quota_align/normalized/quota_filtered_anchors.tsv"),
    )
    parser.add_argument("--te", type=Path, default=Path("analysis/prunus_pan_te/postprocess_v2/te_composition.tsv"))
    parser.add_argument("--te-family", type=Path, default=Path("analysis/prunus_pan_te/postprocess_v2/te_exclusive_family_summary.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("0903/figures/figure2_abc"))
    parser.add_argument("--window-size", type=int, default=601)
    return parser.parse_args()


def resolve(root: Path, path: Path) -> Path:
    return path if path.is_absolute() else root / path


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    output_dir = resolve(root, args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fractionation_path = resolve(root, args.anchors)
    te_path = resolve(root, args.te)
    te_family_path = resolve(root, args.te_family)

    fractionation_summary, window_rows, matrix = compute_chromosome_fractionation(fractionation_path, args.window_size)
    gene_rows, _ = load_gene_counts(root)
    if te_family_path.exists():
        te_summary, te_raw = load_te_family_summary(te_path, te_family_path)
        te_method = "exclusive pan-family bases aggregated to Copia/Gypsy/CACTA/hAT/Mutator/PIF/Tc1/Helitron/Other TEs"
    else:
        te_summary, te_raw = load_te_summary(te_path)
        te_method = "normalized pan-TE classes from te_composition.tsv (family summary not found)"

    write_tsv(output_dir / "figure2a_601window_summary.tsv", list(fractionation_summary[0]), fractionation_summary)
    write_tsv(output_dir / "figure2a_601window_values.tsv", list(window_rows[0]), window_rows)
    write_tsv(output_dir / "figure2b_gene_counts.tsv", list(gene_rows[0]), gene_rows)
    write_tsv(output_dir / "figure2c_te_summary.tsv", list(te_summary[0]), te_summary)

    manifest_rows = [
        {"item": "figure", "value": "figure2_abc_prunus"},
        {"item": "panel_a_anchor_input", "value": str(fractionation_path)},
        {"item": "panel_a_method", "value": f"complete {args.window_size}-gene rolling windows; target presence / {args.window_size} source PZI-anchored genes"},
        {"item": "panel_a_scope", "value": "all eight PZI-A and eight PZI-B haplotype chromosomes; rows labelled C1A-C8A and C1B-C8B"},
        {"item": "panel_b_input", "value": "gene features counted from PZI/PSP/PPA/PCAM/PPE GFF3 annotations"},
        {"item": "panel_c_input", "value": str(te_path)},
        {"item": "panel_c_family_input", "value": str(te_family_path) if te_family_path.exists() else "not available"},
        {"item": "panel_c_method", "value": te_method + "; PZI/PSP A1+A2 and B1+B2 averaged within subgenome"},
        {"item": "window_size", "value": args.window_size},
        {"item": "exports", "value": "SVG, PDF, 600-dpi TIFF and PNG"},
    ]
    write_tsv(output_dir / "figure2_abc_manifest.tsv", ["item", "value"], manifest_rows)
    (output_dir / "figure2_abc_methods.txt").write_text(
        "Panel a follows the reference study's 601-gene-window procedure. For each PZI haplotype chromosome (C1A-C8A and C1B-C8B), "
        "source genes with a quota-filtered GTR anchor were ordered by genomic coordinate. For each complete window containing 300 upstream, "
        "the central gene and 300 downstream source genes, the density is the proportion with a corresponding anchor in each target track. "
        "The rows therefore cover all 16 PZI haplotype chromosomes and are not a replacement for a complete ancestral 22-block map. Panel b counts gene features in the source GFF3 files. "
        "Panel c uses the normalized Prunus pan-TE composition table; TE base pairs are averaged across A1/A2 or B1/B2 for PZI and PSP, "
        "and the right panel divides each class by the corresponding mean haplotype genome size. TE families were aggregated to the "
        "reference categories Copia, Gypsy, CACTA, hAT, Mutator, PIF, Tc1, Helitron and Other TEs using exclusive pan-family base pairs.\n",
        encoding="utf-8",
    )
    (output_dir / "figure2_abc_legend.txt").write_text(
        "Figure 2 | Syntenic gene retention, gene number and transposable-element composition across Prunus genomes. "
        "a, Syntenic gene retention density calculated in complete 601-gene windows (300 genes on either side of the central gene) "
        "for all eight PZI-A and eight PZI-B haplotype chromosomes. Values are mean retained-gene proportions across all complete windows; "
        "numbers in cells are percentages. b, Number of gene features in the corresponding GFF3 annotations. c, Left, stacked TE length; "
        "right, the same classes expressed as a fraction of genome size. PZI-A/B and PSP-A/B summarize the mean of A1/A2 or B1/B2, respectively. "
        "TE classes were assigned using the normalized Prunus pan-TE library and aggregated to the reference TE classes.\n",
        encoding="utf-8",
    )
    plot_figure(matrix, fractionation_summary, gene_rows, te_summary, output_dir / "figure2_abc_prunus")
    print(f"Wrote Figure 2a-c outputs to {output_dir}")


if __name__ == "__main__":
    main()
