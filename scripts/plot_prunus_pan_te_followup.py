#!/usr/bin/env python3
"""Publication-oriented plots for the TE/TRF follow-up analyses."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib import transforms


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.75,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

GROUPS = ["PZI", "PSP", "PPA", "PCAM", "PPE"]
PZI_TARGETS = ["pziA1", "pziA2", "pziB1", "pziB2"]
GROUP_COLORS = {"PZI": "#C44E52", "PSP": "#2A9D8F", "PPA": "#4C78A8", "PCAM": "#777777", "PPE": "#B279A2"}
TE_COLORS = {"LTR": "#4C78A8", "TIR/DNA": "#E39C45", "Other": "#8E3B5B"}


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    records = []
    name = None
    chunks = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records


def save(fig, out_dir, stem):
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def panel_label(ax, value):
    ax.text(-0.12, 1.04, value, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def plot_mechanism(stats_dir, out_dir):
    drivers = read_tsv(stats_dir / "te_family_driver_top.tsv")[:15]
    age_rows = read_tsv(stats_dir / "te_divergence_age.tsv")
    solo_rows = read_tsv(stats_dir / "solo_ltr_by_target_lineage.tsv")

    fig = plt.figure(figsize=(7.25, 5.2))
    grid = fig.add_gridspec(2, 2, width_ratios=(1.2, 1.0), hspace=0.58, wspace=0.48)

    ax = fig.add_subplot(grid[0, 0])
    drivers = list(reversed(drivers))
    y = np.arange(len(drivers))
    vals = np.array([float(row["pzi_psp_excess_bp"]) / 1e6 for row in drivers])
    colors = [TE_COLORS["LTR"] if row["class_family"].lower().startswith("ltr")
              else TE_COLORS["TIR/DNA"] if row["class_family"].lower().startswith(("dna", "tir"))
              else TE_COLORS["Other"] for row in drivers]
    ax.barh(y, vals, color=colors, height=0.68, linewidth=0)
    ax.set_yticks(y, [row["pan_family"] for row in drivers], fontsize=5.5)
    ax.set_xlabel("PZI − PSP excess (Mb)"); ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    ax.set_xlim(0, max(vals) * 1.18 if len(vals) else 1)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=TE_COLORS["LTR"]),
                       plt.Rectangle((0, 0), 1, 1, color=TE_COLORS["TIR/DNA"]),
                       plt.Rectangle((0, 0), 1, 1, color=TE_COLORS["Other"])],
              labels=["LTR", "TIR/DNA", "Other TE"], fontsize=5.2, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.02))
    panel_label(ax, "a")

    ax = fig.add_subplot(grid[0, 1])
    # Divergence is weighted by exclusive hit span, so the curve reflects
    # genomic sequence rather than the number of RepeatMasker records.
    bins = np.arange(0, 105, 5.0)
    group_values = {group: np.zeros(len(bins) - 1) for group in GROUPS}
    for row in age_rows:
        age = float(row["estimated_age_mya"])
        if age < 0 or age >= bins[-1]:
            continue
        group = row.get("group", "")
        if group not in group_values:
            continue
        idx = np.searchsorted(bins, age, side="right") - 1
        group_values[group][idx] += float(row.get("exclusive_bp", 0)) / 1e6
    for group in GROUPS:
        ax.step((bins[:-1] + bins[1:]) / 2, group_values[group], where="mid", color=GROUP_COLORS[group], lw=1.4, label=group)
    ax.set_xlabel("Approximate TE age (Mya)"); ax.set_ylabel("Exclusive TE span (Mb per 5-Mya bin)")
    ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True); ax.legend(fontsize=5.3, ncol=2, loc="upper right")
    ax.text(0.01, -0.25, "Age uses posterior mean μ; consensus divergence is approximate.", transform=ax.transAxes, fontsize=5.0, color="#555555")
    panel_label(ax, "b")

    ax = fig.add_subplot(grid[1, 0])
    totals = defaultdict(lambda: [0, 0])
    for row in solo_rows:
        group = row.get("group", "")
        if group in GROUPS:
            totals[group][0] += int(float(row.get("intact_count", 0)))
            totals[group][1] += int(float(row.get("solo_count", 0)))
    ratios = [totals[group][1] / totals[group][0] if totals[group][0] else 0 for group in GROUPS]
    ax.bar(GROUPS, ratios, color=[GROUP_COLORS[g] for g in GROUPS], width=0.68)
    ax.set_ylabel("solo LTR / intact LTR"); ax.set_xlabel("Genomic package group")
    ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    ax.text(0.01, -0.25, "Ratios pool annotated LTR lineages within each group.", transform=ax.transAxes, fontsize=5.0, color="#555555")
    panel_label(ax, "c")

    ax = fig.add_subplot(grid[1, 1])
    # A compact chromosome-location summary: mean TE and TRF fractions by
    # chromosome for the four PZI haplotypes. Window-level heatmaps are saved
    # separately to preserve resolution and avoid unreadable labels here.
    density_path = stats_dir / "te_trf_1mb_density.tsv"
    if density_path.exists():
        density = read_tsv(density_path)
        by_chr = defaultdict(lambda: defaultdict(list))
        for row in density:
            if row["target"] not in PZI_TARGETS:
                continue
            chrom = row["seqid"].split("_")[0]
            by_chr[chrom]["te"].append(float(row["te_fraction"]))
            by_chr[chrom]["trf"].append(float(row["trf_fraction"]))
        chroms = sorted(by_chr, key=lambda x: int(x.replace("chr", "")))
        x = np.arange(len(chroms)); width = 0.34
        te = [np.mean(by_chr[c]["te"]) * 100 for c in chroms]
        trf = [np.mean(by_chr[c]["trf"]) * 100 for c in chroms]
        ax.bar(x - width / 2, te, width, color="#5476B8", label="pan-TE")
        ax.bar(x + width / 2, trf, width, color="#C44E52", label="TRF")
        ax.set_xticks(x, [c.replace("chr", "Chr") for c in chroms]); ax.set_ylabel("Mean 1-Mb density (%)")
        ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True); ax.legend(fontsize=5.4, loc="upper right")
        ax.text(0.01, -0.25, "PZI A1/A2/B1/B2; full window heatmap exported separately.", transform=ax.transAxes, fontsize=5.0, color="#555555")
    else:
        ax.text(0.5, 0.5, "Run 1-Mb TE/TRF density first", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    panel_label(ax, "d")
    save(fig, out_dir, "figure_te_followup_mechanism")


def draw_logo(ax, members, label, baseline=0.0, scale=1.0):
    if not members:
        return
    sequences = [row["motif"] for row in members]
    length = len(sequences[0])
    counts = np.zeros((length, 4), dtype=float)
    letters = "ACGT"
    for row, seq in zip(members, sequences):
        if len(seq) != length:
            continue
        weight = max(1.0, float(row.get("array_count", 1)))
        for i, base in enumerate(seq):
            if base in letters:
                counts[i, letters.index(base)] += weight
    totals = counts.sum(axis=1)
    for i in range(length):
        if totals[i] == 0:
            continue
        freq = counts[i] / totals[i]
        entropy = -sum(p * math.log(p, 2) for p in freq if p > 0)
        info = max(0.0, 2.0 - entropy)
        order = np.argsort(freq)
        y = 0.0
        for j in order:
            if freq[j] <= 0:
                continue
            height = freq[j] * info
            ax.text(i + 1, baseline + y * scale, letters[j], fontsize=3.0 + 3.8 * height, ha="center", va="bottom",
                    color={"A": "#4C78A8", "C": "#2A9D8F", "G": "#E39C45", "T": "#C44E52"}[letters[j]],
                    family="DejaVu Sans", fontweight="bold")
            y += height
    ax.text(0.01, baseline + 1.95 * scale, label, fontsize=4.7, va="top")


def logo_information(sequences):
    letters = "ACGT"
    length = max((len(seq) for seq in sequences), default=0)
    values = []
    for pos in range(length):
        column = [seq[pos] for seq in sequences if pos < len(seq) and seq[pos] in letters]
        if not column:
            values.append(0.0)
            continue
        total = len(column)
        frequencies = [column.count(letter) / total for letter in letters]
        entropy = -sum(freq * math.log(freq, 2) for freq in frequencies if freq > 0)
        # Schneider et al. small-sample correction for a four-letter alphabet.
        correction = (len(letters) - 1) / (2.0 * math.log(2.0) * total)
        values.append(max(0.0, 2.0 - entropy - correction))
    return values


_LOGO_PATHS = {
    letter: TextPath((0, 0), letter, size=1,
                     prop=FontProperties(family="DejaVu Sans", weight="bold"))
    for letter in "ACGT"
}


def add_logo_letter(ax, letter, x, y, height, color):
    """Place a glyph whose geometric height is exactly its information height."""
    if height <= 0:
        return
    path = _LOGO_PATHS[letter]
    bounds = path.get_extents()
    width = 0.86
    scale_x = width / max(bounds.width, 1e-9)
    scale_y = height / max(bounds.height, 1e-9)
    transform = (
        transforms.Affine2D()
        .translate(-bounds.xmin, -bounds.ymin)
        .scale(scale_x, scale_y)
        .translate(x - width / 2.0, y)
        + ax.transData
    )
    ax.add_patch(PathPatch(path, transform=transform, facecolor=color,
                           edgecolor="none", clip_on=True))


def draw_true_logo(ax, sequences, window=80, start=0):
    """Draw an information-content logo from an actual MAFFT alignment."""
    if not sequences:
        return
    start = max(0, min(start, max(0, max(len(seq) for seq in sequences) - 1)))
    length = min(window, max(len(seq) for seq in sequences) - start)
    letters = "ACGT"
    colors = {"A": "#4C78A8", "C": "#2A9D8F", "G": "#E39C45", "T": "#C44E52"}
    for pos in range(length):
        absolute_pos = start + pos
        column = [seq[absolute_pos] for seq in sequences if absolute_pos < len(seq) and seq[absolute_pos] in letters]
        if not column:
            continue
        counts = {letter: column.count(letter) for letter in letters}
        total = sum(counts.values())
        frequencies = {letter: count / total for letter, count in counts.items() if count}
        entropy = -sum(freq * math.log(freq, 2) for freq in frequencies.values())
        # Small-sample correction is deliberately conservative; the logo is a
        # descriptive representation of the sampled array monomers.
        correction = (len(letters) - 1) / (2.0 * math.log(2.0) * total)
        information = max(0.0, 2.0 - entropy - correction)
        y = 0.0
        for letter, freq in sorted(frequencies.items(), key=lambda item: item[1]):
            height = freq * information
            add_logo_letter(ax, letter, pos + 1, y, height, colors[letter])
            y += height
    ax.set_xlim(0.5, length + 0.5)
    ax.set_ylim(0, 2.05)
    ax.set_yticks([0, 1, 2], ["0", "1", "2"], fontsize=4.5)
    ax.set_ylabel("bits", fontsize=5.0)
    ax.set_xticks(np.arange(0, length + 1, 10), [str(i) for i in range(0, length + 1, 10)], fontsize=4.5)
    ax.tick_params(axis="x", length=2, width=0.5)


def plot_supplement(stats_dir, out_dir):
    clusters = read_tsv(stats_dir / "trf_motif_clusters.tsv")[:8]
    members = read_tsv(stats_dir / "trf_motif_members.tsv")
    by_cluster = defaultdict(list)
    for row in members:
        by_cluster[row["cluster_id"]].append(row)
    orth = read_tsv(stats_dir / "orthosite_evidence_summary.tsv")
    nested = read_tsv(stats_dir / "te_nested_fragmentation.tsv") if (stats_dir / "te_nested_fragmentation.tsv").exists() else []

    fig = plt.figure(figsize=(7.25, 5.0))
    grid = fig.add_gridspec(2, 2, hspace=0.60, wspace=0.46)
    ax = fig.add_subplot(grid[0, 0])
    y = np.arange(len(clusters))
    vals = np.array([float(row["bp"]) / 1e6 for row in reversed(clusters)])
    ax.barh(y, vals, color="#C44E52", height=0.65)
    ax.set_yticks(y, [row["cluster_id"] for row in reversed(clusters)], fontsize=5.5)
    ax.set_xlabel("TRF span (Mb)"); ax.grid(axis="x", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True); panel_label(ax, "a")

    true_logo_summary = stats_dir / "trf_true_logo_summary.tsv"
    if true_logo_summary.exists():
        true_rows = read_tsv(true_logo_summary)[:4]
        logo_grid = grid[0, 1].subgridspec(max(1, len(true_rows)), 1, hspace=0.12)
        logo_axes = []
        for i, row in enumerate(true_rows):
            ax = fig.add_subplot(logo_grid[i, 0])
            alignment = read_fasta(stats_dir / "alignments" / row["alignment_fasta"])
            logo_sequences = [sequence for _, sequence in alignment]
            information = logo_information(logo_sequences)
            window = min(80, len(information))
            start = max(range(max(1, len(information) - window + 1)),
                        key=lambda i: sum(information[i:i + window])) if information else 0
            draw_true_logo(ax, logo_sequences, window=window, start=start)
            ax.set_ylabel(row["cluster_id"], fontsize=4.8, rotation=0, labelpad=18, va="center")
            if i == 0:
                ax.set_title("True logo from genome-extracted monomers", fontsize=7)
                panel_label(ax, "b")
            if i < len(true_rows) - 1:
                ax.set_xticklabels([])
            logo_axes.append(ax)
        logo_axes[-1].set_xlabel("Aligned monomer position (bp)", fontsize=5.5)
    else:
        ax = fig.add_subplot(grid[0, 1])
        ax.text(0.5, 0.5, "Run true-monomer extraction first", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        panel_label(ax, "b")

    ax = fig.add_subplot(grid[1, 0])
    by_target = defaultdict(lambda: defaultdict(int))
    for row in orth:
        by_target[row["target"]][row["evidence_class"]] += int(row["locus_count"])
    targets = sorted(by_target)
    bottom = np.zeros(len(targets))
    evidences = sorted({e for target in by_target for e in by_target[target]})
    evidence_labels = {
        "query_haplotype_specific_or_assembly": "assembly-linked",
        "unresolved": "unresolved",
        "shared_pzi_psp": "shared PZI–PSP",
        "pzi_specific_supported": "PZI-specific",
        "pzi_specific_no_outgroup": "PZI-specific (no outgroup)",
    }
    cmap = plt.get_cmap("tab20")
    for i, evidence in enumerate(evidences):
        vals = np.array([by_target[t].get(evidence, 0) for t in targets])
        ax.bar(np.arange(len(targets)), vals, bottom=bottom, color=cmap(i / max(1, len(evidences) - 1)), width=0.7, label=evidence_labels.get(evidence, evidence))
        bottom += vals
    ax.set_xticks(np.arange(len(targets)), targets, rotation=60, ha="right", fontsize=5)
    ax.set_ylabel("Orthosite loci"); ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
    ax.text(0.01, -0.32, "Current flank5k dataset contains PZI-A1 and PZI-B1 query sets only.", transform=ax.transAxes, fontsize=4.7, color="#555555")
    ax.legend(fontsize=4.5, ncol=2, loc="upper left", bbox_to_anchor=(0, 1.02)); panel_label(ax, "c")

    ax = fig.add_subplot(grid[1, 1])
    if nested:
        labels = [row["target"] for row in nested]
        vals = [100 * float(row["nested_overlap_fraction_of_raw_hits"]) for row in nested]
        ax.bar(np.arange(len(labels)), vals, color="#8E3B5B", width=0.68)
        ax.set_xticks(np.arange(len(labels)), labels, rotation=60, ha="right", fontsize=5)
        ax.set_ylabel("Pairwise nested overlap / raw hit span (%)"); ax.grid(axis="y", color="#E7E7E7", lw=0.5); ax.set_axisbelow(True)
        ax.text(0.01, -0.25, "Pairwise diagnostic; not a unique genome fraction.", transform=ax.transAxes, fontsize=4.7, color="#555555")
    else:
        ax.text(0.5, 0.5, "RepeatMasker .out required", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    panel_label(ax, "d")
    save(fig, out_dir, "supp_trf_motif_orthosite_nested")


def plot_density(stats_dir, out_dir):
    rows = read_tsv(stats_dir / "te_trf_1mb_density.tsv")
    targets = PZI_TARGETS
    subsets = {target: [row for row in rows if row["target"] == target] for target in targets}
    chrom_numbers = sorted({int(row["seqid"].split("_")[0].replace("chr", ""))
                            for subset in subsets.values() for row in subset})
    max_windows = {}
    for number in chrom_numbers:
        max_windows[number] = max(
            sum(1 for row in subsets[target] if int(row["seqid"].split("_")[0].replace("chr", "")) == number)
            for target in targets
        )
    offsets = {}
    total = 0
    for number in chrom_numbers:
        offsets[number] = total
        total += max_windows[number]
    te = np.full((len(targets), total), np.nan); trf = np.full((len(targets), total), np.nan)
    for row_idx, target in enumerate(targets):
        for row in subsets[target]:
            number = int(row["seqid"].split("_")[0].replace("chr", ""))
            idx = offsets[number] + int(int(row["start"]) / 1_000_000)
            if idx < total:
                te[row_idx, idx] = float(row["te_fraction"]) * 100
                trf[row_idx, idx] = float(row["trf_fraction"]) * 100
    fig, axes = plt.subplots(2, 1, figsize=(7.25, 3.6), sharex=True, gridspec_kw={"hspace": 0.10})
    images = []
    for ax, matrix, title, cmap in zip(axes, [te, trf], ["pan-TE", "TRF"], ["Blues", "Reds"]):
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad("#E6E6E6")
        im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap_obj, vmin=0, vmax=100)
        images.append(im)
        ax.set_yticks(np.arange(len(targets)), [f"{target}-{title}" for target in targets], fontsize=5.2)
        ax.tick_params(axis="y", length=0)
        for number in chrom_numbers[:-1]:
            ax.axvline(offsets[number] + max_windows[number] - 0.5, color="white", lw=0.7)
    centers = [offsets[number] + max_windows[number] / 2 - 0.5 for number in chrom_numbers]
    axes[-1].set_xticks(centers, [f"Chr{number:02d}" for number in chrom_numbers], fontsize=6)
    axes[-1].set_xlabel("PZI chromosome windows (1 Mb)")
    fig.colorbar(images[0], ax=axes[0], shrink=0.75, pad=0.02, label="pan-TE coverage (%)")
    fig.colorbar(images[1], ax=axes[1], shrink=0.75, pad=0.02, label="TRF coverage (%)")
    save(fig, out_dir, "figure_te_trf_1mb_density_pzi_four_haplotypes")


def plot_16chr_tracks(stats_dir, out_dir, a_target="pziA1", b_target="pziB1", stem="figure_pzi_16chr_te_trf_tracks"):
    """Draw explicit A01-A08 and B01-B08 chromosome tracks.

    Each chromosome has two aligned physical-position tracks: pan-TE coverage
    (blue) and TRF coverage (red).  The default A1/B1 pair is the representative
    diploidized A/B map; A2/B2 can be rendered with the same function as a
    sensitivity/phase comparison.
    """
    rows = read_tsv(stats_dir / "te_trf_1mb_density.tsv")
    by_target = defaultdict(list)
    for row in rows:
        if row["target"] in {a_target, b_target}:
            by_target[row["target"]].append(row)
    if not by_target.get(a_target) or not by_target.get(b_target):
        return

    tracks = []
    for label, target in (("A", a_target), ("B", b_target)):
        by_chr = defaultdict(list)
        for row in by_target[target]:
            by_chr[row["seqid"]].append(row)
        chromosomes = sorted(by_chr, key=lambda x: int(x.split("_")[0].replace("chr", "")))
        for number, seqid in enumerate(chromosomes, 1):
            tracks.append((f"{label}{number:02d}", by_chr[seqid]))

    fig, ax = plt.subplots(figsize=(7.25, 7.0))
    blue = plt.get_cmap("Blues").copy(); blue.set_bad("#F0F0F0")
    red = plt.get_cmap("Reds").copy(); red.set_bad("#F0F0F0")
    norm = mpl.colors.Normalize(vmin=0, vmax=100)
    row_gap = 1.15
    for i, (label, chrom_rows) in enumerate(tracks):
        y_te = len(tracks) - 1 - i
        y_trf = y_te - 0.34
        max_end = max(float(row["end"]) for row in chrom_rows) / 1e6
        for row in chrom_rows:
            start = float(row["start"]) / 1e6
            width = (float(row["end"]) - float(row["start"])) / 1e6
            te_fraction = float(row["te_fraction"]) * 100
            trf_fraction = float(row["trf_fraction"]) * 100
            ax.barh(y_te, width, left=start, height=0.27, color=blue(norm(te_fraction)), linewidth=0)
            ax.barh(y_trf, width, left=start, height=0.27, color=red(norm(trf_fraction)), linewidth=0)
        ax.text(-0.025, (y_te + y_trf) / 2, label, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=6.5, fontweight="bold")
        ax.text(max_end + 1.2, y_te, f"{max_end:.0f} Mb", ha="left", va="center", fontsize=5.0, color="#555555")
        if label == "A08":
            ax.axhline(y_trf - 0.50, color="#444444", lw=0.9)

    ax.set_yticks([])
    ax.set_xlabel("Physical position along chromosome (Mb)")
    ax.set_xlim(0, max(float(row["end"]) for row in rows if row["target"] in {a_target, b_target}) / 1e6 + 12)
    ax.set_ylim(-1.0, len(tracks) + 0.7)
    ax.grid(axis="x", color="#E7E7E7", lw=0.5)
    ax.set_axisbelow(True)
    ax.text(0.0, 1.02, f"A: {a_target}   B: {b_target}", transform=ax.transAxes, fontsize=6.0, color="#555555")
    sm_te = mpl.cm.ScalarMappable(norm=norm, cmap=blue); sm_te.set_array([])
    sm_trf = mpl.cm.ScalarMappable(norm=norm, cmap=red); sm_trf.set_array([])
    cbar_te = fig.colorbar(sm_te, ax=ax, fraction=0.025, pad=0.025)
    cbar_te.set_label("pan-TE coverage per 1-Mb window (%)", fontsize=6)
    cbar_trf = fig.colorbar(sm_trf, ax=ax, fraction=0.025, pad=0.09)
    cbar_trf.set_label("TRF coverage per 1-Mb window (%)", fontsize=6)
    cbar_te.ax.tick_params(labelsize=5); cbar_trf.ax.tick_params(labelsize=5)
    panel_label(ax, "a")
    save(fig, out_dir, stem)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    plot_mechanism(args.stats_dir, args.output_dir)
    plot_supplement(args.stats_dir, args.output_dir)
    if (args.stats_dir / "te_trf_1mb_density.tsv").exists():
        plot_density(args.stats_dir, args.output_dir)
        plot_16chr_tracks(args.stats_dir, args.output_dir, "pziA1", "pziB1",
                          "figure_pzi_16chr_te_trf_tracks_A1B1")
        plot_16chr_tracks(args.stats_dir, args.output_dir, "pziA2", "pziB2",
                          "supp_pzi_16chr_te_trf_tracks_A2B2")
    print(f"Figures written to {args.output_dir}")


if __name__ == "__main__":
    main()
