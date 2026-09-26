#!/usr/bin/env python3
"""Helpers for validating PZI assembly-length expansion with ONT evidence.

The script deliberately labels sequence lacking a whole-genome alignment to
the supplied comparators as "non-aligned candidate sequence". It does not
infer lineage-specific insertion events; that requires orthologous flanks and
an occupied/empty-site analysis.
"""

import argparse
import bisect
import collections
import csv
import statistics
import subprocess
from pathlib import Path


def read_fai(path):
    lengths = collections.OrderedDict()
    with open(path) as handle:
        for line in handle:
            if line.strip():
                fields = line.rstrip("\n").split("\t")
                lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise SystemExit("ERROR: no sequences found in FAI: {}".format(path))
    return lengths


def merge_intervals(intervals, max_gap=0):
    if not intervals:
        return []
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1] + max_gap:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def parse_paf_candidates(args):
    lengths = read_fai(args.fai)
    per_reference = {}
    union_intervals = {seq_id: [] for seq_id in lengths}

    for label_path in args.paf:
        try:
            label, paf_path = label_path.split(":", 1)
        except ValueError:
            raise SystemExit("ERROR: --paf must be label:path, observed {}".format(label_path))
        intervals = {seq_id: [] for seq_id in lengths}
        accepted = 0
        with open(paf_path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                target = fields[5]
                if target not in lengths:
                    continue
                block_len = int(fields[10])
                mapq = int(fields[11])
                if block_len < args.min_paf_block or mapq < args.min_mapq:
                    continue
                start = max(0, int(fields[7]))
                end = min(lengths[target], int(fields[8]))
                if end <= start:
                    continue
                intervals[target].append((start, end))
                union_intervals[target].append((start, end))
                accepted += 1
        per_reference[label] = {
            seq_id: merge_intervals(values, args.merge_gap)
            for seq_id, values in intervals.items()
        }
        print("Accepted PAF blocks for {}: {}".format(label, accepted))

    union_merged = {
        seq_id: merge_intervals(values, args.merge_gap)
        for seq_id, values in union_intervals.items()
    }

    coverage_rows = []
    for label, by_seq in per_reference.items():
        for seq_id, length in lengths.items():
            covered = sum(end - start for start, end in by_seq[seq_id])
            coverage_rows.append([
                label, seq_id, length, covered,
                "{:.6f}".format(covered / length if length else 0),
                len(by_seq[seq_id]),
            ])
    write_tsv(
        args.coverage_tsv,
        ["comparator", "pzi_seq", "pzi_seq_bp", "aligned_bp", "aligned_fraction", "merged_blocks"],
        coverage_rows,
    )

    candidates = []
    for seq_id, length in lengths.items():
        cursor = 0
        for start, end in union_merged[seq_id]:
            if start - cursor >= args.min_candidate_len:
                candidates.append((seq_id, cursor, start))
            cursor = max(cursor, end)
        if length - cursor >= args.min_candidate_len:
            candidates.append((seq_id, cursor, length))

    with open(args.candidates_bed, "w") as handle:
        for seq_id, start, end in candidates:
            handle.write("{}\t{}\t{}\t{}\n".format(seq_id, start, end, "{}_{}_{}".format(seq_id, start, end)))

    candidate_bp = sum(end - start for _, start, end in candidates)
    write_tsv(
        args.candidate_summary_tsv,
        ["metric", "value"],
        [
            ["pzi_total_bp", sum(lengths.values())],
            ["non_aligned_candidate_count", len(candidates)],
            ["non_aligned_candidate_bp", candidate_bp],
            ["non_aligned_candidate_fraction", "{:.6f}".format(candidate_bp / sum(lengths.values()))],
            ["definition", "complement_of_union_of_filtered_whole_genome_alignment_blocks"],
            ["interpretation", "candidate_sequence_only_not_lineage_specific_insertions"],
        ],
    )


def depth_windows(args):
    lengths = read_fai(args.fai)
    window = args.window_size
    out_rows = []
    summary = []
    current_seq = None
    current_window_start = 0
    depth_sum = 0
    bases = 0
    covered = 0

    def finish_window(seq_id, start, observed_bases, observed_sum, observed_covered):
        end = min(start + window, lengths[seq_id])
        expected = end - start
        missing = expected - observed_bases
        total_sum = observed_sum
        total_covered = observed_covered
        if missing < 0:
            raise SystemExit("ERROR: depth input exceeds expected window length")
        out_rows.append([
            seq_id, start, end, expected, total_sum,
            "{:.6f}".format(total_sum / expected if expected else 0),
            total_covered,
            "{:.6f}".format(total_covered / expected if expected else 0),
        ])

    def finish_sequence(seq_id, start, observed_bases, observed_sum, observed_covered):
        if seq_id is None:
            return
        finish_window(seq_id, start, observed_bases, observed_sum, observed_covered)
        next_start = start + window
        while next_start < lengths[seq_id]:
            finish_window(seq_id, next_start, 0, 0, 0)
            next_start += window

    for line in args.depth_input:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        seq_id, one_based_pos, depth = fields[0], int(fields[1]), int(fields[2])
        if seq_id not in lengths:
            continue
        pos = one_based_pos - 1
        if current_seq != seq_id:
            finish_sequence(current_seq, current_window_start, bases, depth_sum, covered)
            current_seq = seq_id
            current_window_start = 0
            depth_sum = 0
            bases = 0
            covered = 0
        while pos >= current_window_start + window:
            finish_window(current_seq, current_window_start, bases, depth_sum, covered)
            current_window_start += window
            depth_sum = 0
            bases = 0
            covered = 0
        depth_sum += depth
        bases += 1
        if depth > 0:
            covered += 1
    finish_sequence(current_seq, current_window_start, bases, depth_sum, covered)

    by_seq = collections.defaultdict(lambda: [0, 0.0, 0, 0])
    for row in out_rows:
        seq_id, _, _, bp, depth_sum_value, _, covered_bp, _ = row
        value = by_seq[seq_id]
        value[0] += int(bp)
        value[1] += int(depth_sum_value)
        value[2] += int(covered_bp)
        value[3] += 1
    for seq_id in lengths:
        bp, total_depth, covered_bp, n_windows = by_seq[seq_id]
        summary.append([
            seq_id, bp, "{:.6f}".format(total_depth / bp if bp else 0),
            covered_bp, "{:.6f}".format(covered_bp / bp if bp else 0), n_windows,
        ])

    write_tsv(
        args.windows_tsv,
        ["seq_id", "start", "end", "window_bp", "depth_sum", "mean_depth", "covered_bp", "covered_fraction"],
        out_rows,
    )
    write_tsv(
        args.summary_tsv,
        ["seq_id", "bp", "mean_depth", "covered_bp", "covered_fraction", "window_count"],
        summary,
    )


def cigar_ref_span(cigar):
    number = ""
    span = 0
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        if not number:
            return 0
        length = int(number)
        if char in "MDN=X":
            span += length
        number = ""
    return span if not number else 0


def read_bed(path):
    by_seq = collections.defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            start, end = int(fields[1]), int(fields[2])
            if end > start:
                by_seq[fields[0]].append((start, end, fields[3] if len(fields) > 3 else "{}_{}_{}".format(fields[0], start, end)))
    return {seq_id: sorted(values) for seq_id, values in by_seq.items()}


def load_window_depth(path):
    values = collections.defaultdict(list)
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            values[row["seq_id"]].append((int(row["start"]), int(row["end"]), float(row["mean_depth"])))
    return values


def candidate_support(args):
    candidates = read_bed(args.candidates_bed)
    depth = load_window_depth(args.depth_windows)
    all_depth = [item[2] for rows in depth.values() for item in rows if item[2] > 0]
    if not all_depth:
        raise SystemExit("ERROR: no positive-depth windows in {}".format(args.depth_windows))
    global_median = statistics.median(all_depth)

    boundaries = {}
    events = {}
    for seq_id, rows in candidates.items():
        events[seq_id] = []
        for index, (start, end, name) in enumerate(rows):
            events[seq_id].append((start, index, "left"))
            events[seq_id].append((end, index, "right"))
        boundaries[seq_id] = sorted(events[seq_id])
    positions = {seq_id: [item[0] for item in rows] for seq_id, rows in boundaries.items()}
    support = {(seq_id, index, side): set() for seq_id, rows in candidates.items() for index in range(len(rows)) for side in ("left", "right")}

    command = ["samtools", "view", "-F", str(args.exclude_flags), "-q", str(args.min_mapq), args.bam]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        seq_id, start = fields[2], int(fields[3]) - 1
        if seq_id not in boundaries or start < 0:
            continue
        span = cigar_ref_span(fields[5])
        if span <= 0:
            continue
        end = start + span
        left = bisect.bisect_left(positions[seq_id], start + args.flank_bp)
        right = bisect.bisect_right(positions[seq_id], end - args.flank_bp)
        for boundary, index, side in boundaries[seq_id][left:right]:
            support[(seq_id, index, side)].add(fields[0])
    if process.wait() != 0:
        raise SystemExit("ERROR: samtools view failed for {}".format(args.bam))

    rows = []
    accepted = []
    for seq_id, intervals in candidates.items():
        for index, (start, end, name) in enumerate(intervals):
            overlapping_depth = [value for win_start, win_end, value in depth.get(seq_id, []) if win_end > start and win_start < end]
            region_depth = statistics.mean(overlapping_depth) if overlapping_depth else 0
            ratio = region_depth / global_median if global_median else 0
            left_count = len(support[(seq_id, index, "left")])
            right_count = len(support[(seq_id, index, "right")])
            both_count = len(support[(seq_id, index, "left")] & support[(seq_id, index, "right")])
            row = [
                name, seq_id, start, end, end - start,
                left_count, right_count, both_count,
                "{:.6f}".format(region_depth),
                "{:.6f}".format(global_median),
                "{:.6f}".format(ratio),
                "yes" if (both_count >= args.min_spanning_reads and args.min_depth_ratio <= ratio <= args.max_depth_ratio) else "no",
            ]
            rows.append(row)
            if row[-1] == "yes":
                accepted.append(row)
    header = [
        "candidate_id", "seq_id", "start", "end", "candidate_bp",
        "left_boundary_spanning_reads", "right_boundary_spanning_reads", "both_boundaries_spanning_reads",
        "mean_overlapping_window_depth", "global_window_depth_median", "depth_ratio_to_global_median", "passes_ont_support",
    ]
    write_tsv(args.output_tsv, header, rows)
    write_tsv(args.passing_tsv, header, accepted)
    write_tsv(
        args.summary_tsv,
        ["metric", "value"],
        [
            ["candidate_count", len(rows)],
            ["candidate_bp", sum(int(row[4]) for row in rows)],
            ["same_read_crosses_both_boundaries_and_depth_pass_count", len(accepted)],
            ["same_read_crosses_both_boundaries_and_depth_pass_bp", sum(int(row[4]) for row in accepted)],
            ["min_spanning_reads_per_boundary", args.min_spanning_reads],
            ["flank_bp", args.flank_bp],
            ["mapq", args.min_mapq],
            ["caveat", "shared_tetraploid_ONT_reads_can_cross_map_unless_haplotype_binned_reads_are_used"],
        ],
    )


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_visualizations(args):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.patches import Patch
    except ImportError as exc:
        raise SystemExit("ERROR: matplotlib is required for visualization: {}".format(exc))

    windows = read_tsv(args.depth_windows)
    candidates = read_bed(args.candidates_bed)
    coverage = read_tsv(args.coverage_tsv)
    support = read_tsv(args.support_tsv)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    palette = {"pass": "#2b7a78", "fail": "#b2473e", "candidate": "#e8c45b", "depth": "#34495e"}
    by_seq = collections.defaultdict(list)
    for row in windows:
        by_seq[row["seq_id"]].append(row)

    depth_pdf = str(prefix) + ".ont_depth_candidates.pdf"
    depth_png = str(prefix) + ".ont_depth_candidates.png"
    ordered_seqs = list(by_seq)
    with PdfPages(depth_pdf) as pdf:
        ncols = 2
        nrows = 2
        for page_start in range(0, len(ordered_seqs), ncols * nrows):
            seqs = ordered_seqs[page_start:page_start + ncols * nrows]
            fig, axes = plt.subplots(nrows, ncols, figsize=(11.5, 7.5), squeeze=False)
            for ax in axes.flat:
                ax.set_visible(False)
            for ax, seq_id in zip(axes.flat, seqs):
                ax.set_visible(True)
                rows = by_seq[seq_id]
                x = [(int(row["start"]) + int(row["end"])) / 2 / 1e6 for row in rows]
                y = [float(row["mean_depth"]) for row in rows]
                ax.plot(x, y, color=palette["depth"], linewidth=0.7)
                for start, end, _ in candidates.get(seq_id, []):
                    ax.axvspan(start / 1e6, end / 1e6, color=palette["candidate"], alpha=0.45, linewidth=0)
                ax.set_title(seq_id, fontsize=9)
                ax.set_xlabel("Position (Mb)", fontsize=8)
                ax.set_ylabel("ONT depth", fontsize=8)
                ax.tick_params(labelsize=7)
                ax.spines[["top", "right"]].set_visible(False)
            fig.suptitle(args.title + " - ONT depth and non-aligned candidate intervals", fontsize=11)
            fig.legend([Patch(facecolor=palette["candidate"], alpha=0.45)], ["Non-aligned candidate interval"], loc="lower center", ncol=1, frameon=False, fontsize=8)
            fig.tight_layout(rect=(0, 0.06, 1, 0.94))
            pdf.savefig(fig, dpi=300)
            if page_start == 0:
                fig.savefig(depth_png, dpi=300)
            plt.close(fig)

    comparators = sorted({row["comparator"] for row in coverage})
    sequences = sorted({row["pzi_seq"] for row in coverage})
    matrix = []
    lookup = {(row["comparator"], row["pzi_seq"]): float(row["aligned_fraction"]) for row in coverage}
    for comparator in comparators:
        matrix.append([lookup.get((comparator, seq_id), 0.0) for seq_id in sequences])
    fig_width = max(8, len(sequences) * 0.65)
    fig, ax = plt.subplots(figsize=(fig_width, max(3.5, len(comparators) * 0.7 + 1.5)))
    image = ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(sequences)))
    ax.set_xticklabels(sequences, rotation=90, fontsize=7)
    ax.set_yticks(range(len(comparators)))
    ax.set_yticklabels(comparators, fontsize=8)
    ax.set_xlabel("PZI chromosomes")
    ax.set_title(args.title + " - comparator alignment coverage")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    colorbar.set_label("Aligned fraction of PZI chromosome", fontsize=8)
    fig.tight_layout()
    fig.savefig(str(prefix) + ".comparator_coverage_heatmap.png", dpi=300)
    fig.savefig(str(prefix) + ".comparator_coverage_heatmap.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for label, color in [("yes", palette["pass"]), ("no", palette["fail"])]:
        rows = [row for row in support if row["passes_ont_support"] == label]
        ax.scatter(
            [int(row["candidate_bp"]) / 1e3 for row in rows],
            [float(row["depth_ratio_to_global_median"]) for row in rows],
            s=[max(12, min(90, int(row["both_boundaries_spanning_reads"]) * 6)) for row in rows],
            alpha=0.7,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            label="ONT-supported" if label == "yes" else "Not supported",
        )
    ax.axhspan(args.min_depth_ratio, args.max_depth_ratio, color="#d9e8e7", alpha=0.65, zorder=0)
    ax.axhline(1.0, color="#68737d", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Candidate interval length (kb)")
    ax.set_ylabel("Depth / genome-wide median depth")
    ax.set_title(args.title + " - ONT support for candidate intervals")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(str(prefix) + ".candidate_ont_support.png", dpi=300)
    fig.savefig(str(prefix) + ".candidate_ont_support.pdf")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    candidates = subparsers.add_parser("candidates", help="Call non-aligned PZI candidate intervals from PAF alignments.")
    candidates.add_argument("--fai", required=True)
    candidates.add_argument("--paf", action="append", required=True, help="Comparator label and PAF as label:path. PZI must be the PAF target.")
    candidates.add_argument("--candidates-bed", required=True)
    candidates.add_argument("--coverage-tsv", required=True)
    candidates.add_argument("--candidate-summary-tsv", required=True)
    candidates.add_argument("--min-paf-block", type=int, default=5000)
    candidates.add_argument("--min-mapq", type=int, default=20)
    candidates.add_argument("--merge-gap", type=int, default=2000)
    candidates.add_argument("--min-candidate-len", type=int, default=10000)
    candidates.set_defaults(func=parse_paf_candidates)

    depth = subparsers.add_parser("depth", help="Summarize samtools depth -aa output in fixed windows.")
    depth.add_argument("--fai", required=True)
    depth.add_argument("--depth-input", type=argparse.FileType("r"), default="-")
    depth.add_argument("--windows-tsv", required=True)
    depth.add_argument("--summary-tsv", required=True)
    depth.add_argument("--window-size", type=int, default=100000)
    depth.set_defaults(func=depth_windows)

    support = subparsers.add_parser("support", help="Count ONT primary alignments crossing candidate interval boundaries.")
    support.add_argument("--bam", required=True)
    support.add_argument("--candidates-bed", required=True)
    support.add_argument("--depth-windows", required=True)
    support.add_argument("--output-tsv", required=True)
    support.add_argument("--passing-tsv", required=True)
    support.add_argument("--summary-tsv", required=True)
    support.add_argument("--flank-bp", type=int, default=1000)
    support.add_argument("--min-spanning-reads", type=int, default=3)
    support.add_argument("--min-mapq", type=int, default=20)
    support.add_argument("--exclude-flags", type=int, default=2308, help="Unmapped + secondary + supplementary SAM flags.")
    support.add_argument("--min-depth-ratio", type=float, default=0.5)
    support.add_argument("--max-depth-ratio", type=float, default=1.5)
    support.set_defaults(func=candidate_support)

    visualize = subparsers.add_parser("visualize", help="Create QC plots from depth, candidate and ONT-support tables.")
    visualize.add_argument("--depth-windows", required=True)
    visualize.add_argument("--candidates-bed", required=True)
    visualize.add_argument("--coverage-tsv", required=True)
    visualize.add_argument("--support-tsv", required=True)
    visualize.add_argument("--out-prefix", required=True)
    visualize.add_argument("--title", required=True)
    visualize.add_argument("--min-depth-ratio", type=float, default=0.5)
    visualize.add_argument("--max-depth-ratio", type=float, default=1.5)
    visualize.set_defaults(func=make_visualizations)
    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
