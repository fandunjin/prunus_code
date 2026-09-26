#!/usr/bin/env python3
"""Audit high-contribution pan-TE components against tandem-repeat evidence."""

import argparse
import bisect
import collections
import csv
import math
import statistics
from pathlib import Path


CORE_TARGETS = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
]
GROUP_TARGETS = {
    "PZI": CORE_TARGETS[:4],
    "PSP": CORE_TARGETS[4:],
}


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0].split("#", 1)[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks)
    return records


def selected_driver_rows(path, max_components, cumulative_threshold):
    rows = [row for row in read_tsv(path) if int(row.get("positive_excess_bp") or 0) > 0]
    rows.sort(key=lambda row: int(row["positive_rank"]))
    selected = []
    for row in rows:
        if len(selected) >= max_components:
            break
        selected.append(row)
        if float(row["cumulative_positive_excess_fraction"]) >= cumulative_threshold:
            break
    return selected


def parse_tesorter(path):
    if not path or not path.is_file() or path.stat().st_size == 0:
        return {}
    result = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fields = line.split()
            if len(fields) < 4:
                continue
            family = fields[0].split("#", 1)[0]
            result[family] = {
                "tesorter_order": fields[1],
                "tesorter_superfamily": fields[2],
                "tesorter_lineage": fields[3],
                "tesorter_complete": fields[4] if len(fields) > 4 else "",
                "tesorter_domains": fields[6] if len(fields) > 6 else "",
            }
    return result


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [(start, end) for start, end in merged]


def interval_length(intervals):
    return sum(end - start for start, end in intervals)


def read_tandem_bed(path):
    result = collections.defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            result[fields[0]].append((start, end))
    return {seqid: merge_intervals(values) for seqid, values in result.items()}


def build_interval_index(intervals_by_seq):
    index = {}
    for seqid, intervals in intervals_by_seq.items():
        starts = [start for start, _ in intervals]
        index[seqid] = (intervals, starts)
    return index


def overlap_length(start, end, seq_index):
    if not seq_index:
        return 0
    intervals, starts = seq_index
    cursor = bisect.bisect_left(starts, end) - 1
    overlap = 0
    while cursor >= 0 and intervals[cursor][1] > start:
        left, right = intervals[cursor]
        overlap += max(0, min(end, right) - max(start, left))
        cursor -= 1
    return overlap


def parse_rmout(path, selected, library_lengths, full_hit_fraction):
    by_family_seq = collections.defaultdict(lambda: collections.defaultdict(list))
    hit_lengths = collections.defaultdict(list)
    divergences = collections.defaultdict(list)
    full_hits = collections.Counter()
    hit_counts = collections.Counter()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            family = fields[9]
            if family not in selected:
                continue
            try:
                start = int(fields[5]) - 1
                end = int(fields[6])
                divergence = float(fields[1])
            except ValueError:
                continue
            if end <= start:
                continue
            length = end - start
            by_family_seq[family][fields[4]].append((start, end))
            hit_lengths[family].append(length)
            divergences[family].append(divergence)
            hit_counts[family] += 1
            if length >= library_lengths.get(family, math.inf) * full_hit_fraction:
                full_hits[family] += 1
    merged = {
        family: {seqid: merge_intervals(intervals) for seqid, intervals in seqs.items()}
        for family, seqs in by_family_seq.items()
    }
    return merged, hit_lengths, divergences, full_hits, hit_counts


def summarize_target(target, rmout, tandem_bed, selected, lengths, full_hit_fraction):
    tandem_index = build_interval_index(read_tandem_bed(tandem_bed))
    merged, hit_lengths, divergences, full_hits, hit_counts = parse_rmout(
        rmout, selected, lengths, full_hit_fraction
    )
    rows = []
    for family in selected:
        seqs = merged.get(family, {})
        union_bp = sum(interval_length(values) for values in seqs.values())
        tandem_bp = sum(
            overlap_length(start, end, tandem_index.get(seqid))
            for seqid, intervals in seqs.items()
            for start, end in intervals
        )
        chromosome_bp = [interval_length(values) for values in seqs.values()]
        rows.append({
            "target": target,
            "pan_family": family,
            "hit_count": hit_counts[family],
            "union_bp": union_bp,
            "tandem_overlap_bp": tandem_bp,
            "tandem_overlap_fraction": f"{tandem_bp / union_bp:.8f}" if union_bp else "",
            "median_hit_length_bp": (
                f"{statistics.median(hit_lengths[family]):.2f}" if hit_lengths[family] else ""
            ),
            "full_length_hit_fraction": (
                f"{full_hits[family] / hit_counts[family]:.8f}" if hit_counts[family] else ""
            ),
            "median_repeatmasker_divergence_pct": (
                f"{statistics.median(divergences[family]):.6f}" if divergences[family] else ""
            ),
            "largest_sequence_bp_fraction": (
                f"{max(chromosome_bp) / union_bp:.8f}" if union_bp and chromosome_bp else ""
            ),
        })
    return rows


def numeric_mean(rows, field):
    values = [float(row[field]) for row in rows if row.get(field) not in (None, "")]
    return sum(values) / len(values) if values else math.nan


def command_select(args):
    drivers = selected_driver_rows(args.driver_table, args.max_components, args.cumulative_threshold)
    sequences = read_fasta(args.pan_library)
    missing = [row["pan_family"] for row in drivers if row["pan_family"] not in sequences]
    if missing:
        raise SystemExit("ERROR: selected components missing from pan library: " + ", ".join(missing))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(drivers[0]) if drivers else []
    write_tsv(args.output_dir / "selected_components.tsv", fields, drivers)
    with open(args.output_dir / "selected_components.fa", "w", encoding="utf-8") as handle:
        for row in drivers:
            family = row["pan_family"]
            handle.write(f">{family}#{row['class_family']}\n")
            sequence = sequences[family]
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")
    print(f"Selected {len(drivers)} pan-TE components")


def command_trf_to_bed(args):
    sequence = None
    rows = []
    with open(args.input, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Sequence:"):
                sequence = line.split(":", 1)[1].strip().split()[0]
                continue
            if line.startswith("@"):
                sequence = line[1:].strip().split()[0]
                continue
            fields = line.split()
            if sequence is None or len(fields) < 14 or not fields[0].isdigit():
                continue
            try:
                start = int(fields[0]) - 1
                end = int(fields[1])
                period = int(fields[2])
                copies = float(fields[3])
                consensus = int(fields[4])
                score = int(fields[7])
                entropy = float(fields[12])
            except ValueError:
                continue
            if end <= start or score < args.min_score:
                continue
            rows.append({
                "seqid": sequence,
                "start": start,
                "end": end,
                "period_bp": period,
                "copy_number": f"{copies:.3f}",
                "consensus_size_bp": consensus,
                "score": score,
                "entropy": f"{entropy:.3f}",
                # Field 14 in TRF -ngs output is the consensus repeat unit.
                # It is retained for satellite/minisatellite follow-up analysis.
                "consensus_motif": fields[13].upper() if len(fields) > 13 else "",
            })
    write_tsv(
        args.output,
        ["seqid", "start", "end", "period_bp", "copy_number", "consensus_size_bp", "score", "entropy", "consensus_motif"],
        rows,
    )
    print(f"Parsed {len(rows)} TRF intervals")


def command_audit(args):
    drivers = read_tsv(args.selected_components)
    selected = [row["pan_family"] for row in drivers]
    driver_by_family = {row["pan_family"]: row for row in drivers}
    lengths = {row["pan_family"]: int(row["library_sequence_length_bp"]) for row in drivers}
    tesorter = parse_tesorter(args.tesorter_cls)
    target_rows = []
    missing = []
    for target in CORE_TARGETS:
        rmout = args.mask_root / target / f"{target}.fa.out"
        tandem_bed = args.tandem_root / target / f"{target}.trf.tsv"
        if not rmout.is_file() or not tandem_bed.is_file():
            missing.append(target)
            continue
        target_rows.extend(summarize_target(
            target, rmout, tandem_bed, selected, lengths, args.full_hit_fraction
        ))
    if missing and not args.allow_missing_targets:
        raise SystemExit("ERROR: missing RepeatMasker/TRF inputs for: " + ", ".join(missing))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    by_family_target = {(row["pan_family"], row["target"]): row for row in target_rows}
    summary = []
    decisions = []
    for family in selected:
        driver = driver_by_family[family]
        record = dict(driver)
        for group, targets in GROUP_TARGETS.items():
            rows = [by_family_target[(family, target)] for target in targets if (family, target) in by_family_target]
            prefix = group.lower()
            for field in (
                "union_bp", "tandem_overlap_fraction", "full_length_hit_fraction",
                "median_hit_length_bp", "largest_sequence_bp_fraction",
            ):
                value = numeric_mean(rows, field)
                record[f"{prefix}_mean_{field}"] = "" if math.isnan(value) else f"{value:.8f}"
        record.update(tesorter.get(family, {
            "tesorter_order": "", "tesorter_superfamily": "", "tesorter_lineage": "",
            "tesorter_complete": "", "tesorter_domains": "",
        }))
        tandem_fraction = float(record.get("pzi_mean_tandem_overlap_fraction") or 0)
        low_complexity = int(driver.get("library_low_complexity_flag") or 0)
        short_high_copy = lengths[family] < args.short_component_bp and int(driver["pzi_mean_bp"]) >= args.high_copy_bp
        if low_complexity:
            status = "exclude_low_complexity"
        elif tandem_fraction >= args.exclude_tandem_fraction:
            status = "exclude_tandem_candidate"
        elif short_high_copy:
            status = "manual_structure_review"
        else:
            status = "pass_preliminary"
        record["short_high_copy_flag"] = int(short_high_copy)
        record["automated_status"] = status
        record["recommended_for_orthosite"] = int(status == "pass_preliminary")
        summary.append(record)
        decisions.append({
            "pan_family": family,
            "positive_rank": driver["positive_rank"],
            "class_family": driver["class_family"],
            "library_sequence_length_bp": lengths[family],
            "positive_excess_bp": driver["positive_excess_bp"],
            "pzi_mean_tandem_overlap_fraction": record.get("pzi_mean_tandem_overlap_fraction", ""),
            "automated_status": status,
            "manual_decision": "",
            "curated_family": "",
            "manual_evidence": "",
            "notes": "",
        })
    target_fields = [
        "target", "pan_family", "hit_count", "union_bp", "tandem_overlap_bp",
        "tandem_overlap_fraction", "median_hit_length_bp", "full_length_hit_fraction",
        "median_repeatmasker_divergence_pct", "largest_sequence_bp_fraction",
    ]
    summary_fields = list(summary[0]) if summary else []
    decision_fields = list(decisions[0]) if decisions else []
    write_tsv(args.output_dir / "component_audit_by_target.tsv", target_fields, target_rows)
    write_tsv(args.output_dir / "component_audit_summary.tsv", summary_fields, summary)
    write_tsv(args.output_dir / "component_decisions.tsv", decision_fields, decisions)
    print(f"Audited {len(summary)} components across {len(set(row['target'] for row in target_rows))} targets")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    select = subparsers.add_parser("select", help="select and extract high-contribution components")
    select.add_argument("--driver-table", required=True, type=Path)
    select.add_argument("--pan-library", required=True, type=Path)
    select.add_argument("--output-dir", required=True, type=Path)
    select.add_argument("--max-components", type=int, default=50)
    select.add_argument("--cumulative-threshold", type=float, default=0.60)
    select.set_defaults(func=command_select)

    trf = subparsers.add_parser("trf-to-bed", help="convert TRF -ngs output to a tabular BED-like file")
    trf.add_argument("--input", required=True, type=Path)
    trf.add_argument("--output", required=True, type=Path)
    trf.add_argument("--min-score", type=int, default=50)
    trf.set_defaults(func=command_trf_to_bed)

    audit = subparsers.add_parser("audit", help="summarize component and tandem-repeat overlap")
    audit.add_argument("--selected-components", required=True, type=Path)
    audit.add_argument("--mask-root", required=True, type=Path)
    audit.add_argument("--tandem-root", required=True, type=Path)
    audit.add_argument("--output-dir", required=True, type=Path)
    audit.add_argument("--tesorter-cls", type=Path)
    audit.add_argument("--full-hit-fraction", type=float, default=0.80)
    audit.add_argument("--exclude-tandem-fraction", type=float, default=0.50)
    audit.add_argument("--short-component-bp", type=int, default=500)
    audit.add_argument("--high-copy-bp", type=int, default=1_000_000)
    audit.add_argument("--allow-missing-targets", action="store_true")
    audit.set_defaults(func=command_audit)
    return parser


def main():
    args = build_parser().parse_args()
    if hasattr(args, "cumulative_threshold") and not 0 < args.cumulative_threshold <= 1:
        raise SystemExit("ERROR: --cumulative-threshold must be in (0, 1]")
    args.func(args)


if __name__ == "__main__":
    main()
