#!/usr/bin/env python3
"""Summarize non-redundant TRF tandem repeats across standardized genomes.

TRF calls overlap frequently.  Every base is therefore assigned once, to the
overlapping call with the highest TRF score, then longest span.  The resulting
``exclusive_union_bp`` is suitable for cross-genome composition comparisons;
raw call spans are retained separately for diagnostic purposes.
"""

import argparse
import collections
import csv
import gzip
import json
import math
import statistics
from pathlib import Path


TARGET_ORDER = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]
GROUP = {
    **{x: "PZI" for x in TARGET_ORDER[:4]},
    **{x: "PSP" for x in TARGET_ORDER[4:8]},
    "ppaA": "PPA", "ppaB": "PPA", "pcam": "PCAM", "ppe": "PPE",
}
PERIOD_CLASSES = ["1", "2-6", "7-20", "21-100", "101-500", "501-2000", ">2000"]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fasta_lengths(path):
    path = Path(path)
    fai = Path(str(path) + ".fai")
    if fai.is_file() and fai.stat().st_size:
        result = {}
        with open(fai, encoding="utf-8") as handle:
            for line in handle:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    result[f[0]] = int(f[1])
        if result:
            return result
    opener = gzip.open if str(path).endswith(".gz") else open
    result, name, size = {}, None, 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    result[name] = size
                name = line[1:].split()[0]
                size = 0
            else:
                size += len(line.strip())
    if name is not None:
        result[name] = size
    return result


def period_class(period):
    if period <= 1:
        return "1"
    if period <= 6:
        return "2-6"
    if period <= 20:
        return "7-20"
    if period <= 100:
        return "21-100"
    if period <= 500:
        return "101-500"
    if period <= 2000:
        return "501-2000"
    return ">2000"


def shortest_unit(seq):
    for width in range(1, len(seq) + 1):
        if len(seq) % width == 0 and seq == seq[:width] * (len(seq) // width):
            return seq[:width]
    return seq


def canonical_motif(value):
    motif = "".join(base for base in (value or "").upper() if base in "ACGT")
    if not motif:
        return "unknown"
    motif = shortest_unit(motif)
    # Long TRF consensus strings are not motifs in the satellite-DNA sense;
    # retaining their first 100 bp avoids quadratic rotation work and labels
    # them as long consensus units in output tables.
    if len(motif) > 100:
        return f"long_{len(motif)}bp:{motif[:30]}"
    complement = str.maketrans("ACGT", "TGCA")
    reverse = motif.translate(complement)[::-1]
    candidates = [motif[i:] + motif[:i] for i in range(len(motif))]
    candidates.extend(reverse[i:] + reverse[:i] for i in range(len(reverse)))
    return min(candidates)


def parse_trf(path, min_score=50):
    by_seq = collections.defaultdict(list)
    raw_by_class = collections.Counter()
    raw_span_by_class = collections.Counter()
    raw_count = 0
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"seqid", "start", "end", "period_bp", "score", "consensus_motif"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"TRF table has missing columns: {path}")
        for row in reader:
            try:
                start, end = int(row["start"]), int(row["end"])
                period, score = int(row["period_bp"]), int(row["score"])
                copies = float(row.get("copy_number", ""))
            except ValueError:
                continue
            if end <= start or score < min_score:
                continue
            cls = period_class(period)
            record = {
                "start": start, "end": end, "period_class": cls, "period_bp": period,
                "score": score, "copy_number": copies, "motif": canonical_motif(row.get("consensus_motif", "")),
                "raw_motif": row.get("consensus_motif", ""), "entropy": row.get("entropy", ""),
            }
            by_seq[row["seqid"]].append(record)
            raw_by_class[cls] += 1
            raw_span_by_class[cls] += end - start
            raw_count += 1
    return by_seq, raw_count, raw_by_class, raw_span_by_class


def resolve_sequence(records, genome_length, window_bp):
    """Allocate every overlapping TRF segment once and accumulate summaries."""
    events = collections.defaultdict(list)
    for index, row in enumerate(records):
        events[row["start"]].append((1, index))
        events[row["end"]].append((-1, index))
    active = set()
    class_bp, motif_bp = collections.Counter(), collections.Counter()
    total_bp = 0
    window_bp_counts = collections.Counter()
    last = None
    for position in sorted(events):
        if last is not None and position > last and active:
            best = min(active, key=lambda i: (-records[i]["score"], -(records[i]["end"] - records[i]["start"]), records[i]["motif"]))
            row, width = records[best], position - last
            total_bp += width
            class_bp[row["period_class"]] += width
            motif_bp[(row["period_class"], row["period_bp"], row["motif"])] += width
            left, right = last, position
            while left < right:
                window = left // window_bp
                bound = min(right, (window + 1) * window_bp)
                window_bp_counts[window] += bound - left
                left = bound
        for sign, index in events[position]:
            if sign == 1:
                active.add(index)
            else:
                active.discard(index)
        last = position
    return total_bp, class_bp, motif_bp, window_bp_counts


def target_rows(manifest, tandem_root, window_bp, top_motifs, top_arrays, min_score):
    genome_rows, class_rows, motif_by_target, top_array_rows, window_rows, qc_rows = [], [], {}, [], [], []
    for row in manifest:
        target = row["target"]
        genome = Path(row["mask_input_fasta"])
        trf = tandem_root / target / f"{target}.trf.tsv"
        if not genome.is_file() or not trf.is_file():
            raise FileNotFoundError(f"missing standardized genome/TRF table for {target}: {genome}, {trf}")
        sizes = fasta_lengths(genome)
        calls, raw_count, raw_by_class, raw_span_by_class = parse_trf(trf, min_score=min_score)
        genome_bp = sum(sizes.values())
        total_union, total_raw_span = 0, 0
        exclusive_classes, exclusive_motifs = collections.Counter(), collections.Counter()
        windows = {}
        for seqid, length in sizes.items():
            value, classes, motifs, win = resolve_sequence(calls.get(seqid, []), length, window_bp)
            total_union += value
            exclusive_classes.update(classes)
            exclusive_motifs.update(motifs)
            windows[seqid] = win
        total_raw_span = sum(raw_span_by_class.values())
        genome_rows.append({
            "target": target, "group": GROUP[target], "genome_bp": genome_bp,
            "trf_call_count": raw_count, "raw_call_span_bp": total_raw_span,
            "exclusive_union_bp": total_union,
            "genome_fraction": f"{total_union / genome_bp:.8f}" if genome_bp else "0",
            "raw_span_to_union_ratio": f"{total_raw_span / total_union:.6f}" if total_union else "",
        })
        for cls in PERIOD_CLASSES:
            class_rows.append({
                "target": target, "group": GROUP[target], "period_class": cls,
                "raw_call_count": raw_by_class[cls], "raw_call_span_bp": raw_span_by_class[cls],
                "exclusive_union_bp": exclusive_classes[cls],
                "genome_fraction": f"{exclusive_classes[cls] / genome_bp:.8f}" if genome_bp else "0",
            })
        ranked_motifs = sorted(exclusive_motifs.items(), key=lambda x: (-x[1], x[0]))[:top_motifs]
        motif_by_target[target] = ranked_motifs
        for seqid, values in windows.items():
            for number in range(math.ceil(sizes[seqid] / window_bp)):
                start = number * window_bp
                end = min(sizes[seqid], start + window_bp)
                bp = values[number]
                window_rows.append({
                    "target": target, "group": GROUP[target], "seqid": seqid,
                    "start": start, "end": end, "window_bp": end - start,
                    "tandem_union_bp": bp, "tandem_fraction": f"{bp / (end - start):.8f}" if end > start else "0",
                })
        arrays = []
        for seqid, values in calls.items():
            for item in values:
                arrays.append((item["end"] - item["start"], seqid, item))
        for length, seqid, item in sorted(arrays, key=lambda x: (-x[0], x[1], x[2]["start"]))[:top_arrays]:
            top_array_rows.append({
                "target": target, "group": GROUP[target], "seqid": seqid,
                "start": item["start"], "end": item["end"], "array_length_bp": length,
                "period_class": item["period_class"], "period_bp": item["period_bp"],
                "copy_number": f"{item['copy_number']:.3f}", "long_array_ge_10kb": int(length >= 10_000),
                "score": item["score"], "canonical_motif": item["motif"], "raw_consensus_motif": item["raw_motif"],
                "entropy": item["entropy"],
            })
        qc_rows.append({"target": target, "genome_contigs": len(sizes), "trf_sequences": len(calls), "genome_bp": genome_bp, "exclusive_union_bp": total_union, "class_sum_bp": sum(exclusive_classes.values()), "class_sum_matches_union": int(sum(exclusive_classes.values()) == total_union)})
    motif_rows = []
    for target, entries in motif_by_target.items():
        for (cls, period, motif), value in entries:
            motif_rows.append({"target": target, "group": GROUP[target], "period_class": cls, "period_bp": period, "canonical_motif": motif, "exclusive_union_bp": value})
    return genome_rows, class_rows, motif_rows, top_array_rows, window_rows, qc_rows


def group_summary(genome_rows, class_rows, motif_rows):
    genome_by_group = collections.defaultdict(list)
    for row in genome_rows:
        genome_by_group[row["group"]].append(row)
    group_genome = []
    for group in ("PZI", "PSP", "PPA", "PCAM", "PPE"):
        rows = genome_by_group[group]
        vals = [int(x["exclusive_union_bp"]) for x in rows]
        fracs = [float(x["genome_fraction"]) for x in rows]
        group_genome.append({"group": group, "genome_count": len(rows), "mean_tandem_union_bp": round(statistics.mean(vals)), "min_tandem_union_bp": min(vals), "max_tandem_union_bp": max(vals), "mean_genome_fraction": f"{statistics.mean(fracs):.8f}", "min_genome_fraction": f"{min(fracs):.8f}", "max_genome_fraction": f"{max(fracs):.8f}"})
    class_by_group = collections.defaultdict(list)
    for row in class_rows:
        class_by_group[(row["group"], row["period_class"])].append(int(row["exclusive_union_bp"]))
    group_class = []
    for group in ("PZI", "PSP", "PPA", "PCAM", "PPE"):
        for cls in PERIOD_CLASSES:
            vals = class_by_group[(group, cls)]
            group_class.append({"group": group, "period_class": cls, "genome_count": len(vals), "mean_exclusive_union_bp": round(statistics.mean(vals)), "min_exclusive_union_bp": min(vals), "max_exclusive_union_bp": max(vals)})
    motif_values = collections.defaultdict(dict)
    for row in motif_rows:
        motif_values[(row["period_class"], row["period_bp"], row["canonical_motif"])][row["target"]] = int(row["exclusive_union_bp"])
    group_motif = []
    groups = {group: [r["target"] for r in genome_by_group[group]] for group in genome_by_group}
    for (cls, period, motif), target_values in motif_values.items():
        out = {"period_class": cls, "period_bp": period, "canonical_motif": motif}
        for group, targets in groups.items():
            values = [target_values.get(target, 0) for target in targets]
            out[f"{group.lower()}_mean_bp"] = round(statistics.mean(values))
        out["pzi_minus_psp_mean_bp"] = out.get("pzi_mean_bp", 0) - out.get("psp_mean_bp", 0)
        group_motif.append(out)
    group_motif.sort(key=lambda x: x["pzi_minus_psp_mean_bp"], reverse=True)
    return group_genome, group_class, group_motif


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mask-manifest", type=Path, required=True)
    parser.add_argument("--tandem-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--window-bp", type=int, default=1_000_000)
    parser.add_argument("--top-motifs-per-target", type=int, default=100)
    parser.add_argument("--top-arrays-per-target", type=int, default=100)
    parser.add_argument("--min-score", type=int, default=50,
                        help="Retain TRF calls at or above this score; permits a no-rescan sensitivity analysis.")
    args = parser.parse_args()
    if args.window_bp <= 0:
        raise SystemExit("ERROR: --window-bp must be positive")
    manifest = read_tsv(args.mask_manifest)
    lookup = {row["target"]: row for row in manifest}
    missing = [target for target in TARGET_ORDER if target not in lookup]
    if missing:
        raise SystemExit("ERROR: mask manifest is missing " + ", ".join(missing))
    ordered = [lookup[target] for target in TARGET_ORDER]
    if args.min_score < 0:
        raise SystemExit("ERROR: --min-score must be non-negative")
    outputs = target_rows(
        ordered, args.tandem_root, args.window_bp, args.top_motifs_per_target,
        args.top_arrays_per_target, args.min_score,
    )
    genome_rows, class_rows, motif_rows, array_rows, window_rows, qc_rows = outputs
    group_genome, group_class, group_motif = group_summary(genome_rows, class_rows, motif_rows)
    fields = ["target", "group", "genome_bp", "trf_call_count", "raw_call_span_bp", "exclusive_union_bp", "genome_fraction", "raw_span_to_union_ratio"]
    write_tsv(args.output_dir / "tandem_genome_summary.tsv", fields, genome_rows)
    write_tsv(args.output_dir / "tandem_period_summary.tsv", ["target", "group", "period_class", "raw_call_count", "raw_call_span_bp", "exclusive_union_bp", "genome_fraction"], class_rows)
    write_tsv(args.output_dir / "tandem_top_motifs_by_target.tsv", ["target", "group", "period_class", "period_bp", "canonical_motif", "exclusive_union_bp"], motif_rows)
    write_tsv(args.output_dir / "tandem_top_arrays.tsv", ["target", "group", "seqid", "start", "end", "array_length_bp", "period_class", "period_bp", "copy_number", "long_array_ge_10kb", "score", "canonical_motif", "raw_consensus_motif", "entropy"], array_rows)
    write_tsv(args.output_dir / "tandem_window_density.tsv", ["target", "group", "seqid", "start", "end", "window_bp", "tandem_union_bp", "tandem_fraction"], window_rows)
    write_tsv(args.output_dir / "tandem_qc.tsv", ["target", "genome_contigs", "trf_sequences", "genome_bp", "exclusive_union_bp", "class_sum_bp", "class_sum_matches_union"], qc_rows)
    write_tsv(args.output_dir / "tandem_group_genome_summary.tsv", ["group", "genome_count", "mean_tandem_union_bp", "min_tandem_union_bp", "max_tandem_union_bp", "mean_genome_fraction", "min_genome_fraction", "max_genome_fraction"], group_genome)
    write_tsv(args.output_dir / "tandem_group_period_summary.tsv", ["group", "period_class", "genome_count", "mean_exclusive_union_bp", "min_exclusive_union_bp", "max_exclusive_union_bp"], group_class)
    write_tsv(args.output_dir / "tandem_group_motif_summary.tsv", ["period_class", "period_bp", "canonical_motif", "pzi_mean_bp", "psp_mean_bp", "ppa_mean_bp", "pcam_mean_bp", "ppe_mean_bp", "pzi_minus_psp_mean_bp"], group_motif)
    manifest = {"definition": "TRF intervals resolved to a non-redundant union by highest score then longest span", "period_classes": PERIOD_CLASSES, "window_bp": args.window_bp, "min_trf_score": args.min_score, "top_motifs_per_target": args.top_motifs_per_target, "warning": "TRF period classes are operational sequence classes, not automatic centromere/satellite annotations."}
    with open(args.output_dir / "tandem_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print(json.dumps({"targets": len(genome_rows), "union_bp_pzi_mean": round(statistics.mean(int(r["exclusive_union_bp"]) for r in genome_rows if r["group"] == "PZI")), "output_dir": str(args.output_dir)}))


if __name__ == "__main__":
    main()
