#!/usr/bin/env python3
"""Create mutually exclusive genome partitions for the tandem-repeat main figure.

The standardized pan-TE RepeatMasker and TRF tracks overlap extensively in
PZI.  This program assigns every base exactly once for the genome-size panel:

    TRF tandem array | pan-TE mask outside TRF | other sequence.

The partition describes sequence architecture, not a claim that all
TRF-overlapping sequence has a non-TE evolutionary origin.
"""

import argparse
import collections
import csv
import gzip
import json
import statistics
from pathlib import Path


TARGET_ORDER = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]
GROUPS = ["PZI", "PSP", "PPA", "PCAM", "PPE"]


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
    fai = Path(str(path) + ".fai")
    if fai.is_file() and fai.stat().st_size:
        result = {}
        with open(fai, encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    result[fields[0]] = int(fields[1])
        if result:
            return result
    opener = gzip.open if str(path).endswith(".gz") else open
    result, name, length = {}, None, 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    result[name] = length
                name, length = line[1:].split()[0], 0
            else:
                length += len(line.strip())
    if name is not None:
        result[name] = length
    return result


def merge(intervals):
    if not intervals:
        return []
    result = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not result or start > result[-1][1]:
            result.append([start, end])
        else:
            result[-1][1] = max(result[-1][1], end)
    return [(start, end) for start, end in result]


def interval_length(intervals):
    return sum(end - start for start, end in intervals)


def intersection_length(left, right):
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return total


def trf_union(path, min_score):
    by_seq = collections.defaultdict(list)
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                start, end = int(row["start"]), int(row["end"])
                score = int(row["score"])
            except (KeyError, ValueError):
                continue
            if end > start and score >= min_score:
                by_seq[row["seqid"]].append((start, end))
    return {seq: merge(values) for seq, values in by_seq.items()}


def repeatmasker_union(path):
    by_seq = collections.defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            try:
                start, end = int(fields[5]) - 1, int(fields[6])
            except ValueError:
                continue
            if end > start:
                by_seq[fields[4]].append((start, end))
    return {seq: merge(values) for seq, values in by_seq.items()}


def group_summary(rows):
    result = []
    for group in GROUPS:
        subset = [r for r in rows if r["group"] == group]
        if not subset:
            continue
        record = {"group": group, "genome_packages": len(subset)}
        for field in ("genome_bp", "trf_tandem_bp", "pan_te_non_tandem_bp", "other_bp", "repeatmasker_overlap_bp"):
            values = [int(r[field]) for r in subset]
            record[f"mean_{field}"] = round(statistics.mean(values))
            record[f"min_{field}"] = min(values)
            record[f"max_{field}"] = max(values)
        for field in ("trf_tandem_fraction", "pan_te_non_tandem_fraction", "other_fraction"):
            values = [float(r[field]) for r in subset]
            record[f"mean_{field}"] = f"{statistics.mean(values):.8f}"
            record[f"min_{field}"] = f"{min(values):.8f}"
            record[f"max_{field}"] = f"{max(values):.8f}"
        result.append(record)
    return result


def contrast_summary(group_rows, reference):
    lookup = {r["group"]: r for r in group_rows}
    if reference not in lookup:
        raise SystemExit(f"ERROR: reference group {reference} is absent")
    result = []
    ref = lookup[reference]
    for group in GROUPS:
        if group == reference or group not in lookup:
            continue
        other = lookup[group]
        values = {}
        for label, field in (("trf_tandem", "mean_trf_tandem_bp"), ("pan_te_non_tandem", "mean_pan_te_non_tandem_bp"), ("other", "mean_other_bp")):
            values[label] = int(ref[field]) - int(other[field])
        total = int(ref["mean_genome_bp"]) - int(other["mean_genome_bp"])
        result.append({
            "reference_group": reference,
            "comparator_group": group,
            "reference_mean_genome_bp": ref["mean_genome_bp"],
            "comparator_mean_genome_bp": other["mean_genome_bp"],
            "total_size_difference_bp": total,
            "trf_tandem_difference_bp": values["trf_tandem"],
            "pan_te_non_tandem_difference_bp": values["pan_te_non_tandem"],
            "other_difference_bp": values["other"],
            "components_sum_bp": sum(values.values()),
            "components_match_total": int(sum(values.values()) == total),
            "interpretation": "descriptive_group_mean_difference; phased_haplotypes_are_not_independent_replicates",
        })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mask-manifest", type=Path, required=True)
    parser.add_argument("--tandem-root", type=Path, required=True)
    parser.add_argument("--tandem-genome-summary", type=Path,
                        help="Optional non-redundant TRF summary used as an exact union QC at the same score threshold.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-trf-score", type=int, default=50)
    parser.add_argument("--reference-group", default="PZI")
    args = parser.parse_args()
    if args.min_trf_score < 0:
        raise SystemExit("ERROR: --min-trf-score must be non-negative")

    manifest = {r["target"]: r for r in read_tsv(args.mask_manifest)}
    missing = [target for target in TARGET_ORDER if target not in manifest]
    if missing:
        raise SystemExit("ERROR: mask manifest is missing " + ", ".join(missing))
    expected_trf = {}
    if args.tandem_genome_summary:
        expected_trf = {r["target"]: int(r["exclusive_union_bp"]) for r in read_tsv(args.tandem_genome_summary)}

    rows, qc = [], []
    for target in TARGET_ORDER:
        item = manifest[target]
        fasta = Path(item["mask_input_fasta"])
        rmout = Path(item["repeatmasker_out"])
        trf = args.tandem_root / target / f"{target}.trf.tsv"
        for path in (fasta, rmout, trf):
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"ERROR: missing required input for {target}: {path}")
        lengths = fasta_lengths(fasta)
        trf_by_seq = trf_union(trf, args.min_trf_score)
        rm_by_seq = repeatmasker_union(rmout)
        genome_bp = sum(lengths.values())
        tandem_bp = te_union_bp = overlap_bp = 0
        unknown_trf_sequences = sorted(set(trf_by_seq) - set(lengths))
        unknown_rm_sequences = sorted(set(rm_by_seq) - set(lengths))
        for seqid in lengths:
            tandem = trf_by_seq.get(seqid, [])
            te = rm_by_seq.get(seqid, [])
            tandem_bp += interval_length(tandem)
            te_union_bp += interval_length(te)
            overlap_bp += intersection_length(tandem, te)
        te_non_tandem_bp = te_union_bp - overlap_bp
        other_bp = genome_bp - tandem_bp - te_non_tandem_bp
        valid = other_bp >= 0 and tandem_bp + te_non_tandem_bp + other_bp == genome_bp
        summary_bp = expected_trf.get(target)
        summary_match = "" if summary_bp is None else int(summary_bp == tandem_bp)
        rows.append({
            "target": target, "group": item["group"], "genome_bp": genome_bp,
            "trf_tandem_bp": tandem_bp,
            "pan_te_non_tandem_bp": te_non_tandem_bp,
            "other_bp": other_bp,
            "repeatmasker_union_bp": te_union_bp,
            "repeatmasker_overlap_bp": overlap_bp,
            "trf_tandem_fraction": f"{tandem_bp / genome_bp:.8f}",
            "pan_te_non_tandem_fraction": f"{te_non_tandem_bp / genome_bp:.8f}",
            "other_fraction": f"{other_bp / genome_bp:.8f}",
            "repeatmasker_overlap_of_trf_fraction": f"{overlap_bp / tandem_bp:.8f}" if tandem_bp else "0",
            "partition_definition": "TRF_union | RepeatMasker_union_outside_TRF | other",
        })
        qc.append({
            "target": target, "genome_bp": genome_bp, "trf_tandem_bp": tandem_bp,
            "repeatmasker_union_bp": te_union_bp, "repeatmasker_overlap_bp": overlap_bp,
            "pan_te_non_tandem_bp": te_non_tandem_bp, "other_bp": other_bp,
            "partition_sum_bp": tandem_bp + te_non_tandem_bp + other_bp,
            "partition_matches_genome": int(valid),
            "tandem_summary_bp": "" if summary_bp is None else summary_bp,
            "tandem_summary_matches_exact_union": summary_match,
            "trf_sequence_ids_absent_from_fasta": len(unknown_trf_sequences),
            "repeatmasker_sequence_ids_absent_from_fasta": len(unknown_rm_sequences),
        })

    groups = group_summary(rows)
    contrasts = contrast_summary(groups, args.reference_group)
    fields = [
        "target", "group", "genome_bp", "trf_tandem_bp", "pan_te_non_tandem_bp", "other_bp",
        "repeatmasker_union_bp", "repeatmasker_overlap_bp", "trf_tandem_fraction",
        "pan_te_non_tandem_fraction", "other_fraction", "repeatmasker_overlap_of_trf_fraction",
        "partition_definition",
    ]
    write_tsv(args.output_dir / "tandem_main_partition.tsv", fields, rows)
    write_tsv(args.output_dir / "tandem_main_group_partition.tsv", list(groups[0]) if groups else [], groups)
    write_tsv(args.output_dir / "tandem_main_excess_contrasts.tsv", list(contrasts[0]) if contrasts else [], contrasts)
    write_tsv(args.output_dir / "tandem_main_partition_qc.tsv", list(qc[0]) if qc else [], qc)
    metadata = {
        "definition": "Every base is assigned exactly once to TRF tandem union, pan-TE RepeatMasker union outside TRF, or other sequence.",
        "trf_min_score": args.min_trf_score,
        "reference_group": args.reference_group,
        "warning": "TRF--RepeatMasker overlap is retained inside the TRF category to prevent double counting; this is a sequence-architecture partition, not a TE-origin assignment.",
    }
    with open(args.output_dir / "tandem_main_partition_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    print(json.dumps({"targets": len(rows), "all_partition_qc_pass": all(int(r["partition_matches_genome"]) for r in qc), "output_dir": str(args.output_dir)}))


if __name__ == "__main__":
    main()
