#!/usr/bin/env python3
"""Summarize soloLTRseeker outputs without treating haplotypes as replicates."""

import argparse
import collections
import csv
import math
import statistics
from pathlib import Path


TARGET_ORDER = [
    "pziA1", "pziB1",
    "pspA1", "pspB1",
    "ppaA", "ppaB", "pcam", "ppe", "gtr", "rch",
]
TARGET_GROUP = {
    "pziA1": "PZI", "pziB1": "PZI",
    "pspA1": "PSP", "pspB1": "PSP",
    "ppaA": "PPA", "ppaB": "PPA", "pcam": "PCAM", "ppe": "PPE",
    "gtr": "GTR", "rch": "RCH",
}


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_attributes(raw):
    result = {}
    for item in raw.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key.strip().lower()] = value.strip()
    return result


def normalize_superfamily(value):
    value = value.replace("LTR/", "").replace("_LTR_retrotransposon", "")
    value = value.replace("LTR_retrotransposon", "unknown").strip("/_")
    lower = value.lower()
    if "gypsy" in lower:
        return "Gypsy"
    if "copia" in lower:
        return "Copia"
    return "unknown"


def normalized_lineage(superfamily, lineage):
    if not lineage or lineage.upper() == "NA" or lineage == "unknown":
        return superfamily
    return f"{superfamily}/{lineage}"


def parse_ltr_lengths(path):
    counts = collections.Counter()
    raw = 0
    filtered = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 9:
                continue
            raw += 1
            query_id = fields[0]
            if fields[6] != "PASS" or fields[7] != "PASS":
                continue
            filtered += 1
            annotation = query_id.split("#", 1)[1] if "#" in query_id else "unknown/NA"
            classification, lineage = (annotation.split("/", 1) + ["NA"])[:2]
            superfamily = normalize_superfamily(classification)
            counts[normalized_lineage(superfamily, lineage)] += 1
    return counts, raw, filtered


def parse_solo_gff(path, min_length=0):
    counts = collections.Counter()
    total = 0
    zero_mismatch_tsd = 0
    removed = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            length = int(fields[4]) - int(fields[3]) + 1
            if length < min_length:
                removed += 1
                continue
            attrs = parse_attributes(fields[8])
            superfamily = normalize_superfamily(attrs.get("classification", "unknown"))
            lineage = normalized_lineage(superfamily, attrs.get("lineage", "NA"))
            counts[lineage] += 1
            total += 1
            tsd = attrs.get("tsd", "")
            if tsd.endswith("_0"):
                zero_mismatch_tsd += 1
    return counts, total, zero_mismatch_tsd, removed


def parse_hit_count(path):
    values = {}
    if not path.is_file():
        return values
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 4:
                try:
                    values[fields[-2]] = int(fields[-1])
                except ValueError:
                    pass
    return values


def find_input(root, target, stable_name, suffix):
    stable = root / target / stable_name
    if stable.is_file():
        return stable
    matches = sorted((root / target).glob(suffix))
    if len(matches) != 1:
        raise FileNotFoundError(f"{target}: expected one {suffix}, found {len(matches)}")
    return matches[0]


def mean(values):
    return sum(values) / len(values) if values else math.nan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sololtr-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-solo-length", type=int, default=0,
                        help="Minimum soloLTR length in bp; shorter calls are removed before counting (default: 0 = keep all)")
    parser.add_argument("--allow-missing-targets", action="store_true")
    args = parser.parse_args()

    rows = []
    qc_rows = []
    missing = []
    for target in TARGET_ORDER:
        try:
            lengths = find_input(args.sololtr_root, target, f"{target}.LTRRT_lengths.txt", "*_LTRRT_lengths.txt")
            solo = find_input(args.sololtr_root, target, f"{target}.soloLTR.gff3", "*_soloLTR.gff3")
        except FileNotFoundError:
            missing.append(target)
            continue
        hit_count_candidates = sorted((args.sololtr_root / target).rglob("*_hit_count.txt"))
        hit_count = parse_hit_count(hit_count_candidates[0]) if hit_count_candidates else {}
        intact_counts, intact_raw, intact_filtered = parse_ltr_lengths(lengths)
        solo_counts, solo_total, zero_mismatch, solo_removed = parse_solo_gff(solo, args.min_solo_length)
        for lineage in sorted(set(intact_counts) | set(solo_counts)):
            intact = intact_counts[lineage]
            solo_count = solo_counts[lineage]
            rows.append({
                "target": target,
                "group": TARGET_GROUP[target],
                "lineage": lineage,
                "superfamily": lineage.split("/", 1)[0],
                "intact_count": intact,
                "solo_count": solo_count,
                "solo_to_intact_ratio": f"{solo_count / intact:.8f}" if intact else "",
                "deletion_index": f"{solo_count / (solo_count + intact):.8f}" if solo_count + intact else "",
            })
        qc_rows.append({
            "target": target,
            "group": TARGET_GROUP[target],
            "raw_intact_records": intact_raw,
            "quality_filtered_intact_records": intact_filtered,
            "final_solo_ltr_records": solo_total,
            "zero_mismatch_tsd_solo_records": zero_mismatch,
            "zero_mismatch_tsd_fraction": f"{zero_mismatch / solo_total:.8f}" if solo_total else "",
            "min_solo_length_bp": args.min_solo_length,
            "solo_records_removed_by_length_filter": solo_removed,
            "pipeline_raw_fl_ltrrt": hit_count.get("raw_fl_LTRRT", ""),
            "pipeline_pass_fl_ltrrt": hit_count.get("PASS_fl_LTRRT", ""),
            "pipeline_tg_ca": hit_count.get("TG..CA", ""),
            "pipeline_cd_hit_queries": hit_count.get("cd-hit", ""),
            "pipeline_final_tsd": hit_count.get("TSD", ""),
        })
    if missing and not args.allow_missing_targets:
        raise SystemExit("ERROR: missing soloLTRseeker outputs for: " + ", ".join(missing))

    by_group_lineage = collections.defaultdict(list)
    for row in rows:
        by_group_lineage[(row["group"], row["lineage"])].append(row)
    group_rows = []
    for (group, lineage), values in sorted(by_group_lineage.items()):
        intact = [int(row["intact_count"]) for row in values]
        solo = [int(row["solo_count"]) for row in values]
        ratios = [float(row["solo_to_intact_ratio"]) for row in values if row["solo_to_intact_ratio"]]
        deletion = [float(row["deletion_index"]) for row in values if row["deletion_index"]]
        group_rows.append({
            "group": group,
            "lineage": lineage,
            "superfamily": lineage.split("/", 1)[0],
            "genome_count": len(values),
            "mean_intact_count": f"{mean(intact):.8f}",
            "min_intact_count": min(intact),
            "max_intact_count": max(intact),
            "mean_solo_count": f"{mean(solo):.8f}",
            "min_solo_count": min(solo),
            "max_solo_count": max(solo),
            "mean_solo_to_intact_ratio": f"{mean(ratios):.8f}" if ratios else "",
            "min_solo_to_intact_ratio": f"{min(ratios):.8f}" if ratios else "",
            "max_solo_to_intact_ratio": f"{max(ratios):.8f}" if ratios else "",
            "mean_deletion_index": f"{mean(deletion):.8f}" if deletion else "",
        })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.output_dir / "sololtr_by_target_lineage.tsv",
        ["target", "group", "lineage", "superfamily", "intact_count", "solo_count", "solo_to_intact_ratio", "deletion_index"],
        rows,
    )
    write_tsv(
        args.output_dir / "sololtr_group_summary.tsv",
        [
            "group", "lineage", "superfamily", "genome_count", "mean_intact_count", "min_intact_count",
            "max_intact_count", "mean_solo_count", "min_solo_count", "max_solo_count",
            "mean_solo_to_intact_ratio", "min_solo_to_intact_ratio", "max_solo_to_intact_ratio",
            "mean_deletion_index",
        ],
        group_rows,
    )
    write_tsv(
        args.output_dir / "sololtr_qc.tsv",
        list(qc_rows[0]) if qc_rows else [],
        qc_rows,
    )
    print(f"Summarized {len(rows)} target-lineage combinations")


if __name__ == "__main__":
    main()
