#!/usr/bin/env python3
"""Prepare and classify orthologous insertion sites for curated pan-TE components."""

import argparse
import bisect
import collections
import csv
import hashlib
import math
import random
from pathlib import Path


INCLUDE_DECISIONS = {"include", "keep", "yes", "pass"}
EXCLUDE_DECISIONS = {"exclude", "drop", "no", "fail"}


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def selected_components(path, max_components):
    rows = read_tsv(path)
    selected = []
    for row in sorted(rows, key=lambda value: int(value["positive_rank"])):
        decision = row.get("manual_decision", "").strip().lower()
        if decision in EXCLUDE_DECISIONS:
            continue
        include = decision in INCLUDE_DECISIONS or (
            not decision and row.get("automated_status") == "pass_preliminary"
        )
        if include:
            row = dict(row)
            row["curated_family"] = row.get("curated_family") or row["pan_family"]
            selected.append(row)
        if len(selected) >= max_components:
            break
    if not selected:
        raise SystemExit(
            "ERROR: no components are eligible; review component_decisions.tsv and set manual_decision=include"
        )
    return selected


def read_fai(path):
    lengths = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def parse_rmout_candidates(path, selected, full_hit_fraction):
    selected_by_family = {row["pan_family"]: row for row in selected}
    result = collections.defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            family = fields[9]
            row = selected_by_family.get(family)
            if row is None or row["class_family"].startswith("LTR/"):
                continue
            try:
                start = int(fields[5]) - 1
                end = int(fields[6])
                divergence = float(fields[1])
                score = int(fields[0])
                library_length = int(row["library_sequence_length_bp"])
            except ValueError:
                continue
            if end <= start or end - start < library_length * full_hit_fraction:
                continue
            result[family].append({
                "seqid": fields[4], "start0": start, "end0": end,
                "strand": fields[8], "source_type": "repeatmasker_component",
                "score": score, "divergence": divergence,
            })
    return result


def parse_intact_candidates(path, query_target, selected, min_overlap_fraction):
    selected_by_family = {row["pan_family"]: row for row in selected}
    result = collections.defaultdict(list)
    for row in read_tsv(path):
        family = row["pan_family"]
        selected_row = selected_by_family.get(family)
        if (
            row["target"] != query_target or selected_row is None
            or not selected_row["class_family"].startswith("LTR/")
            or row["assignment_method"] != "maximum_coordinate_overlap"
            or float(row["overlap_fraction"]) < min_overlap_fraction
        ):
            continue
        result[family].append({
            "seqid": row["seqid"],
            "start0": int(row["start"]) - 1,
            "end0": int(row["end"]),
            "strand": row["strand"],
            "source_type": "intact_ltr_rt",
            "score": int(row["overlap_bp"]),
            "divergence": float(row["ltr_divergence"]) * 100,
        })
    return result


def stable_seed(seed, family):
    digest = hashlib.sha256(f"{seed}:{family}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def choose_loci(candidates, contig_lengths, flank_bp, max_loci, seed):
    eligible = [
        row for row in candidates
        if row["seqid"] in contig_lengths
        and row["start0"] >= flank_bp
        and row["end0"] + flank_bp <= contig_lengths[row["seqid"]]
    ]
    rng = random.Random(seed)
    rng.shuffle(eligible)
    chosen = []
    occupied = collections.defaultdict(list)
    for row in eligible:
        interval = (row["start0"] - flank_bp, row["end0"] + flank_bp)
        if any(not (interval[1] <= left or interval[0] >= right) for left, right in occupied[row["seqid"]]):
            continue
        chosen.append(row)
        occupied[row["seqid"]].append(interval)
        if len(chosen) >= max_loci:
            break
    return chosen


def command_prepare(args):
    selected = selected_components(args.decisions, args.max_components)
    contig_lengths = read_fai(args.fai)
    rm_candidates = parse_rmout_candidates(args.rmout, selected, args.full_hit_fraction)
    intact_candidates = parse_intact_candidates(
        args.intact_ltr_table, args.query_target, selected, args.min_intact_overlap_fraction
    )
    loci = []
    beds = []
    for component in selected:
        family = component["pan_family"]
        candidates = intact_candidates[family] if component["class_family"].startswith("LTR/") else rm_candidates[family]
        chosen = choose_loci(
            candidates, contig_lengths, args.flank_bp, args.max_loci_per_component,
            stable_seed(args.seed, family),
        )
        for index, candidate in enumerate(sorted(chosen, key=lambda row: (row["seqid"], row["start0"])), 1):
            site_id = f"{args.query_target}__{family}__{index:05d}"
            loci.append({
                "site_id": site_id,
                "query_target": args.query_target,
                "pan_family": family,
                "curated_family": component["curated_family"],
                "class_family": component["class_family"],
                "positive_rank": component["positive_rank"],
                "seqid": candidate["seqid"],
                "start0": candidate["start0"],
                "end0": candidate["end0"],
                "length_bp": candidate["end0"] - candidate["start0"],
                "strand": candidate["strand"],
                "source_type": candidate["source_type"],
                "source_score": candidate["score"],
                "source_divergence_pct": f"{candidate['divergence']:.6f}",
                "flank_bp": args.flank_bp,
            })
            beds.append((candidate["seqid"], candidate["start0"] - args.flank_bp, candidate["start0"], site_id + "|L"))
            beds.append((candidate["seqid"], candidate["end0"], candidate["end0"] + args.flank_bp, site_id + "|R"))
    if not loci:
        raise SystemExit("ERROR: no loci passed boundary and completeness filters")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "query_loci.tsv", list(loci[0]), loci)
    with open(args.output_dir / "query_flanks.bed", "w", encoding="utf-8") as handle:
        for row in beds:
            handle.write("\t".join(map(str, row)) + "\n")
    print(f"Prepared {len(loci)} loci from {len(set(row['pan_family'] for row in loci))} components")


def parse_paf(path, min_mapq, min_qcov, min_identity, uniqueness_ratio):
    by_query = collections.defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            try:
                qlen = int(fields[1])
                qstart = int(fields[2])
                qend = int(fields[3])
                nmatch = int(fields[9])
                alen = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue
            qcov = (qend - qstart) / qlen if qlen else 0
            identity = nmatch / alen if alen else 0
            if mapq < min_mapq or qcov < min_qcov or identity < min_identity:
                continue
            by_query[fields[0]].append({
                "strand": fields[4], "target_seqid": fields[5],
                "target_start0": int(fields[7]), "target_end0": int(fields[8]),
                "nmatch": nmatch, "alignment_bp": alen, "mapq": mapq,
                "query_coverage": qcov, "identity": identity,
            })
    result = {}
    for query, hits in by_query.items():
        hits.sort(key=lambda row: (row["nmatch"], row["mapq"], row["identity"]), reverse=True)
        best = hits[0]
        best["unique"] = not (
            len(hits) > 1 and hits[1]["nmatch"] >= best["nmatch"] * uniqueness_ratio
        )
        result[query] = best
    return result


def parse_target_repeats(path, component_to_curated):
    intervals = collections.defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            try:
                start = int(fields[5]) - 1
                end = int(fields[6])
            except ValueError:
                continue
            family = fields[9]
            intervals[fields[4]].append({
                "start": start, "end": end, "pan_family": family,
                "curated_family": component_to_curated.get(family, family),
            })
    index = {}
    for seqid, values in intervals.items():
        values.sort(key=lambda row: row["start"])
        starts = [row["start"] for row in values]
        prefix = []
        maximum = 0
        for row in values:
            maximum = max(maximum, row["end"])
            prefix.append(maximum)
        index[seqid] = (values, starts, prefix)
    return index


def overlapping_repeats(seq_index, start, end):
    if not seq_index:
        return []
    values, starts, prefix = seq_index
    cursor = bisect.bisect_left(starts, end) - 1
    result = []
    while cursor >= 0:
        if prefix[cursor] <= start:
            break
        row = values[cursor]
        overlap = min(end, row["end"]) - max(start, row["start"])
        if overlap > 0:
            result.append((row, overlap))
        cursor -= 1
    return result


def classify_site(locus, left, right, repeat_index, empty_max_gap, max_gap, min_te_overlap):
    base = {
        "left_mapq": left.get("mapq", "") if left else "",
        "right_mapq": right.get("mapq", "") if right else "",
        "left_identity": f"{left['identity']:.6f}" if left else "",
        "right_identity": f"{right['identity']:.6f}" if right else "",
        "target_seqid": "", "target_interval_start0": "", "target_interval_end0": "",
        "target_gap_bp": "", "matched_te_overlap_bp": 0,
    }
    if not left or not right:
        return "unmapped_flank", base
    if not left["unique"] or not right["unique"]:
        return "nonunique_flank", base
    if left["target_seqid"] != right["target_seqid"] or left["strand"] != right["strand"]:
        return "discordant_flanks", base
    if left["strand"] == "+":
        start = left["target_end0"]
        end = right["target_start0"]
    else:
        start = right["target_end0"]
        end = left["target_start0"]
    gap = end - start
    base.update({
        "target_seqid": left["target_seqid"],
        "target_interval_start0": start,
        "target_interval_end0": end,
        "target_gap_bp": gap,
    })
    if gap < -empty_max_gap:
        return "overlapping_or_rearranged_flanks", base
    if gap > max_gap:
        return "large_structural_gap", base
    overlaps = overlapping_repeats(repeat_index.get(left["target_seqid"]), start, end)
    matched = sum(
        overlap for row, overlap in overlaps
        if row["curated_family"] == locus["curated_family"]
    )
    base["matched_te_overlap_bp"] = matched
    if matched >= min_te_overlap:
        return "occupied", base
    if gap <= empty_max_gap:
        return "empty", base
    if overlaps:
        return "other_repeat_in_gap", base
    return "unresolved_gap", base


def command_summarize(args):
    loci = read_tsv(args.loci)
    decisions = read_tsv(args.decisions)
    component_to_curated = {
        row["pan_family"]: (row.get("curated_family") or row["pan_family"])
        for row in decisions
    }
    targets = read_tsv(args.target_manifest)
    classifications = []
    for target_row in targets:
        mappings = parse_paf(
            Path(target_row["paf"]), args.min_mapq, args.min_query_coverage,
            args.min_identity, args.uniqueness_ratio,
        )
        repeat_index = parse_target_repeats(Path(target_row["repeatmasker_out"]), component_to_curated)
        for locus in loci:
            left = mappings.get(locus["site_id"] + "|L")
            right = mappings.get(locus["site_id"] + "|R")
            state, details = classify_site(
                locus, left, right, repeat_index, args.empty_max_gap,
                args.max_occupied_gap, args.min_te_overlap_bp,
            )
            classifications.append({
                **locus,
                "comparison_target": target_row["target"],
                "comparison_role": target_row["role"],
                "site_state": state,
                **details,
            })

    by_site = collections.defaultdict(list)
    for row in classifications:
        by_site[row["site_id"]].append(row)
    evidence_rows = []
    for locus in loci:
        values = by_site[locus["site_id"]]
        within = [row["site_state"] for row in values if row["comparison_role"] == "within_pzi"]
        psp = [row["site_state"] for row in values if row["comparison_role"] == "psp"]
        outgroups = [row["site_state"] for row in values if row["comparison_role"] == "outgroup"]
        within_occupied = any(value == "occupied" for value in within)
        if within_occupied and any(value == "occupied" for value in psp):
            evidence = "shared_pzi_psp"
        elif within_occupied and psp and all(value == "empty" for value in psp) and any(value == "empty" for value in outgroups):
            evidence = "pzi_specific_supported"
        elif within_occupied and psp and all(value == "empty" for value in psp):
            evidence = "pzi_specific_no_outgroup"
        elif any(value == "empty" for value in within):
            evidence = "query_haplotype_specific_or_assembly"
        else:
            evidence = "unresolved"
        evidence_rows.append({
            **locus,
            "evidence_class": evidence,
            "within_pzi_states": ",".join(within),
            "psp_states": ",".join(psp),
            "outgroup_states": ",".join(outgroups),
        })

    counter = collections.defaultdict(collections.Counter)
    locus_count = collections.Counter()
    metadata = {}
    for row in evidence_rows:
        key = (row["query_target"], row["pan_family"], row["curated_family"], row["class_family"])
        counter[key][row["evidence_class"]] += 1
        locus_count[key] += 1
        metadata[key] = row
    summary = []
    evidence_order = [
        "pzi_specific_supported", "pzi_specific_no_outgroup", "shared_pzi_psp",
        "query_haplotype_specific_or_assembly", "unresolved",
    ]
    for key in sorted(counter, key=lambda value: (int(metadata[value]["positive_rank"]), value[0])):
        counts = counter[key]
        total = locus_count[key]
        resolved = counts["pzi_specific_supported"] + counts["pzi_specific_no_outgroup"] + counts["shared_pzi_psp"]
        summary.append({
            "query_target": key[0], "pan_family": key[1], "curated_family": key[2],
            "class_family": key[3], "positive_rank": metadata[key]["positive_rank"],
            "tested_loci": total, **{name: counts[name] for name in evidence_order},
            "resolved_loci": resolved,
            "pzi_specific_fraction_of_resolved": (
                f"{(counts['pzi_specific_supported'] + counts['pzi_specific_no_outgroup']) / resolved:.8f}"
                if resolved else ""
            ),
        })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "orthosite_classification.tsv", list(classifications[0]), classifications)
    write_tsv(args.output_dir / "orthosite_site_evidence.tsv", list(evidence_rows[0]), evidence_rows)
    write_tsv(args.output_dir / "orthosite_family_summary.tsv", list(summary[0]), summary)
    print(f"Classified {len(evidence_rows)} query loci against {len(targets)} genomes")


def command_combine(args):
    evidence = []
    summaries = []
    for source_dir in args.source_dir:
        evidence.extend(read_tsv(source_dir / "orthosite_site_evidence.tsv"))
        summaries.extend(read_tsv(source_dir / "orthosite_family_summary.tsv"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "orthosite_site_evidence.tsv", list(evidence[0]), evidence)
    write_tsv(args.output_dir / "orthosite_family_summary.tsv", list(summaries[0]), summaries)
    print(f"Combined {len(evidence)} loci and {len(summaries)} component summaries")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--query-target", required=True)
    prepare.add_argument("--rmout", required=True, type=Path)
    prepare.add_argument("--fai", required=True, type=Path)
    prepare.add_argument("--decisions", required=True, type=Path)
    prepare.add_argument("--intact-ltr-table", required=True, type=Path)
    prepare.add_argument("--output-dir", required=True, type=Path)
    prepare.add_argument("--max-components", type=int, default=10)
    prepare.add_argument("--max-loci-per-component", type=int, default=200)
    prepare.add_argument("--flank-bp", type=int, default=1500)
    prepare.add_argument("--full-hit-fraction", type=float, default=0.80)
    prepare.add_argument("--min-intact-overlap-fraction", type=float, default=0.50)
    prepare.add_argument("--seed", type=int, default=20260815)
    prepare.set_defaults(func=command_prepare)

    summarize = subparsers.add_parser("summarize")
    summarize.add_argument("--loci", required=True, type=Path)
    summarize.add_argument("--target-manifest", required=True, type=Path)
    summarize.add_argument("--decisions", required=True, type=Path)
    summarize.add_argument("--output-dir", required=True, type=Path)
    summarize.add_argument("--min-mapq", type=int, default=30)
    summarize.add_argument("--min-query-coverage", type=float, default=0.80)
    summarize.add_argument("--min-identity", type=float, default=0.80)
    summarize.add_argument("--uniqueness-ratio", type=float, default=0.95)
    summarize.add_argument("--empty-max-gap", type=int, default=100)
    summarize.add_argument("--max-occupied-gap", type=int, default=50000)
    summarize.add_argument("--min-te-overlap-bp", type=int, default=50)
    summarize.set_defaults(func=command_summarize)

    combine = subparsers.add_parser("combine")
    combine.add_argument("--source-dir", required=True, action="append", type=Path)
    combine.add_argument("--output-dir", required=True, type=Path)
    combine.set_defaults(func=command_combine)
    return parser


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
