#!/usr/bin/env python3
"""Post-process standardized Prunus pan-TE results for figures and audits."""

import argparse
import bisect
import collections
import csv
import heapq
import json
import math
import statistics
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


TARGET_ORDER = [
    "pziA1", "pziA2", "pziB1", "pziB2",
    "pspA1", "pspA2", "pspB1", "pspB2",
    "ppaA", "ppaB", "pcam", "ppe",
]

TARGET_META = {
    "pziA1": ("PZI", "A", "PZI-A"),
    "pziA2": ("PZI", "A", "PZI-A"),
    "pziB1": ("PZI", "B", "PZI-B"),
    "pziB2": ("PZI", "B", "PZI-B"),
    "pspA1": ("PSP", "A", "PSP-A"),
    "pspA2": ("PSP", "A", "PSP-A"),
    "pspB1": ("PSP", "B", "PSP-B"),
    "pspB2": ("PSP", "B", "PSP-B"),
    "ppaA": ("PPA", "A", "PPA-A"),
    "ppaB": ("PPA", "B", "PPA-B"),
    "pcam": ("PCAM", "", "PCAM"),
    "ppe": ("PPE", "", "PPE"),
}

COMPOSITION_ORDER = [
    "LTR/Gypsy", "LTR/Copia", "Other LTR", "TIR/DNA",
    "Helitron", "LINE", "SINE", "Other TE", "Non-TE",
]


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def mean(values):
    return sum(values) / len(values) if values else 0.0


def classify_composition(class_family):
    value = class_family.strip()
    lower = value.lower()
    if lower.startswith("ltr/gypsy"):
        return "LTR/Gypsy"
    if lower.startswith("ltr/copia"):
        return "LTR/Copia"
    if lower.startswith("ltr"):
        return "Other LTR"
    if "helitron" in lower or lower.startswith("rc/"):
        return "Helitron"
    if lower.startswith("dna") or lower.startswith("tir") or lower.startswith("mite"):
        return "TIR/DNA"
    if lower.startswith("line"):
        return "LINE"
    if lower.startswith("sine"):
        return "SINE"
    return "Other TE"


def validate_targets(rows, label):
    observed = {row["target"] for row in rows}
    missing = [target for target in TARGET_ORDER if target not in observed]
    if missing:
        raise ValueError(f"{label} is missing targets: {', '.join(missing)}")


def build_composition(total_rows, exclusive_results):
    validate_targets(total_rows, "total summary")
    totals = {row["target"]: row for row in total_rows}
    result_by_target = {row["target"]: row for row in exclusive_results}
    missing = [target for target in TARGET_ORDER if target not in result_by_target]
    if missing:
        raise ValueError(f"exclusive RepeatMasker results are missing targets: {', '.join(missing)}")

    output = []
    qc = []
    for target in TARGET_ORDER:
        genome_bp = int(totals[target]["genome_bp"])
        repeat_bp = int(totals[target]["repeat_union_bp"])
        result = result_by_target[target]
        exclusive_repeat_bp = result["repeat_union_bp"]
        if exclusive_repeat_bp != repeat_bp:
            raise ValueError(
                f"{target}: exclusive RepeatMasker union ({exclusive_repeat_bp}) does not match "
                f"Step 15 repeat union ({repeat_bp})"
            )
        adjusted = collections.Counter(result["category_bp"])
        adjusted["Non-TE"] = max(0, genome_bp - repeat_bp)
        for category in COMPOSITION_ORDER:
            bp = adjusted[category]
            output.append({
                "target": target,
                "group": TARGET_META[target][0],
                "subgenome": TARGET_META[target][1],
                "display_group": TARGET_META[target][2],
                "category": category,
                "bp": bp,
                "genome_fraction": f"{bp / genome_bp:.8f}" if genome_bp else "0",
                "genome_bp": genome_bp,
                "repeat_union_bp": repeat_bp,
            })
        qc.append({
            "target": target,
            "genome_bp": genome_bp,
            "repeat_union_bp": repeat_bp,
            "raw_hit_sum_bp": result["raw_hit_sum_bp"],
            "raw_hit_sum_to_repeat_union_ratio": (
                f"{result['raw_hit_sum_bp'] / repeat_bp:.8f}" if repeat_bp else "0"
            ),
            "exclusive_assignment_sum_bp": sum(result["category_bp"].values()),
            "exclusive_sum_to_repeat_union_ratio": (
                f"{sum(result['category_bp'].values()) / repeat_bp:.8f}" if repeat_bp else "0"
            ),
        })
    return output, qc


def build_driver_table(family_rows, cumulative_threshold, library_audit):
    values = collections.defaultdict(lambda: collections.defaultdict(int))
    class_by_family = {}
    for row in family_rows:
        key = (row["class_family"], row["pan_family"])
        values[key][row["target"]] = int(row["exclusive_bp"])
        class_by_family[row["pan_family"]] = row["class_family"]

    records = []
    for (class_family, pan_family), target_values in values.items():
        pzi_a = mean([target_values["pziA1"], target_values["pziA2"]])
        pzi_b = mean([target_values["pziB1"], target_values["pziB2"]])
        psp_a = mean([target_values["pspA1"], target_values["pspA2"]])
        psp_b = mean([target_values["pspB1"], target_values["pspB2"]])
        pzi_mean = mean([target_values[target] for target in TARGET_ORDER[:4]])
        psp_mean = mean([target_values[target] for target in TARGET_ORDER[4:8]])
        context = [target_values[target] for target in ("ppaA", "ppaB", "pcam", "ppe")]
        delta_a = pzi_a - psp_a
        delta_b = pzi_b - psp_b
        delta_mean = pzi_mean - psp_mean
        positive_excess = max(0.0, delta_mean)
        records.append({
            "class_family": class_family,
            "pan_family": pan_family,
            "pzi_a_mean_bp": round(pzi_a),
            "psp_a_mean_bp": round(psp_a),
            "delta_a_bp": round(delta_a),
            "pzi_b_mean_bp": round(pzi_b),
            "psp_b_mean_bp": round(psp_b),
            "delta_b_bp": round(delta_b),
            "pzi_mean_bp": round(pzi_mean),
            "psp_mean_bp": round(psp_mean),
            "delta_pzi_psp_mean_bp": round(delta_mean),
            "positive_excess_bp": round(positive_excess),
            "ppa_pcam_ppe_median_bp": round(statistics.median(context)),
            "expanded_in_both_subgenomes": int(delta_a > 0 and delta_b > 0),
            "consistent_across_pzi_haplotypes": int(
                target_values["pziA1"] > psp_a
                and target_values["pziA2"] > psp_a
                and target_values["pziB1"] > psp_b
                and target_values["pziB2"] > psp_b
            ),
            "library_sequence_length_bp": library_audit.get(pan_family, {}).get("sequence_length_bp", ""),
            "library_gc_fraction": library_audit.get(pan_family, {}).get("gc_fraction", ""),
            "library_base_entropy_bits": library_audit.get(pan_family, {}).get("base_entropy_bits", ""),
            "library_low_complexity_flag": library_audit.get(pan_family, {}).get("low_complexity_flag", ""),
            **{f"{target}_bp": target_values[target] for target in TARGET_ORDER},
        })

    records.sort(key=lambda row: (row["positive_excess_bp"], row["pan_family"]), reverse=True)
    total_positive = sum(row["positive_excess_bp"] for row in records)
    cumulative = 0.0
    rank = 0
    for row in records:
        if row["positive_excess_bp"] > 0:
            rank += 1
            previous_fraction = cumulative / total_positive if total_positive else 0.0
            cumulative += row["positive_excess_bp"]
            row["positive_rank"] = rank
            row["cumulative_positive_excess_fraction"] = f"{cumulative / total_positive:.8f}"
            row["driver_selected"] = int(previous_fraction < cumulative_threshold)
        else:
            row["positive_rank"] = ""
            row["cumulative_positive_excess_fraction"] = ""
            row["driver_selected"] = 0
    return records


def build_mapping_qc(mapping_dir):
    rows = []
    for target in TARGET_ORDER:
        path = mapping_dir / f"{target}.original_to_pan_family.tsv"
        if not path.is_file():
            raise FileNotFoundError(path)
        records = read_tsv(path)
        mapped = [record for record in records if record["mapping_status"] == "mapped"]
        identities = [float(record["identity"]) for record in mapped]
        rows.append({
            "target": target,
            "group": TARGET_META[target][0],
            "original_family_count": len(records),
            "mapped_family_count": len(mapped),
            "unmapped_family_count": len(records) - len(mapped),
            "mapped_fraction": f"{len(mapped) / len(records):.8f}" if records else "0",
            "median_mapped_identity": f"{statistics.median(identities):.8f}" if identities else "",
            "mapped_identity_ge_0.80_fraction": (
                f"{sum(value >= 0.80 for value in identities) / len(identities):.8f}" if identities else ""
            ),
        })
    return rows


def parse_rmout_hits(path):
    hits_by_seq = collections.defaultdict(list)
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            try:
                score = int(fields[0])
                divergence = float(fields[1])
                start = int(fields[5]) - 1
                end = int(fields[6])
            except ValueError:
                continue
            if end <= start or not math.isfinite(divergence):
                continue
            hits_by_seq[fields[4]].append({
                "start": start,
                "end": end,
                "score": score,
                "divergence": divergence,
                "family": fields[9],
                "class_family": fields[10],
            })
    return hits_by_seq


def assign_exclusive_segments(hits, bin_width):
    """Assign every covered base to exactly one RepeatMasker hit."""
    events = []
    raw = collections.Counter()
    for index, hit in enumerate(hits):
        span = hit["end"] - hit["start"]
        hit["span"] = span
        events.append((hit["start"], index, 1))
        events.append((hit["end"], index, -1))
        key = (hit["class_family"], hit["family"])
        raw[(key, "hit_count")] += 1
        raw[(key, "raw_hit_bp")] += span

    events.sort(key=lambda value: value[0])
    active = [False] * len(hits)
    heap = []
    family_bp = collections.Counter()
    category_bp = collections.Counter()
    divergence_bp = collections.Counter()
    repeat_union_bp = 0
    previous = None
    cursor = 0

    while cursor < len(events):
        position = events[cursor][0]
        while heap and not active[heap[0][-1]]:
            heapq.heappop(heap)
        if previous is not None and position > previous and heap:
            winner = hits[heap[0][-1]]
            length = position - previous
            family_key = (winner["class_family"], winner["family"])
            family_bp[family_key] += length
            category_bp[classify_composition(winner["class_family"])] += length
            bin_start = math.floor(winner["divergence"] / bin_width) * bin_width
            divergence_bp[(winner["class_family"], winner["family"], bin_start)] += length
            repeat_union_bp += length

        batch = []
        while cursor < len(events) and events[cursor][0] == position:
            batch.append(events[cursor])
            cursor += 1
        for _, index, action in batch:
            if action == -1:
                active[index] = False
        for _, index, action in batch:
            if action != 1:
                continue
            active[index] = True
            hit = hits[index]
            score_density = hit["score"] / hit["span"]
            priority = (
                -score_density,
                hit["divergence"],
                -hit["score"],
                -hit["span"],
                hit["family"],
                hit["class_family"],
                index,
            )
            heapq.heappush(heap, priority)
        previous = position

    return family_bp, category_bp, divergence_bp, repeat_union_bp, raw


def parse_rmout_one(task):
    target, group, path, bin_width, substitution_rate = task
    hits_by_seq = parse_rmout_hits(path)
    family_bp = collections.Counter()
    category_bp = collections.Counter()
    divergence_bp = collections.Counter()
    raw = collections.Counter()
    repeat_union_bp = 0
    for hits in hits_by_seq.values():
        seq_family, seq_category, seq_divergence, seq_union, seq_raw = assign_exclusive_segments(
            hits, bin_width
        )
        family_bp.update(seq_family)
        category_bp.update(seq_category)
        divergence_bp.update(seq_divergence)
        raw.update(seq_raw)
        repeat_union_bp += seq_union

    divergence_rows = []
    raw_divergence_hits = collections.Counter()
    for hits in hits_by_seq.values():
        for hit in hits:
            bin_start = math.floor(hit["divergence"] / bin_width) * bin_width
            raw_divergence_hits[(hit["class_family"], hit["family"], bin_start)] += 1
    for class_family, family, bin_start in sorted(divergence_bp):
        midpoint = bin_start + bin_width / 2
        age_mya = ""
        if substitution_rate:
            age_mya = midpoint / 100 / (2 * substitution_rate) / 1_000_000
        divergence_rows.append({
            "target": target,
            "group": group,
            "class_family": class_family,
            "pan_family": family,
            "divergence_bin_start_pct": f"{bin_start:.6f}",
            "divergence_bin_mid_pct": f"{midpoint:.6f}",
            "hit_count": raw_divergence_hits[(class_family, family, bin_start)],
            "exclusive_bp": divergence_bp[(class_family, family, bin_start)],
            "estimated_age_mya": f"{age_mya:.8f}" if age_mya != "" else "",
        })

    family_rows = []
    for class_family, family in sorted(family_bp):
        hit_count = raw[((class_family, family), "hit_count")]
        raw_hit_bp = raw[((class_family, family), "raw_hit_bp")]
        family_rows.append({
            "target": target,
            "group": group,
            "class_family": class_family,
            "pan_family": family,
            "raw_hit_count": hit_count,
            "raw_hit_bp": raw_hit_bp,
            "mean_raw_hit_length_bp": f"{raw_hit_bp / hit_count:.4f}" if hit_count else "",
            "exclusive_bp": family_bp[(class_family, family)],
            "exclusive_fraction_of_repeat_union": (
                f"{family_bp[(class_family, family)] / repeat_union_bp:.8f}" if repeat_union_bp else "0"
            ),
        })
    return {
        "target": target,
        "group": group,
        "divergence_rows": divergence_rows,
        "family_rows": family_rows,
        "category_bp": dict(category_bp),
        "repeat_union_bp": repeat_union_bp,
        "raw_hit_sum_bp": sum(
            value for (key, metric), value in raw.items() if metric == "raw_hit_bp"
        ),
    }


def parse_attributes(text):
    result = {}
    for item in text.strip().split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        result[key] = value
    return result


def union_length(intervals):
    if not intervals:
        return 0
    ordered = sorted(intervals)
    start, end = ordered[0]
    total = 0
    for next_start, next_end in ordered[1:]:
        if next_start > end:
            total += end - start
            start, end = next_start, next_end
        else:
            end = max(end, next_end)
    return total + end - start


def parse_independent_gff_one(task):
    target, group, path = task
    intervals = collections.defaultdict(list)
    total_intervals = collections.defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attributes = parse_attributes(fields[8])
            # EDTA structural child features carry Parent; retain only the
            # top-level structural element or homology annotation.
            if "Parent" in attributes:
                continue
            classification = attributes.get("classification", "")
            if not classification:
                continue
            try:
                start = int(fields[3]) - 1
                end = int(fields[4])
            except ValueError:
                continue
            if end <= start:
                continue
            category = classify_composition(classification)
            intervals[(fields[0], category)].append((start, end))
            total_intervals[fields[0]].append((start, end))
    category_bp = collections.Counter()
    for (_, category), values in intervals.items():
        category_bp[category] += union_length(values)
    total_bp = sum(union_length(values) for values in total_intervals.values())
    return target, group, category_bp, total_bp


def load_family_map(path):
    mapping = {}
    for row in read_tsv(path):
        if row["mapping_status"] == "mapped":
            mapping[row["original_family"]] = (row["pan_family"], row["pan_class_family"])
    return mapping


def load_ltr_rm_index(path):
    index = {}
    for seqid, hits in parse_rmout_hits(path).items():
        ltr_hits = [hit for hit in hits if hit["class_family"].lower().startswith("ltr")]
        ltr_hits.sort(key=lambda hit: (hit["start"], hit["end"]))
        starts = [hit["start"] for hit in ltr_hits]
        prefix_max_end = []
        maximum = -1
        for hit in ltr_hits:
            maximum = max(maximum, hit["end"])
            prefix_max_end.append(maximum)
        index[seqid] = (ltr_hits, starts, prefix_max_end)
    return index


def assign_ltr_by_overlap(seq_index, start, end):
    if not seq_index:
        return "", "", 0
    hits, starts, prefix_max_end = seq_index
    cursor = bisect.bisect_left(starts, end) - 1
    overlap_by_family = collections.Counter()
    class_by_family = {}
    while cursor >= 0:
        if prefix_max_end[cursor] <= start:
            break
        hit = hits[cursor]
        overlap = min(end, hit["end"]) - max(start, hit["start"])
        if overlap > 0:
            overlap_by_family[hit["family"]] += overlap
            class_by_family[hit["family"]] = hit["class_family"]
        cursor -= 1
    if not overlap_by_family:
        return "", "", 0
    family, overlap_bp = max(overlap_by_family.items(), key=lambda value: (value[1], value[0]))
    return family, class_by_family[family], overlap_bp


def parse_intact_ltr_one(task):
    target, group, gff_path, map_path, rmout_path, substitution_rate = task
    family_map = load_family_map(map_path)
    ltr_index = load_ltr_rm_index(rmout_path)
    output = []
    with open(gff_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "repeat_region":
                continue
            attributes = parse_attributes(fields[8])
            class_family = attributes.get("classification", "")
            if not class_family.startswith("LTR/") or "ltr_identity" not in attributes:
                continue
            try:
                identity = float(attributes["ltr_identity"])
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            original_family = attributes.get("Name", "")
            start_zero = start - 1
            pan_family, pan_class, overlap_bp = assign_ltr_by_overlap(
                ltr_index.get(fields[0]), start_zero, end
            )
            if pan_family:
                assignment_method = "maximum_coordinate_overlap"
            else:
                pan_family, pan_class = family_map.get(original_family, ("", ""))
                assignment_method = "original_family_map" if pan_family else "unassigned"
            overlap_fraction = overlap_bp / (end - start + 1)
            divergence = max(0.0, 1.0 - identity)
            age_mya = ""
            if substitution_rate:
                age_mya = divergence / (2 * substitution_rate) / 1_000_000
            output.append({
                "target": target,
                "group": group,
                "seqid": fields[0],
                "start": start,
                "end": end,
                "length_bp": end - start + 1,
                "strand": fields[6],
                "original_family": original_family,
                "original_class_family": class_family,
                "pan_family": pan_family,
                "pan_class_family": pan_class,
                "assignment_method": assignment_method,
                "overlap_bp": overlap_bp,
                "overlap_fraction": f"{overlap_fraction:.8f}",
                "ltr_identity": f"{identity:.8f}",
                "ltr_divergence": f"{divergence:.8f}",
                "estimated_age_mya": f"{age_mya:.8f}" if age_mya != "" else "",
            })
    return output


def read_pan_library_audit(path):
    audit = {}
    header = None
    sequence = []

    def store_record(raw_header, chunks):
        if raw_header is None:
            return
        seq = "".join(chunks).upper()
        family = raw_header.split()[0].split("#", 1)[0]
        counts = collections.Counter(base for base in seq if base in "ACGT")
        resolved = sum(counts.values())
        entropy = 0.0
        if resolved:
            for count in counts.values():
                proportion = count / resolved
                entropy -= proportion * math.log2(proportion)
        dominant_fraction = max(counts.values(), default=0) / resolved if resolved else 0.0
        audit[family] = {
            "pan_family": family,
            "raw_header": raw_header,
            "sequence_length_bp": len(seq),
            "resolved_acgt_fraction": f"{resolved / len(seq):.8f}" if seq else "0",
            "gc_fraction": f"{(counts['G'] + counts['C']) / resolved:.8f}" if resolved else "",
            "base_entropy_bits": f"{entropy:.8f}" if resolved else "",
            "dominant_base_fraction": f"{dominant_fraction:.8f}" if resolved else "",
            "low_complexity_flag": int(entropy < 1.2 or dominant_fraction > 0.80) if resolved else 1,
        }

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store_record(header, sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
    store_record(header, sequence)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-dir", required=True, type=Path)
    parser.add_argument("--discovery-dir", required=True, type=Path)
    parser.add_argument("--pan-library", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--cumulative-threshold", type=float, default=0.80)
    parser.add_argument("--divergence-bin-width", type=float, default=1.0)
    parser.add_argument("--substitution-rate", type=float)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    if not 0 < args.cumulative_threshold <= 1:
        parser.error("--cumulative-threshold must be in (0, 1]")
    if args.divergence_bin_width <= 0:
        parser.error("--divergence-bin-width must be positive")
    if args.substitution_rate is not None and args.substitution_rate <= 0:
        parser.error("--substitution-rate must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    total_rows = read_tsv(args.summary_dir / "prunus_pan_te_total_summary.tsv")
    targets_rows = read_tsv(args.summary_dir / "standardized_mask_targets.tsv")
    validate_targets(targets_rows, "standardized mask targets")
    if not args.pan_library.is_file():
        raise FileNotFoundError(args.pan_library)
    library_audit = read_pan_library_audit(args.pan_library)

    mapping_qc = build_mapping_qc(args.summary_dir / "original_to_pan_family")

    rm_tasks = []
    rmout_by_target = {}
    for row in targets_rows:
        rmout_by_target[row["target"]] = row["repeatmasker_out"]
        rm_tasks.append((
            row["target"], row["group"], row["repeatmasker_out"],
            args.divergence_bin_width, args.substitution_rate,
        ))
    with ProcessPoolExecutor(max_workers=max(1, args.threads)) as pool:
        exclusive_results = list(pool.map(parse_rmout_one, rm_tasks))
    rm_rows = [
        row for result in exclusive_results for row in result["divergence_rows"]
    ]
    exclusive_family_rows = [
        row for result in exclusive_results for row in result["family_rows"]
    ]
    composition, composition_qc = build_composition(total_rows, exclusive_results)
    drivers = build_driver_table(
        exclusive_family_rows, args.cumulative_threshold, library_audit
    )

    write_tsv(
        args.output_dir / "te_composition.tsv",
        ["target", "group", "subgenome", "display_group", "category", "bp", "genome_fraction", "genome_bp", "repeat_union_bp"],
        composition,
    )
    driver_fields = [
        "positive_rank", "driver_selected", "class_family", "pan_family",
        "pzi_a_mean_bp", "psp_a_mean_bp", "delta_a_bp",
        "pzi_b_mean_bp", "psp_b_mean_bp", "delta_b_bp",
        "pzi_mean_bp", "psp_mean_bp", "delta_pzi_psp_mean_bp", "positive_excess_bp",
        "cumulative_positive_excess_fraction", "ppa_pcam_ppe_median_bp",
        "expanded_in_both_subgenomes", "consistent_across_pzi_haplotypes",
        "library_sequence_length_bp", "library_gc_fraction",
        "library_base_entropy_bits", "library_low_complexity_flag",
    ] + [f"{target}_bp" for target in TARGET_ORDER]
    write_tsv(args.output_dir / "te_driver_families.tsv", driver_fields, drivers)
    write_tsv(
        args.output_dir / "te_library_mapping_qc.tsv",
        ["target", "group", "original_family_count", "mapped_family_count", "unmapped_family_count", "mapped_fraction", "median_mapped_identity", "mapped_identity_ge_0.80_fraction"],
        mapping_qc,
    )
    write_tsv(
        args.output_dir / "te_composition_qc.tsv",
        ["target", "genome_bp", "repeat_union_bp", "raw_hit_sum_bp",
         "raw_hit_sum_to_repeat_union_ratio", "exclusive_assignment_sum_bp",
         "exclusive_sum_to_repeat_union_ratio"],
        composition_qc,
    )
    write_tsv(
        args.output_dir / "te_exclusive_family_summary.tsv",
        ["target", "group", "class_family", "pan_family", "raw_hit_count",
         "raw_hit_bp", "mean_raw_hit_length_bp", "exclusive_bp",
         "exclusive_fraction_of_repeat_union"],
        exclusive_family_rows,
    )
    write_tsv(
        args.output_dir / "te_pan_library_audit.tsv",
        ["pan_family", "raw_header", "sequence_length_bp", "resolved_acgt_fraction",
         "gc_fraction", "base_entropy_bits", "dominant_base_fraction", "low_complexity_flag"],
        [library_audit[key] for key in sorted(library_audit)],
    )
    write_tsv(
        args.output_dir / "te_repeatmasker_divergence.tsv",
        ["target", "group", "class_family", "pan_family", "divergence_bin_start_pct",
         "divergence_bin_mid_pct", "hit_count", "exclusive_bp", "estimated_age_mya"],
        rm_rows,
    )

    intact_tasks = []
    mapping_dir = args.summary_dir / "original_to_pan_family"
    for target in TARGET_ORDER:
        intact_tasks.append((
            target,
            TARGET_META[target][0],
            args.discovery_dir / target / f"{target}.fa.mod.EDTA.intact.gff3",
            mapping_dir / f"{target}.original_to_pan_family.tsv",
            rmout_by_target[target],
            args.substitution_rate,
        ))
    for task in intact_tasks:
        if not Path(task[2]).is_file():
            raise FileNotFoundError(task[2])
    with ProcessPoolExecutor(max_workers=max(1, args.threads)) as pool:
        intact_nested = list(pool.map(parse_intact_ltr_one, intact_tasks))
    intact_rows = [row for target_rows in intact_nested for row in target_rows]
    write_tsv(
        args.output_dir / "te_intact_ltr.tsv",
        ["target", "group", "seqid", "start", "end", "length_bp", "strand",
         "original_family", "original_class_family", "pan_family", "pan_class_family",
         "assignment_method", "overlap_bp", "overlap_fraction", "ltr_identity",
         "ltr_divergence", "estimated_age_mya"],
        intact_rows,
    )

    independent_tasks = [
        (
            target,
            TARGET_META[target][0],
            args.discovery_dir / target / f"{target}.fa.mod.EDTA.TEanno.gff3",
        )
        for target in TARGET_ORDER
    ]
    for task in independent_tasks:
        if not Path(task[2]).is_file():
            raise FileNotFoundError(task[2])
    with ProcessPoolExecutor(max_workers=max(1, args.threads)) as pool:
        independent_results = list(pool.map(parse_independent_gff_one, independent_tasks))
    pan_by_target_category = collections.Counter()
    for row in composition:
        if row["category"] != "Non-TE":
            pan_by_target_category[(row["target"], row["category"])] = int(row["bp"])
    sensitivity_rows = []
    for target, group, category_bp, independent_total in independent_results:
        for category in [value for value in COMPOSITION_ORDER if value != "Non-TE"]:
            pan_bp = pan_by_target_category[(target, category)]
            standalone_bp = category_bp[category]
            sensitivity_rows.append({
                "target": target,
                "group": group,
                "category": category,
                "pan_union_bp": pan_bp,
                "independent_edta_union_bp": standalone_bp,
                "pan_minus_independent_bp": pan_bp - standalone_bp,
                "pan_to_independent_ratio": f"{pan_bp / standalone_bp:.8f}" if standalone_bp else "",
                "independent_total_union_bp": independent_total,
            })
    write_tsv(
        args.output_dir / "te_pan_vs_independent.tsv",
        ["target", "group", "category", "pan_union_bp", "independent_edta_union_bp", "pan_minus_independent_bp", "pan_to_independent_ratio", "independent_total_union_bp"],
        sensitivity_rows,
    )

    manifest = {
        "formal_targets": TARGET_ORDER,
        "driver_definition": "positive mean PZI-minus-PSP exclusive pan-family bp difference",
        "driver_cumulative_threshold": args.cumulative_threshold,
        "divergence_bin_width_percent": args.divergence_bin_width,
        "substitution_rate_substitutions_per_site_per_year": args.substitution_rate,
        "absolute_age_policy": "blank unless a positive substitution rate is explicitly supplied",
        "composition_policy": (
            "every RepeatMasker-covered base assigned once to the overlapping hit with highest "
            "SW-score density, then lower divergence, higher SW score and longer span"
        ),
        "repeatmasker_age_source": "RepeatMasker .out percent-divergence; not CpG-corrected Kimura distance",
        "intact_ltr_age_source": "1 - EDTA ltr_identity",
        "intact_ltr_pan_family_assignment": (
            "maximum genomic coordinate overlap with standardized LTR RepeatMasker hits; "
            "original-family mapping used only as fallback"
        ),
        "pan_library_audit_policy": (
            "sequence length, GC, ACGT resolution and base entropy; low complexity flagged at "
            "entropy <1.2 bits or dominant base fraction >0.80"
        ),
        "sensitivity_policy": "top-level independent EDTA GFF3 intervals versus standardized pan-library intervals",
    }
    with open(args.output_dir / "postprocess_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    selected = sum(int(row["driver_selected"]) for row in drivers)
    print("=== Prunus pan-TE post-processing completed ===")
    print(f"Selected driver families: {selected}")
    print(f"RepeatMasker divergence rows: {len(rm_rows)}")
    print(f"Intact LTR records: {len(intact_rows)}")
    print(f"Pan-versus-independent rows: {len(sensitivity_rows)}")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
