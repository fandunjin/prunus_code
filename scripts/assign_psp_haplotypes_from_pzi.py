#!/usr/bin/env python3
"""Assign PSP chromosome-scale scaffolds to PZI-compatible A1/A2/B1/B2 names.

The assignment is constrained in two ways. First, the 32 PSP chromosomes are
split into two groups of 16 using their hifiasm hap1/hap2 provenance. Second,
each phase group is matched one-to-one to either the 16 PZI *1 or the 16 PZI
*2 chromosomes. Both global phase orientations are evaluated and the higher
scoring solution is retained.
"""

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


# PZI targets in the revised partition are C1A_A1 ... C8B_B2.  The legacy
# C01_A1 spelling remains accepted only so archived assignments can be read.
# The optional middle A/B token is checked by parse_pzi_name() so that a
# malformed target such as C1A_B1 cannot silently enter the assignment.
PZI_NAME = re.compile(r"^C0?([1-8])_?([AB])(?:_?([AB])?([12]))$", re.IGNORECASE)
HAP_ORDER = ("hapA1", "hapA2", "hapB1", "hapB2")
RC_TABLE = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)


def parse_pzi_name(value):
    """Return (C-group, subgenome, phase) for a revised PZI target."""
    match = PZI_NAME.fullmatch(value)
    if match is None:
        return None
    group, subgenome, phase_subgenome, phase = match.groups()
    subgenome = subgenome.upper()
    if phase_subgenome and phase_subgenome.upper() != subgenome:
        return None
    return int(group), subgenome, int(phase)


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise SystemExit(f"ERROR: duplicate FASTA ID: {name}")
                    records[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise SystemExit(f"ERROR: sequence before FASTA header in {path}")
                chunks.append(line)
        if name is not None:
            if name in records:
                raise SystemExit(f"ERROR: duplicate FASTA ID: {name}")
            records[name] = "".join(chunks)
    return records


def write_fasta(path, rows, records):
    with open(path, "w", encoding="utf-8") as out:
        for row in rows:
            seq = records[row["old_seqid"]]
            if row["orientation"] == "reverse_complement":
                seq = seq.translate(RC_TABLE)[::-1]
            out.write(f">{row['new_seqid']}\n")
            for start in range(0, len(seq), 80):
                out.write(seq[start:start + 80] + "\n")


def interval_union_bp(intervals):
    if not intervals:
        return 0
    total = 0
    start, end = sorted(intervals)[0]
    for next_start, next_end in sorted(intervals)[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def parse_provenance_paf(path, source, chromosome_ids, min_mapq, min_aln_len, store):
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            target = fields[5]
            if target not in chromosome_ids:
                continue
            mapq = int(fields[11])
            aln_len = int(fields[10])
            if mapq < min_mapq or aln_len < min_aln_len:
                continue
            store[target][source].append((int(fields[7]), int(fields[8])))


def infer_phase_groups(records, hap1_paf, hap2_paf, min_mapq, min_aln_len):
    intervals = defaultdict(lambda: {"hap1": [], "hap2": []})
    chromosome_ids = set(records)
    parse_provenance_paf(
        hap1_paf, "hap1", chromosome_ids, min_mapq, min_aln_len, intervals
    )
    parse_provenance_paf(
        hap2_paf, "hap2", chromosome_ids, min_mapq, min_aln_len, intervals
    )

    rows = []
    for seqid, seq in records.items():
        hap1_bp = interval_union_bp(intervals[seqid]["hap1"])
        hap2_bp = interval_union_bp(intervals[seqid]["hap2"])
        total = hap1_bp + hap2_bp
        if total == 0:
            dominant = "unknown"
            dominance = 0.0
        elif hap1_bp >= hap2_bp:
            dominant = "hap1"
            dominance = hap1_bp / total
        else:
            dominant = "hap2"
            dominance = hap2_bp / total
        phase_signal = math.log((hap1_bp + 1.0) / (hap2_bp + 1.0))
        rows.append(
            {
                "old_seqid": seqid,
                "length": len(seq),
                "hap1_covered_bp": hap1_bp,
                "hap2_covered_bp": hap2_bp,
                "dominant_source": dominant,
                "provenance_dominance": dominance,
                "phase_signal": phase_signal,
            }
        )

    # Enforce the biological/assembly constraint of 16 chromosomes from each
    # hifiasm phase while retaining the raw dominance call for review.
    ranked = sorted(rows, key=lambda row: (-row["phase_signal"], row["old_seqid"]))
    hap1_ids = {row["old_seqid"] for row in ranked[:16]}
    for row in rows:
        row["assigned_source"] = "hap1" if row["old_seqid"] in hap1_ids else "hap2"
        flags = []
        if row["dominant_source"] == "unknown":
            flags.append("no_provenance_alignment")
        elif row["assigned_source"] != row["dominant_source"]:
            flags.append("phase_forced_by_16_16_constraint")
        if row["provenance_dominance"] < 0.60:
            flags.append("weak_provenance")
        row["provenance_review"] = ";".join(flags) if flags else "pass"
    return rows


def parse_pzi_paf(path, query_ids, min_mapq, min_aln_len):
    hits = defaultdict(list)
    targets = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            target = fields[5]
            if query not in query_ids or parse_pzi_name(target) is None:
                continue
            mapq = int(fields[11])
            aln_len = int(fields[10])
            if mapq < min_mapq or aln_len < min_aln_len:
                continue
            divergence = None
            identity_method = "matches_per_block"
            for tag in fields[12:]:
                if tag.startswith("de:f:"):
                    divergence = float(tag[5:])
                    identity_method = "1-de"
                    break
            if divergence is None:
                for tag in fields[12:]:
                    if tag.startswith("dv:f:"):
                        divergence = float(tag[5:])
                        identity_method = "1-dv"
                        break
            if divergence is None:
                hit_identity = int(fields[9]) / aln_len if aln_len else 0.0
            else:
                hit_identity = max(0.0, min(1.0, 1.0 - divergence))
            targets.add(target)
            hits[(query, target)].append(
                {
                    "qstart": int(fields[2]),
                    "qend": int(fields[3]),
                    "strand": fields[4],
                    "tstart": int(fields[7]),
                    "tend": int(fields[8]),
                    "matches": int(fields[9]),
                    "aln_len": aln_len,
                    "identity": hit_identity,
                    "identity_method": identity_method,
                }
            )

    expected = {
        f"C{group}{subgenome}_{subgenome}{suffix}"
        for group in range(1, 9)
        for subgenome in "AB"
        for suffix in (1, 2)
    }
    missing = sorted(expected - targets)
    if len(targets) < 24:
        raise SystemExit(
            f"ERROR: retained alignments cover only {len(targets)}/32 PZI targets; "
            "check the PZI names and minimap2 input order. Missing: " + ", ".join(missing)
        )
    if missing:
        print(
            "WARNING: no retained PSP alignment to PZI target(s); these will be "
            "assigned only under the global one-to-one constraint and flagged: "
            + ", ".join(missing)
        )

    scores = {}
    for key, pair_hits in hits.items():
        query_bp = interval_union_bp(
            [(hit["qstart"], hit["qend"]) for hit in pair_hits]
        )
        target_bp = interval_union_bp(
            [(hit["tstart"], hit["tend"]) for hit in pair_hits]
        )
        total_aln = sum(hit["aln_len"] for hit in pair_hits)
        identity = (
            sum(hit["identity"] * hit["aln_len"] for hit in pair_hits) / total_aln
            if total_aln
            else 0.0
        )
        identity_method = ",".join(
            sorted({hit["identity_method"] for hit in pair_hits})
        )
        aligned_bp = min(query_bp, target_bp)
        plus = sum(hit["matches"] for hit in pair_hits if hit["strand"] == "+")
        minus = sum(hit["matches"] for hit in pair_hits if hit["strand"] == "-")
        scores[key] = {
            "score": aligned_bp * identity,
            "aligned_bp": aligned_bp,
            "query_covered_bp": query_bp,
            "target_covered_bp": target_bp,
            "identity": identity,
            "identity_method": identity_method,
            "plus_support": plus,
            "minus_support": minus,
        }
    return scores


def target_names(suffix):
    return [
        f"C{group}{subgenome}_{subgenome}{suffix}"
        for group in range(1, 9)
        for subgenome in "AB"
    ]


def maximum_assignment(queries, targets, pair_scores):
    """Exact maximum-weight 16 x 16 assignment using bitmask dynamic programming."""
    queries = sorted(queries)
    targets = sorted(targets)
    if len(queries) != 16 or len(targets) != 16:
        raise SystemExit(
            f"ERROR: assignment requires 16 queries and 16 targets; got "
            f"{len(queries)} and {len(targets)}"
        )
    current = {0: 0.0}
    parents = []
    for query in queries:
        next_scores = {}
        next_parent = {}
        for mask, cumulative in current.items():
            for target_index, target in enumerate(targets):
                bit = 1 << target_index
                if mask & bit:
                    continue
                new_mask = mask | bit
                value = cumulative + pair_scores.get((query, target), {}).get("score", 0.0)
                if value > next_scores.get(new_mask, -1.0):
                    next_scores[new_mask] = value
                    next_parent[new_mask] = (mask, target_index)
        parents.append(next_parent)
        current = next_scores

    full_mask = (1 << len(targets)) - 1
    total = current[full_mask]
    assignments = []
    mask = full_mask
    for query_index in range(len(queries) - 1, -1, -1):
        previous_mask, target_index = parents[query_index][mask]
        assignments.append((queries[query_index], targets[target_index]))
        mask = previous_mask
    assignments.reverse()
    return assignments, total


def evaluate_orientation(phase_rows, pair_scores, hap1_suffix):
    hap2_suffix = 2 if hap1_suffix == 1 else 1
    hap1_queries = [row["old_seqid"] for row in phase_rows if row["assigned_source"] == "hap1"]
    hap2_queries = [row["old_seqid"] for row in phase_rows if row["assigned_source"] == "hap2"]
    assignment1, score1 = maximum_assignment(hap1_queries, target_names(hap1_suffix), pair_scores)
    assignment2, score2 = maximum_assignment(hap2_queries, target_names(hap2_suffix), pair_scores)
    return {
        "hap1_suffix": hap1_suffix,
        "hap2_suffix": hap2_suffix,
        "assignments": assignment1 + assignment2,
        "hap1_score": score1,
        "hap2_score": score2,
        "total_score": score1 + score2,
    }


def build_assignment_rows(phase_rows, selected, pair_scores, min_assignment_bp):
    phase_by_id = {row["old_seqid"]: row for row in phase_rows}
    chosen = dict(selected["assignments"])
    rows = []
    for old_seqid, new_seqid in chosen.items():
        phase = phase_by_id[old_seqid]
        metrics = pair_scores.get((old_seqid, new_seqid), {})
        allowed_suffix = int(new_seqid[-1])
        allowed_targets = target_names(allowed_suffix)
        ranked = sorted(
            allowed_targets,
            key=lambda target: (
                -pair_scores.get((old_seqid, target), {}).get("score", 0.0),
                target,
            ),
        )
        chosen_rank = ranked.index(new_seqid) + 1
        best_target = ranked[0]
        best_score = pair_scores.get((old_seqid, best_target), {}).get("score", 0.0)
        alternatives = [
            pair_scores.get((old_seqid, target), {}).get("score", 0.0)
            for target in allowed_targets
            if target != new_seqid
        ]
        second_score = max(alternatives) if alternatives else 0.0
        chosen_score = metrics.get("score", 0.0)
        margin = (chosen_score - second_score) / chosen_score if chosen_score else -1.0
        plus = metrics.get("plus_support", 0)
        minus = metrics.get("minus_support", 0)
        orientation = "reverse_complement" if minus > plus else "forward"
        strand_total = plus + minus
        strand_dominance = max(plus, minus) / strand_total if strand_total else 0.0

        parsed = parse_pzi_name(new_seqid)
        if parsed is None:
            raise SystemExit(f"ERROR: invalid PZI target generated by assignment: {new_seqid}")
        c_group, subgenome, suffix = parsed
        hap_set = f"hap{subgenome}{suffix}"
        flags = []
        if phase["provenance_review"] != "pass":
            flags.extend(phase["provenance_review"].split(";"))
        if metrics.get("aligned_bp", 0) < min_assignment_bp:
            flags.append("low_pzi_aligned_bp")
        if metrics.get("identity", 0.0) < 0.80:
            flags.append("low_pzi_identity")
        if strand_dominance < 0.60:
            flags.append("ambiguous_orientation")
        if chosen_rank > 1:
            flags.append("global_assignment_not_row_best")

        rows.append(
            {
                **phase,
                "new_seqid": new_seqid,
                "hap_set": hap_set,
                "c_group": f"C{c_group}",
                "subgenome": subgenome,
                "phase_suffix": suffix,
                "pzi_target": new_seqid,
                "pair_score": chosen_score,
                "aligned_bp": metrics.get("aligned_bp", 0),
                "query_covered_bp": metrics.get("query_covered_bp", 0),
                "target_covered_bp": metrics.get("target_covered_bp", 0),
                "identity": metrics.get("identity", 0.0),
                "identity_method": metrics.get("identity_method", ""),
                "plus_support": plus,
                "minus_support": minus,
                "strand_dominance": strand_dominance,
                "orientation": orientation,
                "row_best_target": best_target,
                "row_best_score": best_score,
                "chosen_target_rank": chosen_rank,
                "second_best_score": second_score,
                "chosen_vs_second_margin": margin,
                "review_flag": ";".join(sorted(set(flags))) if flags else "pass",
            }
        )
    return sorted(rows, key=lambda row: (HAP_ORDER.index(row["hap_set"]), row["new_seqid"]))


def write_tsv(path, rows, fields):
    with open(path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            for key in (
                "provenance_dominance",
                "phase_signal",
                "score",
                "pair_score",
                "identity",
                "strand_dominance",
                "row_best_score",
                "second_best_score",
                "chosen_vs_second_margin",
            ):
                if key in formatted and isinstance(formatted[key], float):
                    formatted[key] = f"{formatted[key]:.6f}"
            writer.writerow(formatted)


def validate_and_write(out_dir, assignment_rows, phase_orientation_margin):
    validations = []

    def add(check, observed, expected, passed):
        validations.append(
            {
                "check": check,
                "observed": observed,
                "expected": expected,
                "status": "pass" if passed else "fail",
            }
        )

    names = [row["new_seqid"] for row in assignment_rows]
    sources = [row["old_seqid"] for row in assignment_rows]
    add("total_chromosomes", len(names), 32, len(names) == 32)
    add("unique_source_chromosomes", len(set(sources)), 32, len(set(sources)) == 32)
    add("unique_new_names", len(set(names)), 32, len(set(names)) == 32)
    for hap_set in HAP_ORDER:
        count = sum(row["hap_set"] == hap_set for row in assignment_rows)
        add(f"{hap_set}_chromosomes", count, 8, count == 8)
    expected_names = {
        f"C{group}{subgenome}_{subgenome}{suffix}"
        for group in range(1, 9)
        for subgenome in "AB"
        for suffix in (1, 2)
    }
    add("complete_C1_C8_A1A2B1B2", len(set(names) & expected_names), 32, set(names) == expected_names)
    review_count = sum(row["review_flag"] != "pass" for row in assignment_rows)
    add("review_flags", review_count, 0, review_count == 0)
    add(
        "global_phase_orientation_margin",
        f"{phase_orientation_margin:.6f}",
        ">=0.020000",
        phase_orientation_margin >= 0.02,
    )
    write_tsv(
        out_dir / "psp_assignment_validation.tsv",
        validations,
        ["check", "observed", "expected", "status"],
    )
    hard_checks = {
        "total_chromosomes",
        "unique_source_chromosomes",
        "unique_new_names",
        "hapA1_chromosomes",
        "hapA2_chromosomes",
        "hapB1_chromosomes",
        "hapB2_chromosomes",
        "complete_C1_C8_A1A2B1B2",
    }
    hard_failures = [
        row
        for row in validations
        if row["check"] in hard_checks and row["status"] == "fail"
    ]
    if hard_failures:
        raise SystemExit("ERROR: structural assignment validation failed; inspect psp_assignment_validation.tsv")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chr32-fasta", required=True)
    parser.add_argument("--hap1-paf", required=True)
    parser.add_argument("--hap2-paf", required=True)
    parser.add_argument("--pzi-paf", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--provenance-min-mapq", type=int, default=10)
    parser.add_argument("--provenance-min-aln-len", type=int, default=10000)
    parser.add_argument("--pzi-min-mapq", type=int, default=0)
    parser.add_argument("--pzi-min-aln-len", type=int, default=10000)
    parser.add_argument("--min-assignment-bp", type=int, default=1000000)
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    genome_dir = out_dir / "genome"
    review_dir = out_dir / "review"
    genome_dir.mkdir(parents=True, exist_ok=True)
    review_dir.mkdir(parents=True, exist_ok=True)

    records = read_fasta(args.chr32_fasta)
    if len(records) != 32:
        raise SystemExit(f"ERROR: expected 32 chromosome records, observed {len(records)}")

    phase_rows = infer_phase_groups(
        records,
        args.hap1_paf,
        args.hap2_paf,
        args.provenance_min_mapq,
        args.provenance_min_aln_len,
    )
    pair_scores = parse_pzi_paf(
        args.pzi_paf,
        set(records),
        args.pzi_min_mapq,
        args.pzi_min_aln_len,
    )
    orientations = [
        evaluate_orientation(phase_rows, pair_scores, 1),
        evaluate_orientation(phase_rows, pair_scores, 2),
    ]
    selected = max(orientations, key=lambda row: row["total_score"])
    alternative_score = min(row["total_score"] for row in orientations)
    phase_orientation_margin = (
        (selected["total_score"] - alternative_score) / selected["total_score"]
        if selected["total_score"]
        else 0.0
    )
    assignments = build_assignment_rows(
        phase_rows, selected, pair_scores, args.min_assignment_bp
    )

    provenance_fields = [
        "old_seqid", "length", "hap1_covered_bp", "hap2_covered_bp",
        "dominant_source", "provenance_dominance", "phase_signal",
        "assigned_source", "provenance_review",
    ]
    write_tsv(review_dir / "psp_chr32_provenance.tsv", phase_rows, provenance_fields)

    pair_rows = []
    all_targets = target_names(1) + target_names(2)
    for query in sorted(records):
        for target in sorted(all_targets):
            metrics = pair_scores.get(
                (query, target),
                {
                    "score": 0.0,
                    "aligned_bp": 0,
                    "query_covered_bp": 0,
                    "target_covered_bp": 0,
                    "identity": 0.0,
                    "identity_method": "",
                    "plus_support": 0,
                    "minus_support": 0,
                },
            )
            pair_rows.append({"old_seqid": query, "pzi_target": target, **metrics})
    write_tsv(
        review_dir / "psp_pzi_pair_scores.tsv",
        pair_rows,
        [
            "old_seqid", "pzi_target", "score", "aligned_bp",
            "query_covered_bp", "target_covered_bp", "identity", "identity_method",
            "plus_support", "minus_support",
        ],
    )

    orientation_rows = []
    for row in orientations:
        orientation_rows.append(
            {
                "hap1_suffix": row["hap1_suffix"],
                "hap2_suffix": row["hap2_suffix"],
                "hap1_score": row["hap1_score"],
                "hap2_score": row["hap2_score"],
                "total_score": row["total_score"],
                "selected": "yes" if row is selected else "no",
                "selected_vs_alternative_margin": phase_orientation_margin,
            }
        )
    write_tsv(
        review_dir / "psp_phase_orientation_scores.tsv",
        orientation_rows,
        [
            "hap1_suffix", "hap2_suffix", "hap1_score", "hap2_score",
            "total_score", "selected", "selected_vs_alternative_margin",
        ],
    )

    assignment_fields = [
        "old_seqid", "new_seqid", "hap_set", "c_group", "subgenome",
        "phase_suffix", "assigned_source", "length", "hap1_covered_bp",
        "hap2_covered_bp", "dominant_source", "provenance_dominance",
        "pzi_target", "pair_score", "aligned_bp", "query_covered_bp",
        "target_covered_bp", "identity", "identity_method", "plus_support", "minus_support",
        "strand_dominance", "orientation", "row_best_target",
        "row_best_score", "chosen_target_rank", "second_best_score",
        "chosen_vs_second_margin", "review_flag",
    ]
    write_tsv(out_dir / "psp_A1A2B1B2.rename_map.tsv", assignments, assignment_fields)
    write_tsv(
        review_dir / "psp_assignment_review.tsv",
        [row for row in assignments if row["review_flag"] != "pass"],
        assignment_fields,
    )

    for hap_set in HAP_ORDER:
        hap_rows = [row for row in assignments if row["hap_set"] == hap_set]
        write_fasta(genome_dir / f"psp_{hap_set}.fa", hap_rows, records)
    combined_rows = sorted(assignments, key=lambda row: row["new_seqid"])
    write_fasta(genome_dir / "psp_A1A2B1B2.renamed.chr32.fa", combined_rows, records)
    validate_and_write(out_dir, assignments, phase_orientation_margin)

    review_count = sum(row["review_flag"] != "pass" for row in assignments)
    print(f"Selected global orientation: hap1 -> suffix {selected['hap1_suffix']}")
    print(f"Selected total alignment score: {selected['total_score']:.0f}")
    print(f"Global phase-orientation margin: {phase_orientation_margin:.6f}")
    if phase_orientation_margin < 0.02:
        print(
            "WARNING: suffix 1/2 is weakly resolved across species; treat A1/A2 "
            "and B1/B2 as operational phase labels rather than distinct ancestry."
        )
    print(f"Chromosomes carrying review flags: {review_count}")
    print(f"Output directory: {out_dir}")


if __name__ == "__main__":
    main()
