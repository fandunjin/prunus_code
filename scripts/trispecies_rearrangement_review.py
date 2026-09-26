#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


PCAM_MAIN_CHR = [f"CM11428{i}.1" for i in range(0, 8)]
PPA_MAIN_CHR = [f"SUPER_{i}" for i in range(1, 17)]


def natural_key(text: str):
    parts = []
    buf = ""
    is_digit = False
    for ch in text:
        if ch.isdigit() != is_digit and buf:
            parts.append(int(buf) if is_digit else buf)
            buf = ch
            is_digit = ch.isdigit()
        else:
            buf += ch
            is_digit = ch.isdigit()
    if buf:
        parts.append(int(buf) if is_digit else buf)
    return parts


def read_fasta_headers(path: Path):
    headers = []
    opener = path.open
    if path.suffix == ".gz":
        import gzip
        opener = lambda: gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    with opener() as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_paf(path: Path, min_mapq: int, min_aln_len: int):
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            row = {
                "query": fields[0],
                "query_len": int(fields[1]),
                "query_start": int(fields[2]),
                "query_end": int(fields[3]),
                "strand": fields[4],
                "target": fields[5],
                "target_len": int(fields[6]),
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "matches": int(fields[9]),
                "aln_len": int(fields[10]),
                "mapq": int(fields[11]),
            }
            if row["mapq"] >= min_mapq and row["aln_len"] >= min_aln_len:
                rows.append(row)
    return rows


def coverage_union(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    cur_s, cur_e = intervals[0]
    total = 0
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    total += cur_e - cur_s
    return total


def summarize_pairs(rows):
    by_query = defaultdict(lambda: defaultdict(lambda: {"intervals": [], "aln": 0, "matches": 0, "plus": 0, "minus": 0}))
    qlen = {}
    tlen = {}
    qorder = []
    torder = []
    for row in rows:
        q = row["query"]
        t = row["target"]
        if q not in qorder:
            qorder.append(q)
        if t not in torder:
            torder.append(t)
        qlen[q] = row["query_len"]
        tlen[t] = row["target_len"]
        stats = by_query[q][t]
        stats["intervals"].append((row["query_start"], row["query_end"]))
        stats["aln"] += row["aln_len"]
        stats["matches"] += row["matches"]
        if row["strand"] == "+":
            stats["plus"] += row["aln_len"]
        else:
            stats["minus"] += row["aln_len"]
    summary = []
    best = []
    for q in sorted(qorder, key=natural_key):
        best_t = ""
        best_cov = -1
        for t, stats in by_query.get(q, {}).items():
            cov = coverage_union(stats["intervals"])
            if cov > best_cov:
                best_cov = cov
                best_t = t
        for t in sorted(torder, key=natural_key):
            stats = by_query.get(q, {}).get(t)
            if not stats:
                continue
            cov = coverage_union(stats["intervals"])
            aln = stats["aln"]
            ident = stats["matches"] / aln if aln else 0.0
            orient = "+" if stats["plus"] >= stats["minus"] else "-"
            row = {
                "query_chr": q,
                "target_chr": t,
                "query_length_bp": str(qlen.get(q, "")),
                "target_length_bp": str(tlen.get(t, "")),
                "aligned_bp": str(aln),
                "covered_bp": str(cov),
                "coverage_ratio": f"{(cov / qlen[q]):.6f}" if q in qlen and qlen[q] else "",
                "identity": f"{ident:.6f}",
                "orientation": orient,
                "best_hit_for_query": "yes" if t == best_t else "no",
            }
            summary.append(row)
            if t == best_t:
                best.append(row)
        if not any(r["query_chr"] == q and r["best_hit_for_query"] == "yes" for r in summary):
            empty = {
                "query_chr": q,
                "target_chr": "",
                "query_length_bp": str(qlen.get(q, "")),
                "target_length_bp": "",
                "aligned_bp": "0",
                "covered_bp": "0",
                "coverage_ratio": "0.000000",
                "identity": "0.000000",
                "orientation": "",
                "best_hit_for_query": "yes",
            }
            summary.append(empty)
            best.append(empty)
    return summary, best


def build_axis_colors(summary_rows, pcam_style="pcam"):
    palette = ["#F6AD55", "#68D391", "#63B3ED", "#F687B3", "#4FD1C5", "#ED8936", "#9F7AEA", "#ECC94B", "#F56565", "#4299E1"]
    colors = {}
    used = {}
    idx = 0
    for row in summary_rows:
        if row["best_hit_for_query"] != "yes":
            continue
        label = row["target_chr"]
        query = row["query_chr"]
        if label not in used:
            used[label] = palette[idx % len(palette)]
            idx += 1
        colors[query] = used[label]
        colors[label] = used[label]
    return colors


def classify_pcam_event(pcam_chr, members, pzi_best, pzi_self_pairs, pzi_lengths, pcam_lengths):
    member_set = set(members)
    hits = [row for row in pzi_best if row["best_hit_for_query"] == "yes" and row["target_chr"] == pcam_chr]
    strong_members = [row["query_chr"] for row in hits]
    strong_set = set(strong_members)
    member_count = len(member_set)
    strong_count = len(strong_set)
    pcam_len = pcam_lengths.get(pcam_chr, 0)
    pzi_total = sum(pzi_lengths.get(m, 0) for m in member_set)
    length_ratio = (pzi_total / pcam_len) if pcam_len else 0.0

    self_synteny_bp = 0
    for a in member_set:
        for b in member_set:
            if a != b:
                self_synteny_bp += int(pzi_self_pairs.get((a, b), 0))

    if strong_count >= 4 and member_count <= 4:
        event = "balanced_quartet"
        eej_ncf_rct = "none"
        review = "pass"
        suggestion = "Keep as balanced quartet candidate."
    elif strong_count >= 4 and member_count > 4:
        event = "complex_rearrangement"
        eej_ncf_rct = "NCF_like"
        review = "review"
        suggestion = "More than four PZI chromosomes map to one PCAM chromosome; inspect for nested fusion or unresolved split."
    elif strong_count == 3:
        event = "possible_break_or_split"
        eej_ncf_rct = "EEJ_like"
        review = "review"
        suggestion = "Three-member support suggests breakpoint/split or missing partner; check dotplot breakpoints and length imbalance."
    elif strong_count == 2:
        if member_count == 2:
            event = "possible_fusion_or_nested_fusion"
            eej_ncf_rct = "RCT_like"
            review = "review"
            suggestion = "Two strong members can reflect fusion/nested fusion; check if each partner corresponds to a different ppa chromosome."
        else:
            event = "possible_break_or_split"
            eej_ncf_rct = "NCF_like"
            review = "review"
            suggestion = "Partial support with extra members; inspect for rearrangement junctions and ambiguous secondary hits."
    else:
        event = "complex_rearrangement"
        eej_ncf_rct = "complex"
        review = "review"
        suggestion = "Weak or fragmented support; manually inspect heatmap and ppa corroboration."

    if length_ratio >= 1.5 and strong_count >= 2:
        event = event + "; length_excess"
        suggestion += " Length substantially exceeds the reference chromosome; consider merged/redundant segments."
    elif 0.5 <= length_ratio < 0.8 and strong_count <= 3:
        event = event + "; length_short"
        suggestion += " Length is noticeably shorter than the reference; consider breakpoint, truncation, or missing sequence."

    if self_synteny_bp > 0 and member_count >= 4:
        event = event + "; self_synteny_support"

    return {
        "pcam_chr": pcam_chr,
        "pzi_members": ",".join(sorted(member_set, key=natural_key)),
        "strong_pzi_members": ",".join(sorted(strong_set, key=natural_key)),
        "member_count": str(member_count),
        "strong_member_count": str(strong_count),
        "self_synteny_bp": str(self_synteny_bp),
        "pcam_length_bp": str(pcam_len),
        "pzi_member_total_bp": str(pzi_total),
        "length_ratio_to_pcam": f"{length_ratio:.4f}",
        "event_class": event,
        "eej_ncf_rct_candidate": eej_ncf_rct,
        "review_status": review,
        "review_note": suggestion,
    }


def classify_ppa_links(ppa_summary, pzi_main_chr, allowed_ppa=None):
    by_ppa = defaultdict(list)
    for row in ppa_summary:
        if allowed_ppa is not None and row["query_chr"] not in allowed_ppa:
            continue
        if row["target_chr"] in pzi_main_chr:
            by_ppa[row["query_chr"]].append(row)
    rows = []
    for ppa_chr in sorted(by_ppa, key=natural_key):
        hits = sorted(by_ppa[ppa_chr], key=lambda r: (-float(r["coverage_ratio"] or 0), -float(r["aligned_bp"] or 0)))
        best = hits[0] if hits else None
        second = hits[1] if len(hits) > 1 else None
        if best is None:
            rows.append({
                "ppa_chr": ppa_chr,
                "best_pzi_chr": "",
                "second_pzi_chr": "",
                "best_coverage_ratio": "0.000000",
                "second_coverage_ratio": "0.000000",
                "classification": "unresolved",
                "review_note": "No confident best hit into the 32-chr assembly.",
            })
            continue
        best_cov = float(best["coverage_ratio"] or 0)
        second_cov = float(second["coverage_ratio"] or 0) if second else 0.0
        if best_cov >= 0.5 and (second is None or best_cov >= second_cov * 1.5):
            cls = "one_to_one"
            note = "PPA chromosome is predominantly aligned to a single PZI chromosome."
        elif best_cov >= 0.3 and second and second_cov >= 0.2:
            cls = "one_to_two_or_fusion_like"
            note = "PPA chromosome shows substantial support for two PZI chromosomes; check for fusion or split."
        else:
            cls = "complex_or_fragmented"
            note = "Weak or split support; inspect for rearrangements or incomplete PZI coverage."
        rows.append({
            "ppa_chr": ppa_chr,
            "best_pzi_chr": best["target_chr"],
            "second_pzi_chr": second["target_chr"] if second else "",
            "best_coverage_ratio": best["coverage_ratio"],
            "second_coverage_ratio": second["coverage_ratio"] if second else "0.000000",
            "classification": cls,
            "review_note": note,
        })
    return rows


def build_reverse_heatmap_summary(ppa_summary_rows):
    rows = []
    for row in ppa_summary_rows:
        if not row["query_chr"] or not row["target_chr"]:
            continue
        rows.append({
            "reference_chr": row["target_chr"],
            "query_chr": row["query_chr"],
            "aligned_block_bp": row["aligned_bp"],
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description="Trispecies rearrangement review for pzi/pcam/ppa.")
    parser.add_argument("--pzi-vs-pcam-paf", required=True, type=Path)
    parser.add_argument("--pzi-self-paf", required=True, type=Path)
    parser.add_argument("--ppa-vs-pzi-paf", required=True, type=Path)
    parser.add_argument("--pzi-fasta", required=True, type=Path)
    parser.add_argument("--pcam-fasta", required=True, type=Path)
    parser.add_argument("--ppa-fasta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    parser.add_argument("--ppa-main-count", type=int, default=16)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    pzi_rows = parse_paf(args.pzi_vs_pcam_paf, args.min_mapq, args.min_aln_len)
    pzi_self_rows = parse_paf(args.pzi_self_paf, args.min_mapq, args.min_aln_len)
    ppa_rows = parse_paf(args.ppa_vs_pzi_paf, args.min_mapq, args.min_aln_len)

    pzi_headers = [h.split()[0] for h in read_fasta_headers(args.pzi_fasta)]
    pcam_headers = [h.split()[0] for h in read_fasta_headers(args.pcam_fasta)]
    ppa_headers = [h.split()[0] for h in read_fasta_headers(args.ppa_fasta)]
    ppa_main_headers = [h for h in ppa_headers if h in PPA_MAIN_CHR][: args.ppa_main_count]

    pzi_summary, pzi_best = summarize_pairs(pzi_rows)
    write_tsv(args.out_dir / "pzi_vs_pcam_chr_alignment.tsv", pzi_summary, [
        "query_chr", "target_chr", "query_length_bp", "target_length_bp", "aligned_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_query"
    ])
    write_tsv(args.out_dir / "pzi_vs_pcam_best_hits.tsv", pzi_best, [
        "query_chr", "target_chr", "query_length_bp", "target_length_bp", "aligned_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_query"
    ])

    pcam_lengths = {row["target_chr"]: int(row["target_length_bp"]) for row in pzi_best if row["target_chr"]}
    pzi_lengths = {row["query_chr"]: int(row["query_length_bp"]) for row in pzi_best if row["query_chr"]}
    pzi_self_pairs = {(r["query"], r["target"]): int(r["aln_len"]) for r in pzi_self_rows}

    pcam_to_pzi = defaultdict(list)
    for row in pzi_best:
        if row["target_chr"]:
            pcam_to_pzi[row["target_chr"]].append(row["query_chr"])

    event_rows = []
    for pcam_chr in sorted(pcam_to_pzi, key=natural_key):
        event_rows.append(classify_pcam_event(pcam_chr, pcam_to_pzi[pcam_chr], pzi_best, pzi_self_pairs, pzi_lengths, pcam_lengths))

    write_tsv(args.out_dir / "pzi_pcam_rearrangement_review.tsv", event_rows, [
        "pcam_chr", "pzi_members", "strong_pzi_members", "member_count", "strong_member_count",
        "self_synteny_bp", "pcam_length_bp", "pzi_member_total_bp", "length_ratio_to_pcam",
        "event_class", "eej_ncf_rct_candidate", "review_status", "review_note"
    ])

    write_tsv(args.out_dir / "pzi_pcam_composition_summary.tsv", [
        {
            "pcam_chr": row["pcam_chr"],
            "member_count": row["member_count"],
            "strong_member_count": row["strong_member_count"],
            "event_class": row["event_class"],
            "eej_ncf_rct_candidate": row["eej_ncf_rct_candidate"],
            "review_status": row["review_status"],
            "review_note": row["review_note"],
        } for row in event_rows
    ], [
        "pcam_chr", "member_count", "strong_member_count", "event_class", "eej_ncf_rct_candidate", "review_status", "review_note"
    ])

    ppa_summary, ppa_best = summarize_pairs(ppa_rows)
    write_tsv(args.out_dir / "ppa_vs_pzi_chr_alignment.tsv", ppa_summary, [
        "query_chr", "target_chr", "query_length_bp", "target_length_bp", "aligned_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_query"
    ])
    write_tsv(args.out_dir / "ppa_vs_pzi_best_hits.tsv", ppa_best, [
        "query_chr", "target_chr", "query_length_bp", "target_length_bp", "aligned_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_query"
    ])

    pzi_main_chr = [row["query_chr"] for row in pzi_best if row["target_chr"]]
    ppa_class_rows = classify_ppa_links(ppa_summary, pzi_main_chr, allowed_ppa=set(ppa_main_headers))
    write_tsv(args.out_dir / "ppa_chr_review.tsv", ppa_class_rows, [
        "ppa_chr", "best_pzi_chr", "second_pzi_chr", "best_coverage_ratio", "second_coverage_ratio", "classification", "review_note"
    ])

    ppa_heatmap_rows = build_reverse_heatmap_summary(ppa_summary)
    write_tsv(args.out_dir / "ppa_vs_pzi_chr_alignment.tsv", ppa_summary, [
        "query_chr", "target_chr", "query_length_bp", "target_length_bp", "aligned_bp", "covered_bp",
        "coverage_ratio", "identity", "orientation", "best_hit_for_query"
    ])
    write_tsv(args.out_dir / "ppa_vs_pzi_heatmap_summary.tsv", ppa_heatmap_rows, [
        "reference_chr", "query_chr", "aligned_block_bp"
    ])

    axis_rows = []
    palette = ["#F6AD55", "#68D391", "#63B3ED", "#F687B3", "#4FD1C5", "#ED8936", "#9F7AEA", "#ECC94B", "#F56565", "#4299E1"]
    color_by_pcam = {}
    for idx, row in enumerate(event_rows):
        color = palette[idx % len(palette)]
        color_by_pcam[row["pcam_chr"]] = color
        axis_rows.append({"axis": "both", "chr": row["pcam_chr"], "color": color, "label": row["event_class"]})
        for member in row["pzi_members"].split(","):
            if member:
                axis_rows.append({"axis": "both", "chr": member, "color": color, "label": row["pcam_chr"]})
    for row in ppa_class_rows:
        color = color_by_pcam.get(row["best_pzi_chr"], "#CBD5E0")
        if row["ppa_chr"] in ppa_main_headers:
            axis_rows.append({"axis": "both", "chr": row["ppa_chr"], "color": color, "label": row["classification"]})
    write_tsv(args.out_dir / "axis_colors.tsv", axis_rows, ["axis", "chr", "color", "label"])

    review_rows = []
    for row in event_rows:
        review_rows.append({
            **row,
            "review_note": row["review_note"] + " Inspect dotplot and heatmap color strips for breakpoint, fusion, or excess-span signals.",
        })
    write_tsv(args.out_dir / "review_summary.tsv", review_rows, [
        "pcam_chr", "pzi_members", "strong_pzi_members", "member_count", "strong_member_count",
        "self_synteny_bp", "pcam_length_bp", "pzi_member_total_bp", "length_ratio_to_pcam",
        "event_class", "eej_ncf_rct_candidate", "review_status", "review_note"
    ])

    (args.out_dir / "README_TRI_SPECIES_REARRANGEMENT.md").write_text(
        "# Trispecies rearrangement review\n\n"
        "This folder contains chromosome-level review tables for P. zippeliana, P. campanulata, and P. padus.\n\n"
        "The outputs are provisional review aids only.\n"
        "Axis bands use colors to mark related chromosomes in the plots.\n"
        "Event classes are heuristic review labels, not final evolutionary calls.\n\n"
        "Provisional labels:\n"
        "- balanced_quartet\n"
        "- possible_break_or_split\n"
        "- possible_fusion_or_nested_fusion\n"
        "- complex_rearrangement\n\n"
        "EEJ/NCF/RCT labels are only candidate annotations in the review table.\n",
        encoding="utf-8",
    )

    print(f"Wrote tri-species rearrangement review tables into {args.out_dir}")


if __name__ == "__main__":
    main()
